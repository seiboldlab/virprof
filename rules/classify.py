#!/usr/bin/env python

"""
Given: a set of blast hits, matching contigs to virus DB
Sought: a set of viruses explaining the contigs

Each hit has
 - an overlap region (sstart, send, qstart, qend)
 - identity, score, bitscore, evalue

The evalue is the prop of observing match randomly
Bitscore = (lambda * score - ln K) / ln 2
Evalue = |query| * |db| * 2 ^ -bitscore
       = K * |query| * |db| * e^(-lambda*score)

Blastn has 2 / -3 mismatch and 5/2 gap open/extend
That means we have
Lambda, K,    H,     Alpha, Beta, Theta
0.625,  0.41, 0.78,  0.8,   -2,   99
(empircal data from Altschule et al)
"""



# Aggregating hits
# We have
# - n hits
# - each matching part of query to part of subject with score
# We want
# - Classifications for each sample
# - Parsimony: Minimum number of database sequences containing all contigs
# - ML: Set of database sequences under which contigs most likely
#   - i.e. contigs least likely to be random, i.e. best evalue total
#   - merge hits with same acc/contig when non-overlapping
#   - recalc evalue
# - for each sample

import csv
import os
import sys
import string
import math
from collections import defaultdict
import itertools
import logging
from abc import ABC, abstractmethod
from copy import copy
from typing import NamedTuple, List, Set, Iterable, Tuple, Optional, Iterator, Dict

import click
import tqdm # type: ignore
import ymp.blast  # type: ignore
import networkx as nx
import graph_tool as gt
import graph_tool.topology
import graph_tool.search
import graph_tool.util

class BlastHit(NamedTuple):
    """Base type for a BLAST hit

    This class is only used for type checking.
    """
    qacc: str
    score: float
    sacc: str
    send: int
    sstart: int
    qstart: int
    qend: int
    qlen: int
    pident: float
    length: int
    stitle: str
    staxids: int


class HitChain:
    """Chain of BLAST hits

    This class stores a sequence of BLAST hits from different query
    sequences to the same target sequence. It also provides methods
    to generate such chains from a collection of BLAST hits and to
    select non-overlapping sets.

    Args:
        hits: collection of hits

        chain_penalty: Penalty added to score for each item in the
            list beyond the first.
    """
    # Values for Lambda and K used by BLAST to calculate evalues.
    #
    # These were empirically determined by Altschul et al and are
    # hard-coded into blast.
    BLAST_CONSTANTS = {
        #'megablast': (1., 0.41),... unclear
        # word size 28, gap open 0, gap_extend 0, score 1/-2
        'blastn': (0.625, 0.41),
        # word size 11, gap open 5, gap_extend 2, score 2/-3
    }

    class OverlapException(Exception):
        """Raised if trying to append hit overlapping with
        hits already contained within the hitchain.
        """

    @classmethod
    def calc_log10_evalue(cls, score: float, query_len: int,
                          database_len: int = 1004024614,
                          blast: str = "blastn") -> float:
        """Compute BLAST evalue

        Args:
          score: BLAST score
          query_len: Length of query sequence
          database_len: Length of queried database
          blast: blast setting (only ``blastn`` for now)

        Returns:
          Exponent of E-Value
        """
        blast_l, blast_k = cls.BLAST_CONSTANTS[blast]
        return round((-blast_l * score
                      + math.log(blast_k)
                      + math.log(database_len)
                      + math.log(query_len)) / math.log(10),
                     2)

    @classmethod
    def bitscore(cls, score: float, blast: str = "blastn"):
        """Compute BLAST bitscore

        Args:
          score: BLAST score
          blast: blast setting (only ``blastn`` for now)
        """
        blast_l, blast_k = cls.BLAST_CONSTANTS[blast]
        return (blast_l * score - math.log(blast_k)) / math.log(2)

    def __init__(self,
                 hits: List[BlastHit] = None,
                 chain_penalty: int = 20) -> None:
        self.hits = copy(hits) if hits else list()
        self.chain_penalty = chain_penalty

        self._score: Optional[float] = None
        self._qlen: Optional[int] = None
        self._alen: Optional[int] = None
        self._pident: Optional[float] = None
        self._hash: Optional[int] = None

    def _reset(self) -> None:
        """Reset cached score and length values"""
        self._score = None
        self._qlen = None
        self._alen = None
        self._pident = None
        self._hash = None

    def __len__(self) -> int:
        return len(self.hits)

    def __str__(self) -> str:
        return f"acc={self.sacc} e={self.log10_evalue} l={self.qlen} r={self.sranges_str()}"

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}({str(self)})")

    def __copy__(self) -> "HitChain":
        return self.__class__(self.hits, self.chain_penalty)

    def __hash__(self) -> int:
        if self._hash is None:
            self._hash = hash(hit for hit in self.hits)
        return self._hash

    def __eq__(self, other) -> bool:
        return all(hit == ohit for hit, ohit in zip(self.hits, other.hits))

    def overlaps(self, hit: BlastHit) -> bool:
        """Checks for overlap with ``hit``

        Args:
           hit: The blast hit to check for overlaps
        Returns:
           True if ``hit`` overlaps on query or subject side
        """
        return self.overlaps_subject(hit) or self.overlaps_query(hit)

    def overlaps_subject(self, hit: BlastHit) -> bool:
        """Checks for subject overlap with ``hit``

        This assumes that the HitChain is filled left to
        right. No checks for gaps in the middle are done.

        Args:
          hit: The blast hit to check for overlaps
        Returns:
          True if ``hit`` overlaps on subject side
        """
        return self.send > min(hit.sstart, hit.send)

    def overlaps_query(self, hit: BlastHit) -> bool:
        """Checks for query overlap with ``hit``.

        This will compare ``hit`` to all BlastHits comprising this
        HitChain. If hit has the same query accession and overlaps
        the hit region on the query range on any hit, ``True`` is
        returned.

        Args:
          hit: The blast hit to check for overlaps
        Returns:
          True if ``hit`` overlaps on query side
        """
        return next(self.get_query_overlaps(hit), None) is not None

    def get_query_overlaps(self, hit: BlastHit) -> Iterator[BlastHit]:
        """Finds hits overlapping ``hit`` on query side

        Args:
          hit: The blast hit to check for overlaps
        Returns:
          An iterator over blast hits with query accession matching
          ``hit`` and overlapping query side of ``hit``.
        """
        start, end = sorted((hit.qstart, hit.qend))
        for ohit in self.hits:
            if hit.qacc == ohit.qacc:
                ostart, oend = sorted((ohit.qstart, ohit.qend))
                if not oend < start and not end < ostart:
                    yield ohit

    def append(self, hit: BlastHit) -> None:
        """Append a BLAST hit to the chain

        Args:
          hit: A BLAST hit
        """
        if self.overlaps_query(hit):
            raise self.OverlapException()
        self.hits.append(hit)
        self._reset()

    def prune(self, hits: List[BlastHit], full_sequence=False) -> None:
        """Remove BLAST hits from chain

        If ``full_sequence`` is ``True``, hits are removed if
        matching query accession number, otherwise, overlapping
        hits are removed.

        Args:
           hits: List of hits
           full_sequence: Whether to remove based on acc only
        """
        oldlen = len(self.hits)
        if full_sequence:
            accs = set(q.qacc for q in hits)
            self.hits = [
                hit for hit in self.hits
                if hit.qacc not in accs
            ]
        else:
            ohits = set(
                ohit
                for hit in hits
                for ohit in self.get_query_overlaps(hit)
            )
            self.hits = [
                ohit for ohit in self.hits
                if ohit not in ohits
            ]
        if oldlen != len(self.hits):
            self._reset()

    def trunc(self, sstart: int) -> "HitChain":
        """Truncate HitChain at ``sstart``

        Args:
            sstart: Position in the reference sequence

        Returns: A HitChain comprising the longest prefix of hits from
            this HitChain that does not extend beyond the parameter
            ``start``.
        """
        oldlen = len(self.hits)
        while hits and max(hits[-1].sstart, hits[-1].send) >= sstart:
            hits.pop()
        if oldlen != len(self.hits):
            self._reset()

    @property
    def qlen(self) -> int:
        """Combined query length of all items in the HitChain.

        Note that ``len(hitChain)`` returns the number of items, not
        the length in basepairs.
        """
        if self._qlen is None:
            seen = set()
            qlen = 0
            for hit in self.hits:
                if hit.qacc not in seen:
                    seen.add(hit.qacc)
                    qlen += hit.qlen
            self._qlen = qlen  # placate mypy
        return self._qlen

    @property
    def alen(self) -> int:
        """Combined alignment length of all items in the HitChain.

        Note that ``len(hitChain)`` returns the number of items, not
        the length in basepairs.
        """
        if self._alen is None:
            alen = sum(hit.length for hit in self.hits)
            self._alen = alen  # placate mypy
        return self._alen

    @property
    def score(self) -> float:
        "Aggegated BLAST score for all items in the HitChain"
        if self._score is None:
            self._score = (
                sum(hit.score for hit in self.hits)
                - self.chain_penalty * (len(self) - 1)
            )
        return self._score

    @property
    def pident(self) -> float:
        "Aggregated percent identity for all items in the HitChain"
        if self._pident is None and self.alen:
            matches = sum(hit.pident * hit.length for hit in self.hits)
            pident = matches / self.alen
            self._pident = round(pident,1)
        return self._pident

    @property
    def log10_evalue(self) -> float:
        "Logarithm base 10 of the E-value computed for the combined hits"
        if self.qlen == 0:
            return None
        return self.calc_log10_evalue(self.score, self.qlen)

    @property
    def sacc(self) -> str:
        """The reference accession number (sacc BLAST format field)"""
        return self.hits[0].sacc if self.hits else ""

    @property
    def stitle(self) -> str:
        """The reference title (stitle BLAST format field)"""
        return self.hits[0].stitle if self.hits else ""

    @property
    def staxids(self) -> List[int]:
        return self.hits[0].staxids if self.hits else []

    def staxids_str(self) -> str:
        return ";".join((str(i) for i in self.staxids))

    @property
    def send(self) -> int:
        """Position of the last base covered by the hits in the HitChain

        If the hit chain is empty, the position will be 0.
        """
        if not self.hits:
            return 0
        hit = self.hits[-1]
        return max(hit.send, hit.sstart)

    def sranges(self) -> List[Tuple[int, int]]:
        """A list of ranges covered by the hits in this HitChain"""
        return [(hit.sstart, hit.send)
                if hit.sstart < hit.send
                else (hit.send, hit.sstart)
                for hit in self.hits]

    def sranges_str(self) -> str:
        return ";".join("{}-{}".format(*tpl) for tpl in self.sranges())

    @property
    def qaccs(self) -> Set[str]:
        """The accession numbers for each hit in this HitChain"""
        return list(hit.qacc for hit in self.hits)

    def qranges(self, qacc: str = None) -> List[Tuple[int, int]]:
        """A list of ranges covered by the hits in this HitChain"""
        return [
            (hit.qstart, hit.qend)
            if hit.qstart < hit.qend
            else (hit.qend, hit.qstart)
            for hit in self.hits
            if qacc is None or hit.qacc == qacc
        ]

    def qranges_str(self) -> str:
        return ";".join("{}-{}".format(*tpl) for tpl in self.qranges())

    def _make_chains(self, hitset: List[BlastHit]) -> List["HitChain"]:
        """Generates hit chains for set of hits with single sacc

        Uses the current HitChain as start point.
        """
        # start with empty chain
        chains = [copy(self)]

        # iterate over hits in order of range on target
        for hit in sorted(hitset, key=lambda h: min(h.sstart, h.send)):
            newchains = set()
            sstart = min(hit.sstart, hit.send)
            # iterate over open chains
            for chain in chains:
                try:
                    chain.append(hit)
                except self.OverlapException:
                    chain = copy(chain)
                    chain.prune([hit])
                    newchains.add(chain)
            for chain in newchains:
                chain.append(hit)
                chains.append(chain)
        return chains

    def make_chains(self, hits: Iterable[BlastHit]) -> Iterator["HitChain"]:
        """Generates hit chains for set of hits

        Uses the current HitChain as start point
        """
        hitsets = defaultdict(list)
        for hit in hits:
            hitsets[hit.sacc].append(hit)
        logging.error("Got %i hitsets", len(hitsets))
        return itertools.chain.from_iterable(
            self._make_chains(hits)
            for hits in hitsets.values()
        )

    @classmethod
    def greedy_select_chains(cls, it: Iterable["HitChain"], alt: float=0.9
    ) -> Iterator["HitChain"]:
        """Greedily selects hit chains with best e-value score"""
        chains = list(it)
        while chains:
            best_chain_id = min(range(len(chains)),
                                key=lambda i: chains[i].log10_evalue)
            result = chains.pop(best_chain_id)
            altresult = sorted((
                chain for chain in chains
                if chain != result
                and chain.log10_evalue < alt * result.log10_evalue
                and chain.alen > alt * result.alen
                and chain.qlen > alt * result.qlen
            ), key=lambda res: res.log10_evalue)

            yield result, altresult

            for chain in chains:
                chain.prune(result.hits, full_sequence=True)
            chains = [chain for chain in chains if len(chain) > 0]

    def to_dict(self):
        """Convert class to dict for writing to CSV"""
        return {
            'sacc': self.sacc,
            'stitle': self.stitle,
            'staxids': self.staxids_str(),
            'log_evalue': self.log10_evalue,
            'pident': self.pident,
            'qlen': self.qlen,
            'alen': self.alen,
            'n_frag': len(self.hits),
            'qaccs': ";".join(self.qaccs),
            'qranges': self.qranges_str(),
            'sranges': self.sranges_str(),
        }

    @property
    def fields(self):
        """List keys of ``to_dict()``"""
        return list(self.to_dict().keys())


class CoverageHitChain(HitChain):
    """HitChain also calculating coverage(s)"""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._coverages = None

    def set_coverage(self, coverages):
        self._coverages = coverages

    def __copy__(self) -> "CoverageHitChain":
        cpy = super().__copy__()
        cpy.set_coverage(self._coverages)
        return cpy

    @property
    def numreads(self):
        """Number of reads"""
        if self._coverages is None:
            return -1
        return sum(
            int(cov[qacc]['numreads'])
            for qacc in self.qaccs
            for cov in self._coverages.values()
        )

    def to_dict(self):
        res = super().to_dict()
        res.update({
            'numreads': self.numreads
        })
        return res


class Taxonomy(ABC):
    """Base class holding NCBI taxonomy

    The tree is implemented in subclasses TaxonomyGT and
    TaxonomyNX, using graph-tool and networkx libraries,
    respectively.

    Args:
      path: Path to the place where names.dmp and nodes.dmp
            can be found.
    """
    def __init__(self, path: str) -> None:
        self.path = path
        self.names_fn = os.path.join(self.path, "names.dmp")
        self.nodes_fn = os.path.join(self.path, "nodes.dmp")

        try:
            self.tree = self.load_tree_binary()
        except FileNotFoundError:
            self.tree = self.load_tree()
            self.save_tree_binary()

    @abstractmethod
    def load_tree_binary(self):
        """Load the tree from implementation specific binary"""
        raise NotImplementedError()

    @abstractmethod
    def save_tree_binary(self):
        """Save the tree to implementation specific binary"""
        raise NotImplementedError()

    @abstractmethod
    def load_tree(self):
        """Load the tree from NCBI dmp files"""
        raise NotImplementedError()

    @property
    def root(self) -> int:
        """The ``tax_id`` of the root node

        Returns:
          Always 1
        """
        return 1

    def node_reader(self):
        """Reader for NCBI names.dmp listing nodes and properties

        Returns:
           Generator of tax_id and property dictionary. The
           latter has ``name`` set to the scientific name from
           the NCBI file.
        """
        with open(self.names_fn, "r") as fd:
            reader = csv.DictReader(
                fd, delimiter='|',
                fieldnames=[
                    'tax_id', 'name', 'unique_name', 'name_class'
                ])
            for row in reader:
                if row['name_class'].strip() == 'scientific name':
                    tax_id = int(row['tax_id'])
                    data = dict((('name', row['name'].strip()),))
                    yield tax_id, data

    def edge_reader(self):
        """Reader for NCBI nodes.dmp, listing node and parent ids

        Returns:
           Generator over taxids ``(source, target)``
        """
        with open(self.nodes_fn, "r") as fd:
            reader = csv.DictReader(
                fd, delimiter='|',
                fieldnames = [
                    'tax_id', 'parent_id', 'rank'
                ]
            )
            for row in reader:
                source_node = int(row['parent_id'])
                target_node = int(row['tax_id'])
                yield source_node, target_node

    @abstractmethod
    def get_name(self, tax_id: int) -> str:
        """Get the scientific name of a node

        Args:
          tax_id: ID of taxonomy node
        Returns:
          The NCBI taxonomy scientific name, or "Unknown"
          for ``tax_id``s not found or named in the tree
        """
        raise NotImplementedError()

    @abstractmethod
    def get_lineage(self, tax_id: int) -> str:
        """Get the lineage for a node

        Args:
          tax_id: ID of taxonomy node

        Returns:
          The names of all nodes from (excluding) the root node until
          (including) node given in ``tax_id``. Names are separated by
          semicolons. If the node ``tax_id`` is not found in the tree,
          "Unknown" is returned. Nodes not described are not mentioned
          in the linage string.
        """
        raise NotImplementedError()

    @abstractmethod
    def get_subtree_ids(self, name: str) -> set():
        """Get the ``tax_id``s for node ``name`` and all subnodes

        Args:
          name: Scientific name of node

        Returns:
          Set of ``tax_id``s of the named node and all subnodes.
        """
        raise NotImplementedError()

    def make_filter(self, include: List[str]=None,
                    exclude: List[str]=None):
        """Create a filtering function

        Args:
          include: List of strings matching NCBI taxonomy scientific
            names of nodes. Only tax_ids of nodes below any of the
            listed "clades" will be accepted by the generated filter
            function.

          exclude: List of strings matching NCBI taxonomy scientific
            names of nodes. IDs of nodes below any of the listed
            "clades" will be rejected by the generated filter
            function.

        Returns:
          A function that takes an integer ``tax_id`` and returns
          a boolean indicating whether the tax_id should be included.
        """
        if include is None:
            include_ids = set()
        else:
            include_ids = set(
                tid
                for name in include
                for tid in self.get_subtree_ids(name)
            )
        if exclude is None:
            exclude_ids = set()
        else:
            exclude_ids = set(
                tid
                for name in exclude
                for tid in self.get_subtree_ids(name)
            )
        if not include_ids:
            if not exclude_ids:
                def nullFilter(taxid: int) -> bool:
                    return True
                return nullFilter
            else:
                def excludeFilter(taxid: int) -> bool:
                    return taxid not in exclude_ids
                logging.info("Made filter excluding %i tax_ids",
                             len(exclude_ids))
                return excludeFilter
        else:
            include_ids = include_ids.intersection(exclude_ids)
            def includeFilter(taxid: int) -> bool:
                return taxid in include_ids
            logging.info("Made filter including %i tax_ids",
                         len(include_ids))

            return includeFilter

    def isNull(self):
        return False


class TaxonomyGT(Taxonomy):
    """NCBI Taxonomy class using ``graph-tool`` to manage tree"""
    def load_tree(self) -> gt.Graph:
        tree = gt.load_graph_from_csv(
            self.nodes_fn,
            directed=True, eprop_types=[], string_vals=False,
            ecols=(1,0),
            csv_options={'delimiter':'|'})
        tree.vp.name = tree.new_vertex_property("string")
        for tax_id, data in self.node_reader():
            tree.vp.name[tree.vertex(tax_id)] = data['name']
        return tree

    def load_tree_binary(self) -> gt.Graph:
        return gt.load_graph('taxtree.gt')

    def save_tree_binary(self) -> None:
        self.tree.save('taxtree.gt')

    def get_name(self, taxid: int) -> str:
        try:
            return self.tree.vp.name[self.tree.vertex(taxid)] or 'Unknown'
        except ValueError:
            return 'Unknown'

    def get_lineage(self, taxid: int) -> str:
        try:
            target = self.tree.vertex(taxid)
        except ValueError:
            return 'Unknown'
        nodes, _ = gt.topology.shortest_path(self.tree, self.root, target)
        return '; '.join(self.tree.vp.name[node]
                         for node in nodes[1:])

    def get_subtree_ids(self, name):
        vertices = gt.util.find_vertex(self.tree, self.tree.vp.name, name)
        if len(vertices) != 1:
            return set()
        ids = set(
            int(tgt) for _, tgt in
            gt.search.dfs_iterator(self.tree, vertices[0])
        )
        ids.add(int(vertices[0]))
        return ids


class TaxonomyNX(Taxonomy):
    def load_tree_binary(self) -> nx.DiGraph:
        return nx.read_gpickle('taxtree.pkl')

    def save_tree_binary(self):
        nx.write_gpickle(self.tree, 'taxtree.pkl')

    def load_tree(self):
        tree = nx.DiGraph()
        tree.add_nodes_from(self.node_reader())
        tree.add_edges_from(self.edge_reader())
        return tree

    def get_subtree_ids(self, name):
        for node, data in self.tree.nodes(data=True):
            if data.get('name') == name:
                taxid = node
                break
        else:
            return []

        return [node for node in nx.descendants(self.tree, taxid)] + [taxid]

    def get_lineage(self, taxid):
        if taxid not in self.tree:
            return "Unknown"
        path = nx.shortest_path(self.tree, self.root, taxid)
        return '; '.join(self.tree.nodes[node].get('name')
                         for node in path[1:])

    def get_name(self, taxid):
        if taxid not in self.tree:
            return "Unknown"
        return self.tree.nodes[taxid].get('name')


class TaxonomyNull(Taxonomy):
    def load_tree_binary(self):
        return None

    def save_tree_binary(self):
        pass

    def load_tree(self):
        pass

    def get_name(self, tax_id: int) -> str:
        return "Unknown"

    def get_lineage(self, tax_id: int) -> str:
        return "Unknown"

    def get_subtree_ids(self, name: str) -> Set[int]:
        return set()

    def make_filter(self, include, exclude):
        def nullFilter(taxid: int) -> bool:
            return True
        return nullFilter

    def isNull(self):
        return True


def load_taxonomy(path: str):
    """Makes object of class Taxonomy

    If ``path`` is ``None``, returns an object of
    class `TaxonomyNull`, otherwise instanciates
    `TaxonomyGT`
    """
    if path is None:
        return TaxonomyNull(path)
    return TaxonomyGT(path)


class WordScorer:
    """Creates best-scoring-word summary from sequence titles

    This class takes a set of `HitChain`s and generates a string
    from the ``stitle`` fields containing the sequence description
    using word order, frequency and HitChain score.

    Args:
      stopwords: List of words to omit from output
      maxwordlen: Maximum length of words to include in output
      keepwords: Maximum number of words to include in output
    """

    def __init__(self,
                 stopwords = set([
                     'genome', 'complete', 'sequence', 'strain',
                     'isolate', 'and', 'or', 'sapiens', 'genes',
                     'predicted', 'human', 'homo', 'assembly',
                     'mus', 'musculus', 'virus', 'polyprotein'
                 ]),
                 maxwordlen = 27,
                 orderweight = 2,
                 keepwords = 4):
        self.stopwords = stopwords
        self.maxwordlen = maxwordlen
        self.orderweight = orderweight
        self.keepwords = keepwords(words)

    @staticmethod
    def split_words(chain: HitChain) -> Iterable[str]:
        """Split title of ``chain`` into words"""
        return chain.stitle.split()

    @staticmethod
    def strip_punctuation(words: Iterable[str]) -> Iterable[str]:
        """Remove punctuation from the outside of each word

        Args:
          words: sequence of words
        Returns:
          Generator over words
        """
        for word in words:
            yield word.strip(string.punctuation)

    @staticmethod
    def uncapitalize(words: Iterable[str]) -> Iterable[str]:
        """Remove capitalization for words without internal caps

        Args:
          words: sequence of words
        Returns:
          Generator over words
        """
        for word in words:
            if word[1:].isLower():
                yield word.lower()
            else:
                yield word

    @staticmethod
    def combine_shortwords(words: Iterable[str]) -> Iterable[str]:
        """Combine short words with their predecessor

        Args:
          words: sequence of words
        Returns:
          Generator over words
        """
        stack = []
        for word in words:
            if len(word) <= 2:
                stack += [' '.join((stack[-1], word))]
            else:
                yield from reversed(stack)
                stack  = [word]
        yield from reversed(stack)

    @staticmethod
    def filter_stopwords(words: Iterable[str]) -> Iterable[str]:
        """Filter out words occurring in our stopword list

        Args:
          words: sequence of words
        Returns:
          Generator over words
        """
        for word in words:
            if word.lower() not in self.stopwords:
                yield word

    def score(self, chains: List[HitChain]) -> str:
        """Generate summary using word scoring

        Breaks the title of each passed in HitChain into words,
        removes outside punctuation, removes captitalization, merges
        short words with their preceding word(s) and removes
        (singleton) stopwords.

        Each word occurring first in any input title receives a score
        equal to the sum of the log-evalues (exponent) of the
        respective HitChain. The words and scores are added with each
        subsequent word, dividing the log-evalue by ``orderweight`` as
        we move through the title.

        The resulting words are sorted by score and emitted as
        result. Up to ``keepwords`` words may comprise the
        result. Words composed from multiple words are split for this
        accounting, and no (sub)word may occur twice.
        """
        for chain in chains:
            words = self.split_words(words)
            words = self.strip_punctuation(words)
            words = self.uncapitalize(words)
            words = self.combine_shortwords(words)
            words = self.filter_stopwords(words)

            for idx, word in enumerate(words):
                word_scores[word] += chain.log10_evalue / self.orderweight ** idx

        sorted_words = sorted(word_scores.items(), key=lambda k: k[1])

        result_words = []
        seen = set()
        for word, score in sorted_words:
            parts = [part for part in word.split()
                     if part not in seen]
            seen.update(parts)
            for part in parts:
                result_words.append(part)
            if len(result_words) > self.keepwords:
                break
        return " ".join(result_words[:self.keepwords])

    def score_taxids(self, chains: List[HitChain]) -> int:
        """Selects best scoring taxid from chains

        Each input HitChain may have multiple NCBI taxonomy IDs
        assigned (NCBI taxonomy may assign multiple to a single
        reference sequence). For each found ID, the respective
        log-evalues are summed and the best score selected.
        """
        taxid_scores = defaultdict(float)
        for chain in chains:
            for taxid in chain.staxids:
                taxid_scores[taxid] += chain.log10_evalue
        staxids = [
            taxid for taxid, score
            in sorted(taxid_scores.items(), key=lambda k: k[1])
        ]
        return staxids[0]


def setup_logging():
    """Sets up python logging facility"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S'
    )


def parse_file_args(files):
    """Composes sample dictionary from list of files passed on cmdline

    - If we have only files ending in ``.blast7``, each is a sample of
      its own.
    - If we have an equal number of ``.blast7`` and ``.coverage``
      files, each pair is a sample and the basenames must match.
    - If we have one file ending in ``.blast7`` and a number of
      ``.coverage`` files, we have one sample with multiple coverage
      files (e.g. from a co-assembly over technical replicates or time
      series).

    Returns:
      A dictionary with sample names (file base names) as key and
      pairs of file-descriptor and coverage-file dictionary as
      value. The coverage-file dictionary itself has the coverage file
      basename as keys and the respective file-descriptor as values. E.g.

        ```
           {'sample1': (sample_1_blast7_file, {'coverage1': coverage1_fd})}
        ```
    """
    # Sort file arguments by extension
    cov_files = {}
    blast7_files = {}
    other_files = []
    for fd in files:
        name = os.path.basename(fd.name)
        if name.endswith('.coverage'):
            cov_files[name[:-len('.coverage')]] = fd
        elif name.endswith('.blast7'):
            blast7_files[name[:-len('.blast7')]] = fd
        else:
            other_files.append(name)
    if other_files:
        raise click.BadArgumentUsage(
            f"Unknown file type in arguments: {other_files}"
        )

    # Build result dictionary
    samples = {}
    if cov_files:
        chain_class = CoverageHitChain
        if len(blast7_files) > 1:
            basenames = set(cov_files.keys()) | set(blast7_files.keys())
            missing_files = set(
                name for name in basenames
                if name not in cov_files
                or name not in blast7_files
            )
            if missing_files:
                raise click.BadArgumentUsage(
                    f"Missing either blast7 or coverage file for:"
                    f" {missing_files}"
                )
            for name in basenames:
                samples[name] = (blast7_files[name], {name: cov_files[name]})
            logging.info("Processing %i samples with coverage", len(samples))
        else:
            name = next(iter(blast7_files.keys()))
            samples[name] = (blast7_files[name], cov_files)
            logging.info("Processing single sample with %i coverage files",
                         len(cov_files))
    else:
        chain_class = HitChain
        for name in blast7_files:
            samples[name] = (blast7_files[name], {})
        logging.info("Processing %i samples without coverage", len(samples))
    return chain_class, samples


def group_hits_by_qacc(hits):
    """Groups input hits into lists with same query accession"""
    qacc = None
    group = []
    for hit in hits:
        if hit.qacc != qacc:
            if group:
                yield group
            group = []
            qacc = hit.qacc
        group.append(hit)
    if group:
        yield group


@click.command()
@click.argument('files', type=click.File('r', lazy=True), nargs=-1)
@click.option('--out', '-o', type=click.File('w'), default='-',
              help="Output CSV file")
@click.option('--ncbi-taxonomy', '-t', type=click.Path(),
              help="Path to an unpacked NCBI taxonomy")
@click.option('--include', '-i', multiple=True,
              help="Scientific name of NCBI taxonomy node identifying"
              " subtree to include in output")
@click.option('--exclude', '-e', multiple=True,
              help="Scientific name of NCBI taxonomy node identifying"
              " subtree to exclude from output")
@click.option('--no-standard-excludes', '-E', type=bool,
              help="Do not exclude Human, artificial and unclassified"
              " sequences by default", default=False)
@click.option('--chain-penalty', type=int, default=20,
             help="Cost to BLAST score for concatenating two hits")
@click.option('--num-words', type=int, default=4,
              help="Number of words to add to 'words' field")
def main(files, out, ncbi_taxonomy=None, include=None,
         exclude=None, no_standard_excludes=False,
         chain_penalty=20, num_words=4):
    """Merge and classify contigs based on BLAST search results

    """
    setup_logging()

    if not files:
        logging.info("No files to process")
        return
    chain_class, samples = parse_file_args(files)
    chain_tpl = chain_class(chain_penalty=chain_penalty)

    wordscorer = WordScorer(keepwords=num_words)

    taxonomy = load_taxonomy(ncbi_taxonomy)

    if not no_standard_excludes:
        exclude += (
            'Homo sapiens',
            'Mus musculus',
            'artificial sequences',
            'unclassified sequences',
        )
        # cellular organisms
        # Microviridae
        # Caudovirales

    prefilter = taxonomy.make_filter(exclude=exclude)
    taxfilter = taxonomy.make_filter(include, exclude)

    field_list = ['sample']
    field_list += ['words']
    if not taxonomy.isNull():
        field_list += ['taxname']
    field_list += chain_tpl.fields
    field_list += ['taxid']
    if not taxonomy.isNull():
        field_list += ['lineage']
    logging.info("Output fields: %s", ' '.join(field_list))

    writer = csv.DictWriter(out, field_list)
    writer.writeheader()
    logging.info("Writing to: %s", out.name)

    for sample, (fd, cov_files) in tqdm.tqdm(samples.items()):
        logging.info("Processing %s", sample)

        # Do prefiltering
        reader = group_hits_by_qacc(ymp.blast.reader(fd))
        n_queries = n_filtered = 0
        hits = []
        for hitgroup in reader:
            n_queries += 1
            minscore = max(hit.bitscore for hit in hitgroup) * 0.9
            top_keep = [
                all(prefilter(taxid) for taxid in hit.staxids)
                for hit in hitgroup
                if hit.bitscore > minscore
            ]
            keep_group = len(top_keep)/2 < top_keep.count(True)
            if keep_group:
                hits += hitgroup
            else:
                n_filtered += 1
        logging.info("  %i query sequences (contigs) had hits",
                     n_queries)
        logging.info("  %i matched excludes", n_filtered)
        logging.info("  %i HSPs to process", len(hits))

        # Load Coverage
        cov = {}
        for cov_sample, cov_fd in cov_files.items():
            logging.error("Loading coverage for %s %s", sample, cov_sample)
            cov_reader = csv.DictReader(cov_fd, delimiter="\t")
            cov_data = {row['#rname']:row for row in cov_reader}
            logging.error("Got %i rows", len(cov_data))
            cov[cov_sample] = cov_data
            cov_fd.close()
        if cov:
            chain_tpl.set_coverage(cov)

        # Generate Chains
        chains = chain_tpl.make_chains(hits)
        fd.close()  # not needed any more, keep open fd's low

        # Choose "best" chains
        for chain, altchains in chain_tpl.greedy_select_chains(chains):
            allchains = [chain] + (altchains or [])
            taxid = wordscorer.score_taxids(allchains)
            if not taxfilter(taxid):
                continue

            row = chain.to_dict()
            row['sample'] = sample
            row['taxid'] = taxid
            row['words'] = wordscorer.score(allchains)
            row['lineage'] = taxonomy.get_lineage(taxid)
            row['taxname'] = taxonomy.get_name(taxid)

            writer.writerow(row)


if __name__ == "__main__":
    main()
