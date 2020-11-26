"""Implements joining of blast hits into subject-accession bins"""

import copy
import logging
import math
from collections import defaultdict
from typing import NamedTuple, List, Set, Iterator, Optional, Tuple, Iterable, Dict, Any, Mapping

import tqdm  # type: ignore

from virprof.regionlist import RegionList

LOG = logging.getLogger(__name__)

# Values for Lambda and K used by BLAST to calculate evalues.
#
# These were empirically determined by Altschul et al and are
# hard-coded into BLAST.
BLAST_CONSTANTS = {
    # 'megablast': (1., 0.41),... unclear
    # word size 28, gap open 0, gap_extend 0, score 1/-2
    'blastn': (0.625, 0.41),
    # word size 11, gap open 5, gap_extend 2, score 2/-3
}

def calc_log10_evalue(score: float, query_len: int,
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
    blast_l, blast_k = BLAST_CONSTANTS[blast]
    return round((- blast_l * score
                  + math.log(blast_k)
                  + math.log(database_len)
                  + math.log(query_len))
                 /
                 math.log(10),
                 2)


def bitscore(score: float, blast: str = "blastn") -> float:
    """Compute BLAST bitscore

    Args:
        score: BLAST score
        blast: blast setting (only ``blastn`` for now)
    """
    blast_l, blast_k = BLAST_CONSTANTS[blast]
    return - round((- blast_l * score
                    + math.log(blast_k))
                   /
                   math.log(2),
                   2)


class BlastHit(NamedTuple):
    # pylint: disable=too-few-public-methods
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
    staxids: List[int]
    bitscore: int


class OverlapException(Exception):
    """Raised if trying to append hit overlapping with
    hits already contained within the hitchain.
    """


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
    __slots__ = ("hits", "chain_penalty", "_hash", "_score", "_subject_regions")

    def __init__(self,
                 hits: Optional[List[BlastHit]] = None,
                 chain_penalty: int = 20) -> None:
        #: List of hits sorted by subject start
        self.hits: List[BlastHit] = []
        self.chain_penalty = chain_penalty
        self._hash = None
        self._score = None
        self._subject_regions = RegionList()

        if hits:
            for hit in hits:
                self.add(hit)

    def _reset(self) -> None:
        """Reset cached score and length values"""
        self._score = None
        self._hash = None

    def __len__(self) -> int:
        return len(self.hits)

    def __str__(self) -> str:
        return (
            f"acc={self.sacc} e={self.log10_evalue} "
            f"l={self.qlen} r={self.sranges()}"
        )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({str(self)})"

    def _copy_from_other(self, other):
        self.hits = copy.copy(other.hits)
        self.chain_penalty = other.chain_penalty
        self._hash = other._hash
        self._score = other._score
        self._subject_regions = copy.copy(other._subject_regions)

    def __copy__(self) -> "HitChain":
        cpy = object.__new__(type(self))
        cpy._copy_from_other(self)
        return cpy

    def __hash__(self) -> int:
        if self._hash is None:
            val = self._hash = hash(tuple(self.hits))
            return val
        return self._hash

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, HitChain):
            return NotImplemented
        return self.hits == other.hits

    def add(self, hit: BlastHit) -> None:
        """Add a BLAST hit to the chain

        Args:
          hit: A BLAST hit
        """
        self.hits.append(hit)
        self._subject_regions.add(hit.sstart, hit.send, hit)
        self._reset()

    def prune(self, accs: Set[str]) -> None:
        """Remove BLAST hits from chain

        Removes hits from chain having qacc matching any of ``hits``.

        Args:
           hits: List of hits
        """
        oldlen = len(self.hits)
        oldhits = self.hits
        self.hits = [
            hit for hit in oldhits
            if hit.qacc not in accs
        ]
        if oldlen != len(self.hits):
            if not self.hits:
                self._subject_regions = RegionList()
            else:
                for hit in oldhits:
                    if hit.qacc in accs:
                        self._subject_regions.remove(hit.sstart, hit.send, hit)
            self._reset()

    def trunc(self, sstart: int) -> None:
        """Truncate HitChain at ``sstart``

        Args:
            sstart: Position in the reference sequence

        Returns: A HitChain comprising the longest prefix of hits from
            this HitChain that does not extend beyond the parameter
            ``start``.
        """
        hits = self.hits
        oldlen = len(hits)
        while hits and max(hits[-1].sstart, hits[-1].send) >= sstart:
            hits.pop()
        if oldlen != len(hits):
            self._reset()
            self.hits = hits

    @property
    def qlen(self) -> int:
        """Combined query length of all items in the HitChain.

        Note that ``len(hitChain)`` returns the number of items, not
        the length in basepairs.
        """
        seen = set()
        qlen = 0
        for hit in self.hits:
            if hit.qacc not in seen:
                seen.add(hit.qacc)
                qlen += hit.qlen
        return qlen

    @property
    def slen(self) -> int:
        """Combined non-redundant subject length of all items in the HitChain.

        Note that ``len(hitChain)`` returns the number of items, not
        the length in basepairs.
        """
        return sum(stop - start + 1 for start, stop, data in self._subject_regions if data)

    @property
    def blast_score(self) -> float:
        "Aggegated BLAST score for all items in the HitChain"
        return (
            sum(hit.score for hit in self.hits)
            - self.chain_penalty * (len(self) - 1)
        )

    @property
    def pident(self) -> Optional[float]:
        "Aggregated percent identity for all items in the HitChain"
        if not self.hits:
            return None
        ## FIXME: temporary hac
        matches = 0
        length = 0
        for start, stop, data in self._subject_regions:
            if data:
                matches += (stop - start + 1) * max(hit.pident for hit in data)
                length += (stop - start + 1)
        return round(matches/length, 1)

    @property
    def qpident(self) -> Optional[float]:
        "Aggregated query percent identity for all items in the HitChain"
        if not self.alen:
            return None

        matches = sum(hit.pident * hit.length for hit in self.hits)
        qpident = matches / self.qlen
        return round(qpident, 1)

    @property
    def pidents(self) -> List[float]:
        "Percent identiy for each hit"
        return [hit.pident for hit in self.hits]

    @property
    def log10_evalue(self) -> float:
        "Logarithm base 10 of the E-value computed for the combined hits"
        if self.qlen == 0:
            return 0
        return calc_log10_evalue(self.blast_score, self.qlen)

    @property
    def score(self) -> float:
        "Compute overall score for this hit chain"
        if self._score is None:
            val = self._score = -self.log10_evalue
            return val
        return self._score

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
        """List of NCBI Taxids for subject"""
        return self.hits[0].staxids if self.hits else []

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
        """Ranges of subject sequence covered by the hits in this HitChain

        Returns:
          List of pairs with start and end positions (start < end).
        """
        return [(hit.sstart, hit.send)
                if hit.sstart < hit.send
                else (hit.send, hit.sstart)
                for hit in self.hits]

    def range_reversed(self) -> List[bool]:
        return [hit.sstart > hit.send for hit in self.hits]

    @property
    def qaccs(self) -> List[str]:
        """The accession numbers for each hit in this HitChain

        This list may contain duplicates!
        """
        return [hit.qacc for hit in self.hits]

    def qranges(self, qacc: Optional[str] = None) -> List[Tuple[int, int]]:
        """Ranges of query sequences covered by the hits in this HitChain

        Args:
           qacc: Only get ranges for this query accession

        Returns:
           List of pairs with start and end positions (start < end).
        """
        # These are always qstart <= qend
        if qacc is None:
            return [(hit.qstart, hit.qend) for hit in self.hits]
        return [(hit.qstart, hit.qend) for hit in self.hits
                if hit.qacc == qacc]

    def make_chain_single(self, hitset: List[BlastHit]) -> List["HitChain"]:
        """Generates hit chain for set of hits with single sacc

        Uses the current HitChain as start point.
        """
        hitset.sort(key=lambda h: min(h.sstart, h.send))
        chain = copy.copy(self)
        for hit in hitset:
            chain.add(hit)
        return [chain]

    def make_chains(self, hits: Iterable[BlastHit], **args) -> List["HitChain"]:
        """Generates hit chains for set of hits

        Uses the current HitChain as start point
        """
        hitsets = defaultdict(list)
        for hit in hits:
            hitsets[hit.sacc].append(hit)
        return [
            chain
            for hitset in hitsets.values()
            for chain in self.make_chain_single(hitset, **args)
        ]

    def to_dict(self) -> Dict[str, Any]:
        """Convert class to dict for writing to CSV"""
        def range_str(ranges):
            return ";".join("{}-{}".format(*rang) for rang in ranges)

        return {
            'log_evalue': self.log10_evalue,
            'qlen': self.qlen,
            'slen': self.slen,
            'n_frag': len(self.hits),
            'sacc': self.sacc,
            'stitle': self.stitle,
            'staxids': ";".join((str(i) for i in self.staxids)),
            'pident': self.pident,
            'pidents': ";".join((str(pid) for pid in self.pidents)),
            'qaccs': ";".join(self.qaccs),
            'qranges': range_str(self.qranges()),
            'sranges': range_str(self.sranges()),
            'reversed': ";".join('T' if rev else 'F' for rev in self.range_reversed()),
        }

    @property
    def fields(self) -> List[str]:
        """List keys of ``to_dict()``"""
        return list(self.to_dict().keys())


class CheckOverlaps:
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
        qacc = hit.qacc
        for ohit in filter(lambda h: h.qacc == qacc, self.hits):
            ostart, oend = sorted((ohit.qstart, ohit.qend))
            if not oend < start and not end < ostart:
                yield ohit

    def prune_overlapping(self, hits: List[BlastHit]) -> None:
        """Remove overlapping BLAST hits from chain

        Removes hits from chain having qacc matching any of ``hits``
        and having qstart/qstop overlapping the hit.

        Args:
           hits: List of hits

        """
        oldlen = len(self.hits)
        overlapping = set(
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

    def make_chain_single(self, hitset: List[BlastHit],
                          check_overlap = True) -> List["HitChain"]:
        """Generates hit chains for set of hits with single sacc

        Uses the current HitChain as start point.
        """
        # start with empty chain
        chains = [copy.copy(self)]
        hitset.sort(key=lambda h: min(h.sstart, h.send))
        if len(hitset) > 50:
            hitset = tqdm.tqdm(hitset, desc=hitset[0].sacc)

        # iterate over hits in order of range on target
        for hit in hitset:
            newchains = set()
            sstart = min(hit.sstart, hit.send)
            # iterate over open chains
            for chain in chains:
                if chain.overlaps(hit):
                    newchain = copy.copy(chain)
                    newchain.trunc(sstart)
                    newchain.prune_overlapping([hit])
                    newchains.add(newchain)
                else:
                    chain.add(hit)
            for chain in newchains:
                chain.add(hit)
                chains.append(chain)
        return chains


class CoverageHitChain(HitChain):
    """HitChain also calculating read counts"""
    def __init__(self,
                 hits: Optional[List[BlastHit]] = None,
                 chain_penalty: int = 20) -> None:
        super().__init__(hits, chain_penalty)
        self._coverages: Mapping[str, Mapping[str, Mapping[str, str]]] = {}

    def set_coverage(self, coverages: Mapping[str, Mapping[str, Mapping[str, str]]]) -> None:
        """Set coverage data"""
        self._coverages = coverages
        self._units = list(coverages.keys())
        self._numreads = {
            qacc: [int(coverages[unit][qacc]['numreads'])
                   for unit in self._units]
            for qacc in coverages[self._units[0]]
        }

    def _copy_from_other(self, other):
        super()._copy_from_other(other)
        self._coverages = other._coverages
        self._units = other._units
        self._numreads = other._numreads

    @property
    def numreads(self) -> int:
        """Number of reads"""
        if self._coverages is None:
            return -1
        # Count unique qaccs only
        qaccs = set(self.qaccs)
        return sum(sum(self._numreads[qacc]) for qacc in qaccs)

    def get_numreads(self, qacc) -> List[int]:
        return self._numreads[qacc]

    @property
    def numreadss(self) -> str:
        """Number of reads"""
        if self._coverages is None:
            return ""
        return ";".join(str(sum(self.get_numreads(qacc))) for qacc in self.qaccs)

    def filter_hitgroups(
            self,
            hitgroups: Iterable[List[BlastHit]],
            min_read_count: int
    ) -> Iterable[List[BlastHit]]:
        filtered = 0
        for hitgroup in hitgroups:
            qacc = hitgroup[0].qacc
            if sum(self.get_numreads(qacc)) >= min_read_count:
                yield hitgroup
            else:
                filtered += 1
        LOG.info(f"Removed {filtered} contigs having < {min_read_count} reads")

    def to_dict(self) -> Dict[str, Any]:
        res = super().to_dict()
        res.update({
            'numreads': self.numreads,
            'numreadss': self.numreadss
        })
        return res


def greedy_select_chains(chains: List["HitChain"],
                         alt: float = 0.9) -> Iterator[List["HitChain"]]:
    """Greedily selects hit chains with best e-value score

    Iteratively finds and returns the best scoring chain together with
    chains with similarily good results. After each iteration, query
    sequences used are removed from the candidate hit chains.

    Similarly good chains are those that are within ``alt`` of the
    best scoring chain on evalue, alen and qlen and comprise at least
    90% of the same query sequences.

    """
    while chains:
        best_chain_idx = max(range(len(chains)),
                             key=lambda i: chains[i].score)
        best_chain = chains.pop(best_chain_idx)
        qaccs = set(best_chain.qaccs)

        # find similarly soring chains
        extra_chains = [
            chain for chain in chains
            if chain != best_chain
            and chain.score > alt * best_chain.score
            and chain.log10_evalue < alt * best_chain.log10_evalue
            and chain.slen > alt * best_chain.slen
            and chain.qlen > alt * best_chain.qlen
            and (len(qaccs.intersection(set(chain.qaccs))) > len(qaccs)/10)
        ]
        extra_chains.sort(key=lambda res: res.log10_evalue)

        yield [best_chain] + extra_chains

        # remove used query accs from remaining chains
        qaccs = set(best_chain.qaccs)
        for chain in chains:
            chain.prune(qaccs)
        # remove chains now empty
        chains = [chain for chain in chains if chain]
