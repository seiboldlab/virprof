"""Implements joining of blast hits into subject-accession bins"""

import copy
import logging
import math
from collections import defaultdict
from typing import NamedTuple, List, Iterator, Optional, Tuple, Iterable, Dict, Any, Mapping

import tqdm  # type: ignore

LOG = logging.getLogger(__name__)

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
    # hard-coded into BLAST.
    BLAST_CONSTANTS = {
        # 'megablast': (1., 0.41),... unclear
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
        return round(
            (- blast_l * score
             + math.log(blast_k)
             + math.log(database_len)
             + math.log(query_len))
            /
            math.log(10),
            2)

    @classmethod
    def bitscore(cls, score: float, blast: str = "blastn") -> float:
        """Compute BLAST bitscore

        Args:
          score: BLAST score
          blast: blast setting (only ``blastn`` for now)
        """
        blast_l, blast_k = cls.BLAST_CONSTANTS[blast]
        return - round(
            (- blast_l * score
             + math.log(blast_k))
            /
            math.log(2),
            2)

    def __init__(self,
                 hits: Optional[List[BlastHit]] = None,
                 chain_penalty: int = 20) -> None:
        self.hits = hits.copy() if hits else list()
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
        return (
            f"acc={self.sacc} e={self.log10_evalue} "
            f"l={self.qlen} r={self.sranges_str()}"
        )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({str(self)})"

    def __copy__(self) -> "HitChain":
        return self.__class__(self.hits, self.chain_penalty)

    def __hash__(self) -> int:
        if self._hash is None:
            self._hash = hash(tuple(self.hits))
        return self._hash

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, HitChain):
            return NotImplemented
        return self.hits == other.hits

    def overlaps(self, hit: BlastHit) -> bool:
        """Checks for overlap with ``hit``

        Args:
           hit: The blast hit to check for overlaps
        Returns:
           True if ``hit`` overlaps on query or subject side
        """
        return False #self.overlaps_subject(hit) # or self.overlaps_query(hit)

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

    def append(self, hit: BlastHit) -> None:
        """Append a BLAST hit to the chain

        Args:
          hit: A BLAST hit
        """
        if self.overlaps(hit):
            raise self.OverlapException()
        self.hits.append(hit)
        self._reset()

    def prune(self,
              hits: List[BlastHit],
              full_sequence: bool = False) -> None:
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
    def pident(self) -> Optional[float]:
        "Aggregated percent identity for all items in the HitChain"
        if self._pident is None and self.alen:
            matches = sum(hit.pident * hit.length for hit in self.hits)
            pident = matches / self.alen
            self._pident = round(pident, 1)
        return self._pident

    @property
    def log10_evalue(self) -> float:
        "Logarithm base 10 of the E-value computed for the combined hits"
        if self.qlen == 0:
            return 0
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
        """List of NCBI Taxids for subject"""
        return self.hits[0].staxids if self.hits else []

    def staxids_str(self) -> str:
        """Semi-colon separated string listing NCBI taxids for subject"""
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
        """List of ranges covered by the hits in this HitChain"""
        return [(hit.sstart, hit.send)
                if hit.sstart < hit.send
                else (hit.send, hit.sstart)
                for hit in self.hits]

    def sranges_str(self) -> str:
        """Semicolon separated string of `sranges()`"""
        return ";".join("{}-{}".format(*tpl) for tpl in self.sranges())

    @property
    def qaccs(self) -> List[str]:
        """The accession numbers for each hit in this HitChain"""
        return list(hit.qacc for hit in self.hits)

    def qranges(self, qacc: Optional[str] = None) -> List[Tuple[int, int]]:
        """A list of ranges covered by the hits in this HitChain"""
        return [
            (hit.qstart, hit.qend)
            if hit.qstart < hit.qend
            else (hit.qend, hit.qstart)
            for hit in self.hits
            if qacc is None or hit.qacc == qacc
        ]

    def qranges_str(self) -> str:
        """Semicolon separated string of qranges"""
        return ";".join("{}-{}".format(*tpl) for tpl in self.qranges())

    def _make_chains(self, hitset: List[BlastHit]) -> List["HitChain"]:
        """Generates hit chains for set of hits with single sacc

        Uses the current HitChain as start point.
        """
        # start with empty chain
        chains = [copy.copy(self)]
        hitset.sort(key=lambda h: min(h.sstart, h.send))
        if len(hitset) > 10:
            hitset = tqdm.tqdm(hitset, desc=hitset[0].sacc)

        # iterate over hits in order of range on target
        for hit in hitset:
            newchains = set()
            sstart = min(hit.sstart, hit.send)
            # iterate over open chains
            for chain in chains:
                try:
                    chain.append(hit)
                except self.OverlapException:
                    newchain = copy.copy(chain)
                    newchain.trunc(sstart)
                    newchain.prune([hit])
                    newchains.add(newchain)
            for chain in newchains:
                chain.append(hit)
                chains.append(chain)
        return chains

    def make_chains(self, hits: Iterable[BlastHit]) -> List["HitChain"]:
        """Generates hit chains for set of hits

        Uses the current HitChain as start point
        """
        hitsets = defaultdict(list)
        for hit in hits:
            hitsets[hit.sacc].append(hit)
        LOG.info("Making chains from %i hitsets", len(hitsets))
        return [chain
                for hitset in tqdm.tqdm(hitsets.values(),
                                        desc="Making chains",
                                        disable=None)
                for chain in self._make_chains(hitset)]

    @classmethod
    def greedy_select_chains(
            cls, chains: List["HitChain"],
            alt: float = 0.9) -> Iterator[
                Tuple["HitChain", List["HitChain"]]]:
        """Greedily selects hit chains with best e-value score"""
        total = len(chains)
        pbar = tqdm.tqdm(total=total, desc="Pruning chains")
        while chains:
            best_chain_id = min(range(len(chains)),
                                key=lambda i: chains[i].log10_evalue)
            result = chains.pop(best_chain_id)
            qaccs = set(result.qaccs)
            altresult = sorted((
                chain for chain in chains
                if chain != result
                and chain.log10_evalue < alt * result.log10_evalue
                and chain.alen > alt * result.alen
                and chain.qlen > alt * result.qlen
                and (len(qaccs.intersection(set(chain.qaccs)))
                     > len(qaccs) * alt)
            ), key=lambda res: res.log10_evalue)

            yield result, altresult

            for chain in chains:
                chain.prune(result.hits, full_sequence=True)
            chains = [chain for chain in chains if chain]
            pbar.update(total - len(chains))
            total = len(chains)

    def to_dict(self) -> Dict[str, Any]:
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
    def fields(self) -> List[str]:
        """List keys of ``to_dict()``"""
        return list(self.to_dict().keys())


class CoverageHitChain(HitChain):
    """HitChain also calculating coverage(s)"""
    def __init__(self,
                 hits: Optional[List[BlastHit]] = None,
                 chain_penalty: int = 20) -> None:
        super().__init__(hits, chain_penalty)
        self._coverages: Mapping[str, Mapping[str, Mapping[str, str]]] = {}

    def set_coverage(self, coverages: Mapping[str, Mapping[str, Mapping[str, str]]]) -> None:
        """Set coverage data"""
        self._coverages = coverages

    def __copy__(self) -> "CoverageHitChain":
        cpy = self.__class__(self.hits, self.chain_penalty)
        cpy.set_coverage(self._coverages)
        return cpy

    @property
    def numreads(self) -> int:
        """Number of reads"""
        if self._coverages is None:
            return -1
        return sum(
            int(cov[qacc]['numreads'])
            for qacc in self.qaccs
            for cov in self._coverages.values()
        )

    def to_dict(self) -> Dict[str, Any]:
        res = super().to_dict()
        res.update({
            'numreads': self.numreads
        })
        return res
