"""Helper functions common to cli code"""

import os
import logging

from typing import List, Iterator, Iterable, Callable

from ..blast import BlastHit


LOG = logging.getLogger(__name__)


def group_hits_by_qacc(hits: Iterable[BlastHit]) -> Iterator[List[BlastHit]]:
    """Groups input hits into lists with same query accession

    Assumes that the hits are sorted by qacc (as is the case for BLAST)
    """
    qacc = None
    group: List[BlastHit] = []
    for hit in hits:
        if hit.qacc != qacc:
            if group:
                yield group
            group = []
            qacc = hit.qacc
        group.append(hit)
    if group:
        yield group


def as_file_name(name):
    return (
        name.replace(" ", "_")
        .replace("/", "_")
        .replace(".", "_")
        .replace("__", "_")
        .strip("_")
    )


def get_fnames_from_file(fdes):
    """Fetch newline separated filenames from a file"""
    with open(fdes, "r", encoding="utf-8") as filedes:
        fnames = [line.strip() for line in filedes.readlines()]
    missing = [fname for fname in fnames if not os.path.exists(fname)]
    if missing:
        raise RuntimeError(
            f"Missing files listed in {fdes.name}:\n" f"{', '.join(missing)}\n"
        )
    return fnames


def filter_hits(
    hitgroups: Iterable[List[BlastHit]], prefilter: Callable[[int], bool]
) -> List[List[BlastHit]]:
    """Filter BLAST search results at query level

    The input hitgroups are expected to each comprise HSPs from the
    same query sequence. Query sequences (i.e. the group of HSPs) are
    removed based the return value of the ``prefilter`` function. For
    each query, the taxids of the matched subject sequence are fed to
    the function. If it returns `False` for more than half of the
    matches within 90% bitscore of the best match, the query/hitgroup
    is excluded from the output.

    Args:
      hitgroups: BlastHits grouped by query accession
      prefilter: function flagging NCBI taxids to keep

    Returns:
      Filtered subset of hitgroups

    """
    n_filtered = 0
    result = []
    for hitgroup in hitgroups:
        minscore = max(hit.bitscore for hit in hitgroup) * 0.9
        top_keep = [
            all(prefilter(taxid) for taxid in hit.staxids)
            for hit in hitgroup
            if hit.bitscore > minscore
        ]
        if top_keep.count(True) >= len(top_keep) / 2:
            filtered = [
                hit
                for hit in hitgroup
                if all(prefilter(taxid) for taxid in hit.staxids)
            ]
            result.append(filtered)
        else:
            n_filtered += 1
    LOG.info("Removed %i contigs matching prefilter taxonomy branches", n_filtered)
    return result
