"""Helper functions common to cli code"""

import os

from typing import List, Iterator

from ..blastbin import BlastHit


def group_hits_by_qacc(hits: List[BlastHit]) -> Iterator[List[BlastHit]]:
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
