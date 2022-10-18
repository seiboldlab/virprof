import logging
import gzip

import click

from .. import blast
from ._utils import group_hits_by_qacc
from ..blastbin import HitChain, greedy_select_chains
from ..fasta import filter_fasta


LOG = logging.getLogger(__name__)


@click.command()
@click.option(
    "--in-blast7",
    "in_blast7",
    type=click.File("r"),
    required=True,
    help="input blast7 format file",
)
@click.option(
    "--in-fasta",
    "in_fasta",
    type=click.File("r"),
    required=True,
    help="input blast7 format file",
)
@click.option(
    "--out",
    "-o",
    "outfile",
    type=click.File("w"),
    required=True,
    help="output blast7 format file",
)
@click.option(
    "--min-unaligned-bp", type=int, help="minimum number of unaligned basepairs"
)
@click.option("--max-aligned-bp", type=int, help="maximum number of aligned basepairs")
def cli(
    in_blast7: click.utils.LazyFile,
    in_fasta: click.utils.LazyFile,
    outfile: click.utils.LazyFile,
    min_unaligned_bp: int,
    max_aligned_bp: int,
) -> bool:
    """Filter sequences based on blast hits

    Reads fasta formatted file ``in-fasta`` and removes all sequences for
    which not at least ``min-unaligned-bp`` basepairs remain uncovered by
    blast hits.
    """
    toremove = set()
    LOG.info("Removing sequences having BLAST hits")
    if min_unaligned_bp is not None:
        LOG.info("  covering all but %i bp", min_unaligned_bp)
    if max_aligned_bp is not None:
        LOG.info("  covering no more than %i bp", max_aligned_bp)
    LOG.info("Input FASTA: %s", in_fasta.name)
    LOG.info("Input BLAST: %s", in_blast7.name)
    LOG.info("Outut FASTA: %s", outfile.name)

    LOG.info("Loading Blast7...")
    if in_blast7.name.endswith(".gz"):
        # unzip = read_from_command(["gunzip", "-dc", in_blast7.name])
        unzip = gzip.open(in_blast7.name, "rt")
        if unzip.read(1) == "":
            reader = []
        else:
            unzip.seek(0)
            reader = blast.reader(unzip)
    else:
        reader = blast.reader(in_blast7)
    hitgroups = list(group_hits_by_qacc(reader))
    LOG.info("%i distinct query sequences had matches", len(hitgroups))

    LOG.info("Making Chains...")
    for hitgroup in hitgroups:
        all_chains = HitChain().make_chains(hitgroup)
        best_chain_sets = list(greedy_select_chains(all_chains))
        best_chain_set = best_chain_sets[0]
        best_chain = best_chain_set[0]
        unaligned = best_chain.qlen - best_chain.slen
        if min_unaligned_bp is not None and unaligned < min_unaligned_bp:
            toremove.add(hitgroup[0].qacc)
        if max_aligned_bp is not None and best_chain.slen >= max_aligned_bp:
            toremove.add(hitgroup[0].qacc)

    LOG.info(
        "%i matched sequences had more than %i unaligned bp and will be kept",
        len(hitgroups) - len(toremove),
        min_unaligned_bp,
    )
    LOG.info("%i sequences will be removed", len(toremove))

    LOG.info("Filtering FASTA...")
    nin, nout = filter_fasta(in_fasta, outfile, toremove, remove=True)
    LOG.info("Read %i sequences", nin)
    LOG.info("Wrote %i sequences (%i fewer than read)", nout, nin - nout)
    return True
