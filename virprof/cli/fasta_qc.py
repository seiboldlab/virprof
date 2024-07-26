"""Calculates contig intrinsic QC values"""

import click
import csv
import logging
from statistics import mean, quantiles, StatisticsError
import time

from ..entropy import sequence_entropy, homopolymer_ratio, normalize_sequence
from ..fasta import FastaFile

LOG = logging.getLogger(__name__)


@click.command()
@click.option(
    "--in-fasta",
    "-i",
    type=click.File("r"),
    required=True,
    help="FASTA file containing contigs",
)
@click.option(
    "--out-csv",
    "-o",
    type=click.File("w"),
    required=True,
    help="Output CSV file for scores",
)
@click.option(
    "--entropy-k-sizes",
    "--ek",
    type=str,
    help="Lengths of k-mers to use for entropy calculation. Comma separated",
    default="1,2,3,4,6,10",
)
@click.option(
    "--homopolymer-min-size",
    "--hmin",
    type=int,
    help="Minimum length for a repeat to be considered long homopolymer."
    "Set to 0 to disable.",
    default=5,
)
def cli(in_fasta, out_csv, entropy_k_sizes, homopolymer_min_size):
    """Calculates contig QC values"""
    LOG.info("Executing contig-qc")
    if entropy_k_sizes.strip() == "":
        klens = []
    else:
        klens = [int(klen) for klen in entropy_k_sizes.split(",")]
    field_names = ["acc", "len_all", "len_used"] + [f"entropy{i}" for i in klens]
    if homopolymer_min_size > 0:
        field_names.append("frac_hp")

    LOG.info("Writing results to: %s", out_csv.name)
    out_writer = csv.DictWriter(out_csv, field_names)
    out_writer.writeheader()
    stats = {}

    LOG.info("Loading FastQ '%s' ...", in_fasta.name)
    contigs = FastaFile(in_fasta)

    LOG.info("Counting things ...")
    last_print = time.time()
    for index, acc in enumerate(contigs):
        if time.time() - last_print > 3:
            LOG.info("Processing %i/%i: %s", index, len(contigs), acc)
            last_print = time.time()
        sequence = contigs.get(acc)
        orig_len = len(sequence)
        sequence = normalize_sequence(sequence)
        consider_len = len(sequence)

        row = {"acc": acc, "len_all": orig_len, "len_used": consider_len}
        for i in klens:
            entropy = round(sequence_entropy(sequence, klen=i), 2)
            row[f"entropy{i}"] = entropy
            stats.setdefault(f"entropy{i}", []).append(entropy)
        if homopolymer_min_size > 0:
            frac_hp = homopolymer_ratio(sequence, kmin=homopolymer_min_size)
            row["frac_hp"] = frac_hp
            stats.setdefault("frac_hp", []).append(round(frac_hp, 4))
        out_writer.writerow(row)

    LOG.info("Done")

    for stat in stats:
        click.echo(f"Stat '{stat}':")
        click.echo(f" mean = {mean(stats[stat])}")
        try:
            quant = quantiles(stats[stat], n=100, method="inclusive")
            quants = [round(quant[i], 2) for i in (0, 4, 24, 49, 74, 94, 98)]
            click.echo(f" 1/5/25/50/75/95/99%ile = {quants}")
            click.echo(f" min/max = {min(quants)} {max(quants)}")
        except StatisticsError:
            click.echo(f" (percentiles failed, probably not enough data)")

    return True
