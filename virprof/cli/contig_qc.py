"""Calculates contig intrinsic QC values"""

import click
import csv

from ..entropy import sequence_entropy
from ..fasta import FastaFile

@click.command()
@click.option(
    "--in-fasta", "-i",
    type=click.File("r"),
    required=True,
    help="FASTA file containing contigs"
)
@click.option(
    "--out-csv", "-o",
    type=click.File("w"),
    required=True,
    help="Output CSV file for scores"
)
def cli(in_fasta,out_csv):
    """Calculates contig QC values"""
    out_writer = csv.DictWriter(out_csv, ["acc", "entropy5"])
    out_writer.writeheader()
    contigs = FastaFile(in_fasta)
    for acc in contigs:
        sequence = contigs.get(acc)
        entropy5 = sequence_entropy(sequence, klen = 5)
        out_writer.writerow({"acc": acc, "entropy5": entropy5})
    return True
