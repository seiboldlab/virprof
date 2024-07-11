"""Calculates contig intrinsic QC values"""

import click
import csv
from statistics import mean, quantiles


from ..entropy import sequence_entropy, homopolymer_ratio
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
    klens = [1,2,3,4,5,6,8,10,12,15]
    field_names = ["acc", "frac_hp"] + [f"entropy{i}" for i in klens]
    out_writer = csv.DictWriter(out_csv, field_names)
    out_writer.writeheader()
    contigs = FastaFile(in_fasta)
    stats = {}

    for acc in contigs:
        sequence = contigs.get(acc)
        row = { "acc": acc }
        for i in klens:
            entropy = round(sequence_entropy(sequence, klen = i),2)
            row[f"entropy{i}"] = entropy
            stats.setdefault(f"entropy{i}", []).append(entropy)
        frac_hp = homopolymer_ratio(sequence, kmin = 5)
        row["frac_hp"] = frac_hp
        stats.setdefault("frac_hp", []).append(frac_hp)
        out_writer.writerow(row)

    for stat in stats:
        click.echo(f"Stat '{stat}':")
        click.echo(f" mean = {mean(stats[stat])}")
        quant = quantiles(stats[stat], n=100, method = "inclusive")
        quants = [round(quant[i],2) for i in (0,4,24,49,74,94,98)]
        click.echo(f" 1/5/25/50/75/95/99%ile = {quants}")
        click.echo(f" min/max = {min(quants)} {max(quants)}")

    return True
