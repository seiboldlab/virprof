import csv
import logging
import re

import click

from ._utils import as_file_name, get_fnames_from_file
from ..fasta import FastaFile


LOG = logging.getLogger(__name__)

def load_bins(fnames, filter_lineage, bin_by, out, out_bins):
    # Make sure this file is created even if it's empty
    out_bins.open()

    acc_to_bin = {}
    bin_to_file = {}
    for fname in fnames:
        with open(fname, "r", encoding="utf-8") as filedes:
            for row in csv.DictReader(filedes):
                if not re.match(filter_lineage, row["lineage"]):
                    continue
                bin_name = as_file_name(row[bin_by.lower()])
                fasta_fname = out % bin_name
                if bin_name not in bin_to_file:
                    bin_to_file[bin_name] = FastaFile(open(fasta_fname, "w"), "w")
                    LOG.info("Will write to %s", fasta_fname)
                    out_bins.write(bin_name + "\n")
                acc_to_bin[row["sacc"]] = bin_name
    LOG.info(
        "Found %i bins for %i unique reference accs", len(bin_to_file), len(acc_to_bin)
    )
    return acc_to_bin, bin_to_file


@click.command()
@click.option(
    "--in-call-files",
    type=str,
    required=True,
    help="Calls CSV from blastbin command",
)
@click.option(
    "--in-fasta-files",
    type=str,
    required=True,
    help="FASTA matching",
)
@click.option("--filter-lineage", type=str, help="Filter by lineage prefix")
@click.option(
    "--bin-by",
    type=click.Choice(["Sacc", "Species", "Taxname"]),
    default="sacc",
    help="Field to use for binning",
)
@click.option(
    "--out", type=str, default="out%s.fasta.gz", help="FASTA output containing bins"
)
@click.option("--out-bins", type=click.File("w"), help="Output list of created bins")
def cli(in_call_files, in_fasta_files, filter_lineage, bin_by, out, out_bins):
    """Collects recovered sequences into per-species files"""
    # pylint: disable=too-many-arguments
    LOG.info("Input call files:   %s", in_call_files)
    LOG.info("Input fasta files:  %s", in_fasta_files)
    LOG.info("Filter lineage:     %s", filter_lineage)
    LOG.info("Out pattern:        %s", out)
    LOG.info("Binning by:         % s", bin_by)
    LOG.info("Bin list file:      %s", out_bins.name)

    call_fnames = get_fnames_from_file(in_call_files)
    fasta_fnames = get_fnames_from_file(in_fasta_files)

    LOG.info("Loading %i call files...", len(call_fnames))
    acc_to_bin, bin_to_file = load_bins(
        call_fnames, filter_lineage, bin_by, out, out_bins
    )

    LOG.info("Processing %i fasta files...", len(fasta_fnames))
    for idx, fasta_fname in enumerate(fasta_fnames):
        LOG.info("  (%i/%i) %s", idx + 1, len(fasta_fnames), fasta_fname)
        with open(fasta_fname, "r", encoding="utf-8") as filedes:
            fasta = FastaFile(filedes, "r")
            for acc in fasta:
                _, _, sacc = acc.partition(".")
                sacc, _, _ = sacc.partition("_")
                if sacc not in acc_to_bin:
                    continue
                bin_to_file[acc_to_bin[sacc]].put(
                    acc, fasta.get(acc), comment=f"[bp={fasta.count_bp(acc)}]"
                )

    LOG.info("Closing output files...")
    for filedes in bin_to_file.values():
        filedes.close()
    LOG.info("FINISHED")
