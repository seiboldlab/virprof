"""
CLI for scaffolding (export fasta)
"""

import csv
import logging
import re
import sys

import click

from ._utils import as_file_name
from ..fasta import FastaFile, Btop, scaffold_contigs
from ..regionlist import RegionList


LOG = logging.getLogger(__name__)


@click.command()
@click.option(
    "--in-bins",
    type=click.File("r"),
    required=True,
    help="Bins CSV from blastbin command",
)
@click.option(
    "--in-hits",
    type=click.File("r"),
    required=True,
    help="Hits CSV from blastbin command",
)
@click.option(
    "--in-fasta",
    type=click.File("r"),
    required=True,
    help="FASTA file containing contigs",
)
@click.option(
    "--out", type=str, default="out%s.fasta.gz", help="FASTA output containing bins"
)
@click.option(
    "--out-scaffolds", type=click.File("w"), help="Metadata for scaffolds (CSV)"
)
@click.option("--out-bins", type=click.File("w"), help="Output list of created bins")
@click.option(
    "--bin-by",
    type=click.Choice(["sacc", "species", "taxname"]),
    default="sacc",
    help="Field to use for binning",
)
@click.option(
    "--fasta-id-format",
    type=str,
    default="{bin_name}_{qacc} {bp} bp",
    help="Format for output FASTA header",
)
@click.option("--file-per-bin", is_flag=True, help="Create separate file for each bin")
@click.option("--filter-lineage", type=str, help="Filter by lineage prefix")
@click.option(
    "--scaffold/--no-scaffold",
    is_flag=True,
    default=True,
    help="Do not merge overlapping regions",
)
@click.option(
    "--max-fill-length",
    default=50000,
    type=int,
    help="Break scaffolds if connecting contigs would require inserting more basepairs"
    "than this number.",
)
def cli(
    in_bins,
    in_hits,
    in_fasta,
    out,
    out_bins,
    out_scaffolds,
    bin_by,
    fasta_id_format,
    file_per_bin,
    filter_lineage,
    scaffold,
    max_fill_length,
):
    """Exports blastbin hits in FASTA format"""
    # Check arguments
    if out.count("%s") > 1:
        LOG.error("No more than 1 '%s' allowed in --out")
        sys.exit(1)

    LOG.info("Writing scaffold info to '%s'", out_scaffolds.name)
    meta_writer = csv.DictWriter(
        out_scaffolds,
        [
            "acc",
            "bin",
            "sstart",
            "send",
            "qacc",
            "qstart",
            "qend",
            "qlen",
            "reversed",
            "bp",
            "scaffold",
        ],
    )
    meta_writer.writeheader()

    # Handle multi-file output options
    if file_per_bin:
        if out.count("%s") != 1:
            LOG.error("Must have exactly one '%%s' in --out when using --file-per-bin")
            sys.exit(1)

        if out_bins:
            # Make sure the file is written, even if empty.
            out_bins.write("")

        def update_outfile(name=None):
            outfile = update_outfile.outfile
            if outfile is not None:
                outfile.close()
                outfile.iofile.close()
            if name is None:
                outfile = None
            else:
                outname = as_file_name(name)
                outfile = FastaFile(open(out % outname, "w"), "w")
                if out_bins:
                    out_bins.write(outname + "\n")
                LOG.info("Writing to '%s'", outfile.name)

            update_outfile.outfile = outfile
            return outfile

        update_outfile.outfile = None
    else:
        if out.count("%s"):
            out = out.replace("%s", "")
        outfile = FastaFile(open(out, "w"), "w")
        if out_bins:
            out_bins.write(outfile.name + "\n")
        LOG.info("Writing to '%s'", outfile.name)

        def update_outfile(name=None):
            if name is None:
                outfile.close()
                outfile.iofile.close()
            return outfile

    # Load hits
    LOG.info("Loading hits from '%s'", in_hits.name)
    hits = {}
    for row in csv.DictReader(in_hits):
        regl = hits.setdefault(row["sacc"], RegionList())
        regl.add(
            int(row["sstart"]),
            int(row["send"]),
            (
                row["qacc"],
                Btop(
                    int(row["sstart"]),
                    int(row["send"]),
                    int(row["qstart"]),
                    row["btop"],
                ),
            ),
        )
    LOG.info("Found %i hits", sum(len(hit) for hit in hits.values()))

    # Load and merge bins
    LOG.info("Binning by '%s'", bin_by)
    LOG.info("Loading calls from '%s'", in_bins.name)
    bins = {}
    for row in csv.DictReader(in_bins):
        sacc = row["sacc"]
        bin_name = as_file_name(row[bin_by])
        row["regionlist"] = hits[sacc]
        bins.setdefault(bin_name, []).append(row)
    LOG.info(
        "  found %i bins in %i calls", len(bins), sum(len(bin) for bin in bins.items())
    )

    # Handle filtering options
    if filter_lineage is not None:
        regex = re.compile(filter_lineage)
        LOG.info("Filtering bins")
        bins = {
            bin: rows for bin, rows in bins.items() if regex.match(rows[0]["lineage"])
        }
        LOG.info("  %i bins matched lineage", len(bins))

    # Load FASTA
    LOG.info("Loading FASTA from '%s'...", in_fasta.name)
    contigs = FastaFile(in_fasta)
    LOG.info("  found %i sequences", len(contigs))

    # Write FASTA
    for bin_name, bin_data in bins.items():
        LOG.info("writing bin %s", bin_name)
        outfile = update_outfile(bin_name)
        for call in bin_data:
            reglist = call["regionlist"]

            # Convert call to region list
            if scaffold:
                map_dicts, sequences = scaffold_contigs(
                    reglist, contigs, max_fill_length
                )
            else:
                accs = set(data[0] for _, _, datas in reglist for data in datas)
                sequences = {acc: contigs.get(acc) for acc in accs}
                map_dicts = {
                    acc: [
                        {
                            "sstart": 1,
                            "send": len(contigs.get(acc)),
                            "qacc": 1,
                            "qend": len(contigs.get(acc)),
                            "reversed": False,
                        }
                    ]
                    for acc in accs
                }

            for acc, sequence in sequences.items():
                bp = sum(sequence.count(base) for base in (b"A", b"G", b"C", b"T"))
                outacc, _, comment = fasta_id_format.format(
                    bin_name=bin_name, acc=acc, bp=bp, **call
                ).partition(" ")
                LOG.info("writing sequence %s %s", acc, comment)
                outfile.put(outacc, sequence, comment)
                for row in map_dicts[acc]:
                    row["acc"] = outacc
                    row["bin"] = bin_name
                    meta_writer.writerow(row)

    update_outfile()
    LOG.info("done")
    return True
