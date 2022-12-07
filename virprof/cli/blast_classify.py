"""CLI for classifying individual queries based on BLAST results

"""

import logging
import gzip
import csv
import os

from typing import Union, TextIO

import click

from ..taxonomy import load_taxonomy
from .. import blast
from ._utils import group_hits_by_qacc


LOG = logging.getLogger(__name__)


@click.command()
@click.option(
    "--in-blast7",
    "-b",
    type=click.File("rt"),
    required=True,
    help="BLAST result file in format 7. May be gzipped.",
)
@click.option(
    "--out",
    "-o",
    type=click.File("w"),
    default="-",
    help="Output CSV file containing final calls (one row per bin)",
)
@click.option(
    "--ncbi-taxonomy",
    "-t",
    required=True,
    type=click.Path(),
    help="Path to NCBI taxonomy (tree or raw)",
)
@click.option(
    "--quorum",
    "-q",
    default=0.9,
    type=float,
    help="Fraction of all hits that have to be classified to the"
    "reported depth. I.e., if 0.9, 90% of all input blast hits for a"
    "query must have been classified to species level for a species"
    "level result to be reported.",
)
@click.option(
    "--majority",
    "-m",
    type=float,
    default=0.9,
    help="The majority required for the reported result. If 0.9, the"
    "deepest path at which 90% of all hits agree will be reported.",
)
def cli(
    in_blast7: click.utils.LazyFile,
    out: click.utils.LazyFile,
    ncbi_taxonomy: str,
    quorum: float,
    majority: float,
):
    """Compute LCA classification from BLAST search result"""
    if in_blast7.name.endswith(".gz"):
        fdes: Union[TextIO, click.utils.LazyFile] = gzip.open(in_blast7.name, "rt")
    else:
        fdes = in_blast7
    reader = blast.reader(fdes)
    fieldnames = ["qacc", "score", "taxname", "genus", "lineage", "ranks"]
    writer = csv.DictWriter(out, fieldnames=fieldnames)
    writer.writeheader()
    taxonomy = load_taxonomy(ncbi_taxonomy)
    blast_results = group_hits_by_qacc(reader)
    for hitgroup in blast_results:
        lca = taxonomy.get_lca((hit.staxids[0], hit.bitscore) for hit in hitgroup)
        lca_node, lca_score = lca.classify(quorum=quorum, majority=majority)
        writer.writerow(
            {
                "qacc": hitgroup[0].qacc,
                "score": round(lca_score, 3),
                "taxname": taxonomy.get_name(lca_node),
                "genus": taxonomy.get_rank(lca_node, "genus"),
                "lineage": "; ".join(taxonomy.get_lineage(lca_node)),
                "ranks": "; ".join(taxonomy.get_lineage_ranks(lca_node)),
            }
        )
    if in_blast7.name.endswith(".gz"):
        fdes.close()
