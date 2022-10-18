"""Downloads genomes from Entrez for a given species"""

import logging

import click

from ..entrez import EntrezAPI, GenomeSizes
from ..taxonomy import load_taxonomy

LOG = logging.getLogger(__name__)


def get_genome_size(taxid=None, organism_name=None, entrez=None):
    genome_sizes = GenomeSizes(entrez=entrez)
    if not taxid:
        if not organism_name:
            raise RuntimeError("Need either taxid or organism name")
        LOG.info("Finding taxid for %s", organism_name)
        taxid = genome_sizes.get_taxid(organism_name)
    LOG.info("Determining genome size...")

    genome_method, genome_size = genome_sizes.get(taxid)
    if not genome_size:
        raise RuntimeError("Failed to determine genome size.")
    LOG.info("Found genome size %i using method %s", genome_size, genome_method)
    min_len = genome_size * 0.8
    max_len = genome_size * 1.2
    return min_len, max_len


@click.command()
@click.option("--species", required=True)
@click.option("--outgroup", type=click.Choice(["yes", "no", "only"]), default="no")
@click.option("--auto-len/--no-auto-len", default=True)
@click.option("--min-len", default=500)
@click.option("--max-len", default=50000)
@click.option(
    "--ncbi-taxonomy",
    "-t",
    type=click.Path(),
    help="Path to NCBI taxonomy (tree or raw)",
    required=True,
)
@click.option("--ncbi-api-key", type=str, help="NCBI API Key")
@click.option("--out-accs", type=click.File("w"))
@click.option("--out-fasta", type=click.File("w"))
@click.option("--out-gb", type=click.File("w"))
def cli(
    species,
    outgroup,
    auto_len,
    min_len,
    max_len,
    ncbi_api_key,
    out_accs,
    out_fasta,
    out_gb,
    ncbi_taxonomy,
):
    entrez = EntrezAPI(api_key=ncbi_api_key)

    taxonomy = load_taxonomy(ncbi_taxonomy)
    taxid = taxonomy.get_taxid(species)
    if taxid is None:
        return 1

    if auto_len:
        min_len, max_len = get_genome_size(taxid=taxid, entrez=entrez)

    ids = []
    accs = []
    if outgroup in ("yes", "no"):
        LOG.info(
            "Querying Entrez Nuccore for matching entries with len between %i and %i.",
            min_len,
            max_len,
        )
        query = entrez.search(
            "nucleotide",
            f"(txid{taxid}[Organism:exp])" f"AND ({min_len}:{max_len}[SLEN])",
        )
        summaries = entrez.summary("nucleotide", query)
        LOG.info("Found %i matches", len(summaries))
        ids += [tag.text for tag in summaries.findall(".//Id")]
        accs += [
            tag.text for tag in summaries.findall(".//*[@Name='AccessionVersion']")
        ]

    if outgroup in ("yes", "only"):
        taxids = taxonomy.get_siblings(taxid)
        LOG.info("Found %i sibling species", len(taxids))
        for taxid in taxids:
            LOG.info(
                "  searching for reference for %i (%s)", taxid, taxonomy.get_name(taxid)
            )
            query = entrez.search(
                "nucleotide",
                f"(txid{taxid}[Organism:exp])"
                f'AND ("complete genome") AND (refseq[filter])',
            )
            summaries = entrez.summary("nucleotide", query)
            if summaries.findall("ERROR"):
                LOG.info("    no results found")
                continue
            LOG.info("    found %i matches", len(summaries))
            ids += [next(iter(summaries.findall(".//Id"))).text]
            accs += [
                next(iter(summaries.findall(".//*[@Name='AccessionVersion']"))).text
            ]
            LOG.info("    added first match (%s)", accs[-1])

    if out_accs:
        LOG.info("Writing accession.version to %s", out_accs.name)
        for acc in accs:
            print(acc, file=out_accs)

    if out_fasta:
        LOG.info("Writing FASTA sequences to %s", out_fasta.name)
        sequences = entrez.fetch("nucleotide", ids=ids, rettype="fasta")
        out_fasta.write(sequences)

    if out_gb:
        LOG.info("Writing genbank to %s", out_fasta.name)
        data = entrez.fetch("nucleotide", ids=ids, rettype="gbwithparts")
        out_gb.write(data)
