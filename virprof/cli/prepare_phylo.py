"""
Prepares sequences for phylogenetic analysis
"""

import logging
import re

import click


from ..entrez import EntrezAPI, GenomeSizes
from ..fasta import FastaFile


LOG = logging.getLogger(__name__)


@click.command()
@click.argument(
    "in-fasta",
    type=click.File("r"),
    nargs=-1,
)
@click.option(
    "--species",
    type=str,
    help="Species to fetch sequences for. May be specified multiple times. "
    "If given, species names will not be autodetected from input fasta file names",
    multiple=True,
)
@click.option(
    "--in-filter",
    type=str,
    help="Specify regex to filter FASTA input sequences by name",
)
@click.option(
    "--add-outgroup/--no-add-outgroup",
    is_flag=True,
    default=False,
    show_default=True,
    help="Not implemented",
)
@click.option(
    "--add-references/--no-add-references",
    is_flag=True,
    default=True,
    show_default=True,
    help="Add reference sequences to output",
)
@click.option(
    "--only-references",
    is_flag=True,
    default=False,
    help="Only export reference sequences",
)
@click.option(
    "--add-all-genomes/--no-add-all-genomes",
    is_flag=True,
    default=False,
    show_default=True,
    help="Add all NCBI genomes for selected species to output",
)
@click.option(
    "--include-short/--exclude-short", is_flag=True, default=False, show_default=True
)
@click.option(
    "--include-full/--exclude-full", is_flag=True, default=True, show_default=True
)
@click.option(
    "--out-fasta",
    "-o",
    type=click.File("w"),
    help="Output FASTA file for alignment+treeing",
    required=True,
)
@click.option(
    "--min-bp",
    "override_min_bp",
    type=int,
    help="Minimum number of unambious base pairs. [default: 80% detected genome size]",
)
@click.option(
    "--min-bp",
    "override_max_bp",
    type=int,
    help="Maximum number of unambious base pairs (for reference include only). [default: 120% detected genome size]",
)
@click.option("--ncbi-api-key", type=str, help="NCBI API Key")
def cli(
    in_fasta,
    in_filter,
    species,
    out_fasta,
    add_outgroup,
    add_references,
    only_references,
    add_all_genomes,
    include_short,
    include_full,
    override_min_bp,
    override_max_bp,
    ncbi_api_key,
):
    """
    Prepares sequences for phylogenetic analysis

    Output FASTA files comprise:
    - Sequences from input FASTA files matching minimum/maximum lengths
    - Sequences from Genbank for given species
    - Outgroup sequences from Genbank for given species
    """
    if add_outgroup:
        raise NotImplementedError()

    entrez = EntrezAPI(api_key=ncbi_api_key)
    genome_sizes = GenomeSizes(entrez=entrez)
    out = FastaFile(out_fasta, "w")
    references = set()
    total_count = 0

    if species and in_fasta:
        if len(species) == 1:
            species = [species[0]] * len(in_fasta)
        elif len(species) != len(in_fasta):
            raise click.UsageError(
                "Must have same number of species names and input fasta files"
            )
    elif species:
        in_fasta = [None] * len(species)
    elif in_fasta:
        species = []
        for fafd in in_fasta:
            match = re.search(r"(^|/)([^.]+)\.fasta.gz", fafd.name)
            if not match:
                raise click.UsageError(
                    "Unable to detect organism name from input fasta file. "
                    "Use --species flag to specify the species name(s) explicitly"
                )
            name = match.group(2).replace("_", " ")
            species.append(name)
    else:
        raise click.UsageError("Please specify input file(s) or species name(s)")

    for infile, name in zip(in_fasta, species):
        LOG.info("Processing file %s", infile.name)
        LOG.info("Fetching taxid for %s", name)
        taxid = genome_sizes.get_taxid(name)
        if not taxid:
            raise click.UsageError("Unknown Species")

        LOG.info("Determining genome size for %s (taxid=%s)", name, taxid)
        genome_method, genome_size = genome_sizes.get(taxid)
        if not genome_size:
            raise click.UsageError("Failed to determine genome size.")
        LOG.info("Found genome size %i using method %s", genome_size, genome_method)

        if override_min_bp is None:
            min_bp = int(genome_size * 0.8)
        else:
            min_bp = override_min_bp
        if override_max_bp is None:
            max_bp = int(genome_size * 1.2)
        else:
            max_bp = override_max_bp
        LOG.info("Full output sequences must have >= %i AGCT bases", min_bp)
        LOG.info(
            "Reference sequences must also have between %i and %i bases", min_bp, max_bp
        )

        if infile:
            genomes = FastaFile(infile)
            LOG.info("Found %i sequences in file '%s'", len(genomes), infile.name)
            count = 0
            for acc in genomes:
                sequence = genomes.get(acc)
                comment = genomes.comments.get(acc.encode("utf-8"))
                sample, _, ref_acc = acc.removesuffix("_pilon").partition(".")
                comment += f" [ref={ref_acc}] [species={name}]"
                if in_filter and not re.search(in_filter, sample):
                    LOG.info("Skipping in-filter mismatch %s - %s", sample, comment)
                    continue
                bp = sum(sequence.count(base) for base in (b"A", b"G", b"C", b"T"))
                is_short = bp < min_bp
                if is_short:
                    sample = sample + "_partial"
                if is_short and not include_short:
                    LOG.info("Skipping short %s - %s", sample, comment)
                    continue
                if not is_short and not include_full:
                    LOG.info("Skipping full %s - %s", sample, comment)
                    continue
                LOG.info("Using %s - %s", sample, comment)
                if not only_references:
                    out.put(sample, sequence, comment)
                count += 1
                if add_references:
                    references.add(ref_acc)
            LOG.info("Found %i/%i sequences with desired length", count, len(genomes))
            total_count += count

        if add_all_genomes:
            query = f"(txid{taxid}[Organism:exp]) AND ({min_bp}:{max_bp}[SLEN])"
            LOG.info("Querying Entrez Nuccore for '%s' ...", query)
            result = entrez.search("nucleotide", query)
            summaries = entrez.summary("nucleotide", result)
            LOG.info("Found %i matches", len(summaries))
            references.update(
                tag.text for tag in summaries.findall(".//*[@Name='AccessionVersion']")
            )
            LOG.info("Total accessions to fetch now %i", len(references))

    ref_no_vers = set(acc for acc in references if "." not in acc)
    references = set(acc for acc in references if "." in acc)
    LOG.info("Resolving %i accessions without version...", len(ref_no_vers))
    ref_no_vers = set(
        acc
        for acc in ref_no_vers
        if not any(ref.startswith(acc + ".") for ref in references)
    )
    LOG.info("  .. %i after removing duplicates with versioned accs", len(ref_no_vers))
    if ref_no_vers:
        query = " OR ".join(f"{acc}[ACCN]" for acc in ref_no_vers)
        result = entrez.search("nucleotide", query)
        summaries = entrez.summary("nucleotide", result)
        LOG.info(" .. found  %i matches", len(summaries))
        references.update(
            tag.text for tag in summaries.findall(".//*[@Name='AccessionVersion']")
        )

    LOG.info("Fetching %i remote sequences...", len(references))
    sequences = entrez.fetch("nucleotide", ids=list(references), rettype="fasta")
    ref_fasta = FastaFile(mode="")
    ref_fasta.load_fasta(sequences.encode("ascii").splitlines())
    found = set(iter(ref_fasta))
    if len(found) != len(references):
        print(references - found)
    LOG.info("Got %i sequences from Entrez", len(ref_fasta))
    for acc in ref_fasta:
        out.put(acc, ref_fasta.get(acc), ref_fasta.comments.get(acc.encode("utf-8")))
