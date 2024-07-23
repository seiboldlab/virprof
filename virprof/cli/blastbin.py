"""Merges and classifies BLAST HSPs

A BLAST search of contigs against a reference results in a large
number of HSPs. These match regions of the contigs ("query "sequences)
to regions of the reference ("subject") sequences. For the purpose of
further analysis, we instead desire a minimal set of subjects
explaining the set of contigs. This script attempts to find maximally
scoring chains of HSPs explaining the contigs with a minimal number of
subject sequences.
"""

import logging
import os
import csv
import gzip

from typing import Iterable, List, Tuple, Optional, Callable, Dict, Any

import click

from ..blastbin import (
    BlastHit,
    HitChain,
    CoverageHitChain,
    CoverageFastAQcHitChain,
    FastAQcHitChain,
    greedy_select_chains,
)
from ..entrez import FeatureTables, GenomeSizes
from .. import blast  # type: ignore
from ..wordscore import WordScorer
from ..taxonomy import load_taxonomy
from ._utils import group_hits_by_qacc, filter_hits

LOG = logging.getLogger(__name__)


def prefilter_blast_hits(
    hitgroups: List[List[BlastHit]],
    prefilter: Tuple[str, ...],
    prefilter_contigs: Tuple[str, ...],
    no_standard_excludes: bool,
    taxonomy,
    contig_lcas,
):
    """Applies filtering to input blast hits"""
    # LCA based filtering of entire contigs
    if not no_standard_excludes:
        prefilter_contigs += ("Euteleostomi",)
    if "None" in prefilter_contigs:
        prefilter_contigs = ()
    taxfilter_preclass = taxonomy.make_filter(exclude=prefilter_contigs)
    filtered_hits = [
        hitgroup
        for hitgroup in hitgroups
        if taxfilter_preclass(contig_lcas[hitgroup[0].qacc].classify()[0])
    ]
    LOG.info(
        "Removed %i contigs classified into exclusion taxa",
        len(hitgroups) - len(filtered_hits),
    )

    # Prefiltering of using top taxonomy hits
    if not no_standard_excludes:
        # Prune uninformative hits
        prefilter += (
            "artificial sequences",
            "unclassified sequences",
            "uncultured bacterium",
            "uncultured eukaryote",
            "uncultured fungus",
            "uncultured phage",
            "uncultured virus",
        )
    if "None" in prefilter:
        prefilter = ()
    filtered_hits = filter_hits(filtered_hits, taxonomy.make_filter(exclude=prefilter))

    # Prefiltering based on score
    filtered_hits = prefilter_hits_score(filtered_hits)
    LOG.info(
        "After prefiltering, %i contigs with %i hits remain",
        len(filtered_hits),
        sum(len(hitgroup) for hitgroup in filtered_hits),
    )

    return filtered_hits


def prefilter_hits_score(hitgroups: Iterable[List[BlastHit]]):
    """Remove hits with significantly lower score than best hit"""
    n_filtered = 0
    result = []
    for hitgroup in hitgroups:
        result_group = []
        hit_iter = iter(hitgroup)
        hitsets = [[next(hit_iter)]]
        for hit in hit_iter:
            if hitsets[-1][-1].sacc != hit.sacc:
                hitsets.append([])
            hitsets[-1].append(hit)
        best_score = sum(hit.score for hit in hitsets[0])
        best_pident = sum(hit.pident for hit in hitsets[0]) / len(hitsets[0])
        min_pident = 100 - (100 - best_pident) * 1.5 - 0.5
        min_score = best_score * 0.9
        for hitset in hitsets:
            cur_score = sum(hit.score for hit in hitset)
            cur_pident = sum(hit.pident for hit in hitset) / len(hitset)
            if cur_score > min_score and cur_pident > min_pident:
                result_group.extend(hitset)
            else:
                n_filtered += len(hitset)
        result.append(result_group)
    LOG.info("Removed %i hits with comparatively low score", n_filtered)
    return result


def classify_contigs(hitgroups: Iterable[List[BlastHit]], taxonomy):
    """Classify each input contig"""
    result = {}
    for hitgroup in hitgroups:
        taxids = set()
        for hit in hitgroup:
            if not hit.staxids:
                LOG.warning("Hit to %s on %s had no taxids", hit.sacc, hit.qacc)
                continue
            taxids.add((hit.staxids[0], hit.bitscore))
        lca = taxonomy.get_lca(taxids)
        result[hitgroup[0].qacc] = lca
    return result


def make_genome_size_fetch_function(
    cache_path: str, ncbi_api_key: Optional[str], genome_size: bool, fields: List[str]
) -> Callable[[Dict[str, Any], int], None]:
    """Makes a function to add genome sizes to result if requested"""
    if not genome_size:

        def dummy_fetcher(_row: Dict[str, Any], _taxid: int) -> None:
            return

        return dummy_fetcher

    api = GenomeSizes(cache_path=cache_path, api_key=ncbi_api_key)
    fields.extend(["genome_size_source", "genome_size"])

    def genome_size_fetcher(row: Dict[str, Any], taxid: int) -> None:
        row["genome_size_source"], row["genome_size"] = api.get(taxid)

    return genome_size_fetcher


def make_taxonomy_annotate_function(
    taxonomy, fields: List[str]
) -> Callable[[Dict[str, Any], int], None]:
    """Makes a function to add taxonomy information to result"""
    fields += ["taxname", "species", "lineage", "lineage_ranks"]

    def taxonomy_annotator(row: Dict[str, Any], taxid: int) -> None:
        row.update(
            {
                "lineage": "; ".join(taxonomy.get_lineage(taxid)),
                "lineage_ranks": "; ".join(taxonomy.get_lineage_ranks(taxid)),
                "taxname": taxonomy.get_name(taxid),
                "species": taxonomy.get_rank(taxid, "species"),
            }
        )

    return taxonomy_annotator


def make_hitchain_template(
    in_coverage: Tuple[click.utils.LazyFile, ...],
    in_fastaqc: click.utils.LazyFile,
    chain_penalty: int,
) -> HitChain:
    """Creates hitchain template object

    Return can be a CoverageHitChain if we have coverages, otherwise
    it's a regular HitChain.
    """
    if in_coverage:
        if in_fastaqc:
            klass = CoverageFastAQcHitChain
        else:
            klass = CoverageHitChain
    else:
        if in_fastaqc:
            klass = FastAQcHitChain
        else:
            klass = HitChain

    chain_tpl = klass(chain_penalty=chain_penalty)

    if isinstance(chain_tpl, CoverageHitChain):
        chain_tpl.load_coverage(in_coverage)

    if isinstance(chain_tpl, FastAQcHitChain):
        chain_tpl.load_fastaqc(in_fastaqc)

    return chain_tpl


def make_writer(
    description: str, lfile: click.utils.LazyFile, fields: List[str]
) -> csv.DictWriter:
    """Makes a dict writer and logs about it"""
    LOG.info("Writing %s to: %s", description, lfile.name)
    LOG.info("  output fields: %s", " ".join(fields))
    writer = csv.DictWriter(lfile, fields)
    writer.writeheader()
    return writer


def create_blast_reader(
    in_blast7: click.utils.LazyFile,
) -> Tuple[str, Iterable[BlastHit]]:
    """Open input Blast7 file and determine sample name"""
    fname = os.path.basename(in_blast7.name)
    if fname.endswith(".gz"):
        sample, _, ext = fname[:-3].rpartition(".")
        if ext != "blast7":
            LOG.warning(
                "Parsing of filename '%s' failed. Should end in .blast7.gz", fname
            )
        unzip = gzip.open(in_blast7.name, "rt")
        if unzip.read(1) == "":
            reader: Iterable[BlastHit] = []
        else:
            unzip.seek(0)
            reader = blast.reader(unzip)
    else:
        sample, _, ext = fname.rpartition(".")
        if ext != "blast7":
            LOG.warning("Parsing of filename '%s' failed. Should end in .blast7")
        reader = blast.reader(in_blast7)
    LOG.info("Using sample=%s", sample)
    return sample, reader


@click.command()
@click.option(
    "--in-blast7",
    "-b",
    type=click.File("r"),
    required=True,
    help="BLAST result file in format 7",
)
@click.option(
    "--in-coverage",
    "-c",
    type=click.File("r"),
    multiple=True,
    help="Samtools coverage result files",
)
@click.option(
    "--in-fastaqc",
    type=click.File("r"),
    help="VirProf FastA QC reports (entropies, homopolymers)",
)
@click.option(
    "--out",
    "-o",
    type=click.File("w"),
    default="-",
    help="Output CSV file containing final calls (one row per bin)",
)
@click.option(
    "--out-hits",
    "--oc",
    required=True,
    type=click.File("w"),
    help="Output CSV file containining contig details (one row per hit)",
)
@click.option(
    "--out-features",
    "--of",
    type=click.File("w"),
    help="Output CSV file containing reference feature annotation"
    " (one row per feature)",
)
@click.option(
    "--ncbi-taxonomy",
    "-t",
    type=click.Path(),
    required=True,
    help="Path to NCBI taxonomy (tree or raw)",
)
@click.option(
    "--exclude",
    "-e",
    multiple=True,
    help="Add NCBI taxonomy scientific names to list of taxa"
    " omitted from final results."
    " Use 'None' to disable default (Hominidae)",
)
@click.option(
    "--include",
    "-i",
    multiple=True,
    help="Add NCBI taxonomy scientific names to list of taxa"
    " included in final results.",
)
@click.option(
    "--prefilter-hits",
    multiple=True,
    help="Add NCBI taxonomy scientific names to list of taxa"
    " filtered from input BLAST search results."
    " If within 90% bitscore of the top hit more than half"
    " the hits are removed by this filter, the entire contig"
    " is excluded."
    " Use 'None' to disable default (various artificial,"
    " unclassified and uncultured taxonomy nodes).",
)
@click.option(
    "--prefilter-contigs",
    multiple=True,
    help="Add NCBI taxonomy scientific names to list of taxa"
    " filtered from input contigs as classified by relaxed LCA"
    " on input BLAST results."
    " Use 'None' to disable default (Euteleostomi).",
)
@click.option(
    "--no-standard-excludes",
    "-E",
    type=bool,
    help="Do not exclude Human, artificial and unclassified" " sequences by default",
    default=False,
)
@click.option(
    "--min-read-count",
    type=int,
    default=2,
    help="Exclude contigs with less than this number of reads."
    " Requires --in-coverage to be set to take effect.",
)
@click.option(
    "--chain-penalty",
    type=int,
    default=20,
    help="Cost to BLAST score for concatenating two hits",
)
@click.option(
    "--num-words", type=int, default=4, help="Number of words to add to 'words' field"
)
@click.option(
    "--cache-path",
    type=click.Path(),
    help="Location for caching of remote data",
    default="entrez_cache",
)
@click.option("--ncbi-api-key", type=str, help="NCBI API Key")
@click.option(
    "--genome-size/--no-genome-size",
    default=True,
    help="Obtain genome sizes for bins from NCBI",
)
def cli(
    in_blast7: click.utils.LazyFile,
    in_coverage: Tuple[click.utils.LazyFile, ...],
    in_fastaqc: click.utils.LazyFile,
    out: click.utils.LazyFile,
    out_hits: click.utils.LazyFile,
    include: Tuple[str, ...],
    exclude: Tuple[str, ...],
    prefilter_hits: Tuple[str, ...],
    prefilter_contigs: Tuple[str, ...],
    out_features: click.utils.LazyFile = None,
    ncbi_taxonomy: Optional[str] = None,
    no_standard_excludes: bool = False,
    min_read_count=2,
    chain_penalty: int = 20,
    num_words: int = 4,
    cache_path: str = "entrez_cache",
    ncbi_api_key: Optional[str] = None,
    genome_size: bool = True,
) -> bool:
    # pylint: disable=too-many-arguments
    """Merge and classify contigs based on BLAST search results"""
    LOG.info("Executing blastbin")

    sample, reader = create_blast_reader(in_blast7)
    taxonomy = load_taxonomy(ncbi_taxonomy)
    wordscorer = WordScorer(keepwords=num_words)
    chain_tpl = make_hitchain_template(in_coverage, in_fastaqc, chain_penalty)

    fields = chain_tpl.fields
    add_genome_sizes = make_genome_size_fetch_function(
        cache_path, ncbi_api_key, genome_size, fields
    )
    add_taxonomy_info = make_taxonomy_annotate_function(taxonomy, fields)

    fields = ["sample", "words"] + fields + ["taxid"]

    bin_writer = make_writer("bins", out, fields)
    hits_writer = make_writer(
        "contig hit details", out_hits, ["sample", "sacc"] + chain_tpl.hit_fields
    )

    # Set up feature output writer
    if out_features is not None:
        features = FeatureTables(cache_path=cache_path, api_key=ncbi_api_key)
        feature_writer = make_writer("feature details", out_features, features.fields)
    else:
        features = None

    # Load BLAST hits
    LOG.info("Loading blast result...")
    hitgroups = list(group_hits_by_qacc(reader))
    LOG.info(
        "Found %i contigs with %i hits",
        len(hitgroups),
        sum(len(hitgroup) for hitgroup in hitgroups),
    )

    # Filter by minimum read coverage
    if isinstance(chain_tpl, CoverageHitChain):
        hitgroups = list(chain_tpl.filter_hitgroups(hitgroups, min_read_count))

    # Classify each contig
    LOG.info("Classifying each contig...")
    contig_lcas = classify_contigs(hitgroups, taxonomy)

    # Apply hit/contig prefiltering
    filtered_hits = prefilter_blast_hits(
        hitgroups,
        prefilter_hits,
        prefilter_contigs,
        no_standard_excludes,
        taxonomy,
        contig_lcas,
    )

    # Prepare output filter
    if not no_standard_excludes:
        # Exclude results hitting humans
        exclude += ("Hominidae",)
    if "None" in exclude:
        exclude = ()
    LOG.info("Excluding from final results: %s", " ".join(exclude))
    LOG.info("Including in final results: %s", " ".join(include))
    taxfilter = taxonomy.make_filter(include, exclude)

    # Prepare generator for all chains
    all_chains = chain_tpl.make_chains(
        hit for hitgroup in filtered_hits for hit in hitgroup
    )
    # Select best chains greedily
    best_chains = list(greedy_select_chains(all_chains))

    # Run through and output result for each selected chain
    for chains in best_chains:
        lca = taxonomy.get_lca(taxid for chain in chains for taxid in chain.staxids)
        taxid, _score = lca.classify(quorum=0.9, majority=0.9)
        words = wordscorer.score(chains)

        if not taxfilter(taxid):
            LOG.info("Excluding %s -- %s", taxonomy.get_name(taxid), words)
            continue

        top_chain = chains[0]
        row = top_chain.to_dict()
        row.update({"sample": sample, "taxid": taxid, "words": words})
        add_genome_sizes(row, taxid)
        add_taxonomy_info(row, taxid)

        LOG.info("Found     %s -- %s", row.get("taxname", ""), row["words"])
        bin_writer.writerow(row)

        for row in top_chain.hits_to_dict():
            row.update({"sample": sample, "sacc": top_chain.sacc})
            hits_writer.writerow(row)

        if features is not None:
            ftable = features.get(top_chain.sacc, ["gene/gene", "CDS/product"])
            for row in ftable:
                feature_writer.writerow(row)

    return True
