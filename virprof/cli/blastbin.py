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

from typing import Callable, Iterable, List, Tuple, Optional

import click

from ..blastbin import BlastHit, HitChain, CoverageHitChain, greedy_select_chains
from ..entrez import FeatureTables, GenomeSizes
from .. import blast  # type: ignore
from ..wordscore import WordScorer
from ..taxonomy import load_taxonomy
from ._utils import group_hits_by_qacc

LOG = logging.getLogger(__name__)


def prefilter_hits_taxonomy(
    hitgroups: Iterable[List[BlastHit]], prefilter: Callable[[int], bool]
) -> List[List[BlastHit]]:
    """Filter ``hitgroups`` using ``prefilter`` function.

    The input hitgroups are expected to each comprise HSPs from the
    same query sequence. Query sequences (i.e. the group of HSPs) are
    removed based the return value of the ``prefilter`` function. For
    each query, the taxids of the matched subject sequence are fed to
    the function. If it returns `False` for more than half of the
    matches within 90% bitscore of the best match, the query/hitgroup
    is excluded from the output.

    Args:
      hitgroups: BlastHits grouped by query accession
      prefilter: function flagging NCBI taxids to keep

    Returns:
      Filtered subset of hitgroups

    """
    n_filtered = 0
    result = []
    for hitgroup in hitgroups:
        minscore = max(hit.bitscore for hit in hitgroup) * 0.9
        top_keep = [
            all(prefilter(taxid) for taxid in hit.staxids)
            for hit in hitgroup
            if hit.bitscore > minscore
        ]
        if top_keep.count(True) >= len(top_keep) / 2:
            filtered = [
                hit
                for hit in hitgroup
                if all(prefilter(taxid) for taxid in hit.staxids)
            ]
            result.append(filtered)
        else:
            n_filtered += 1
    LOG.info(f"Removed {n_filtered} contigs matching prefilter taxonomy branches")
    return result


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
    LOG.info(f"Removed {n_filtered} hits with comparatively low score")
    return result


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
    help="Samtools coverage result file",
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
    type=click.File("w"),
    help="Output CSV file containining contig details (one row per hit)",
)
@click.option(
    "--out-features",
    "--of",
    type=click.File("w"),
    help="Output CSV file containing reference feature annotation (one row per feature)"
)
@click.option(
    "--ncbi-taxonomy",
    "-t",
    type=click.Path(),
    help="Path to NCBI taxonomy (tree or raw)",
)
@click.option(
    "--include",
    "-i",
    multiple=True,
    help="Scientific name of NCBI taxonomy node identifying"
    " subtree to include in output",
)
@click.option(
    "--exclude",
    "-e",
    multiple=True,
    help="Scientific name of NCBI taxonomy node identifying"
    " subtree to exclude from output",
)
@click.option(
    "--prefilter",
    "-e",
    multiple=True,
    help="Scientific name of NCBI taxonomy node identifying"
    " subtree to exclude from input",
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
    in_coverage: click.utils.LazyFile,
    out: click.utils.LazyFile,
    out_hits: click.utils.LazyFile,
    include: Tuple[str, ...],
    exclude: Tuple[str, ...],
    prefilter: Tuple[str, ...],
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

    if out_features is not None:
        features = FeatureTables(cache_path=cache_path, api_key=ncbi_api_key)
    else:
        features = None

    if genome_size:
        genome_sizes = GenomeSizes(cache_path=cache_path, api_key=ncbi_api_key)
    else:
        genome_sizes = None

    if not in_coverage:
        chain_tpl = HitChain(chain_penalty=chain_penalty)
    else:
        chain_tpl = CoverageHitChain(chain_penalty=chain_penalty)
        cov = {}
        for cov_fd in in_coverage:
            fname = os.path.basename(cov_fd.name)
            cov_sample, _, ext = fname.rpartition(".")
            if ext != "coverage":
                LOG.warning(
                    "Parsing of filename '%s' failed. Should end in .coverage", fname
                )
            LOG.info("Loading coverage file %s", cov_sample)
            cov_reader = csv.DictReader(cov_fd, delimiter="\t")
            cov_data = {row["#rname"]: row for row in cov_reader}
            cov[cov_sample] = cov_data
            cov_fd.close()
        chain_tpl.set_coverage(cov)

    fname = os.path.basename(in_blast7.name)
    if fname.endswith(".gz"):
        sample, _, ext = fname[:-3].rpartition(".")
        if ext != "blast7":
            LOG.warning(
                "Parsing of filename '%s' failed. Should end in .blast7.gz", fname
            )
        # unzip = read_from_command(["gunzip", "-dc", in_blast7.name])
        unzip = gzip.open(in_blast7.name, "rt")
        if unzip.read(1) == "":
            reader = []
        else:
            unzip.seek(0)
            reader = blast.reader(unzip)
    else:
        sample, _, ext = fname.rpartition(".")
        if ext != "blast7":
            LOG.warning("Parsing of filename '%s' failed. Should end in .blast7")
        reader = blast.reader(in_blast7)
    LOG.info("Using sample=%s", sample)

    wordscorer = WordScorer(keepwords=num_words)
    taxonomy = load_taxonomy(ncbi_taxonomy)

    if not no_standard_excludes:
        prefilter += (
            "artificial sequences",
            "unclassified sequences",
            "uncultured bacterium",
            "uncultured eukaryote",
            "uncultured fungus",
            "uncultured phage",
            "uncultured virus",
        )
        exclude += ("Hominidae",)

    LOG.info(f"Excluding: {exclude}")

    taxfilter_pre = taxonomy.make_filter(exclude=prefilter)
    taxfilter = taxonomy.make_filter(include, exclude)

    LOG.info("Writing bins to: %s", out.name)
    fields = ["sample", "words"] + chain_tpl.fields + ["taxid"]
    if genome_sizes is not None:
        fields += ["genome_size_source", "genome_size"]
    if not taxonomy.is_null():
        fields += ["taxname", "species", "lineage", "lineage_ranks"]
    LOG.info("  output fields: %s", " ".join(fields))
    bin_writer = csv.DictWriter(out, fields)
    bin_writer.writeheader()

    LOG.info("Writing contig hit details to: %s", out_hits.name)
    fields = ["sample", "sacc"] + chain_tpl.hit_fields
    LOG.info("  output fields: %s", " ".join(fields))
    hits_writer = csv.DictWriter(out_hits, fields)
    hits_writer.writeheader()

    if features:
        LOG.info("Writing feature details to: %s", out_features.name)
        fields = features.fields
        LOG.info("  output fields: %s", " ".join(fields))
        feature_writer = csv.DictWriter(out_features, fields)
        feature_writer.writeheader()

    hitgroups = list(group_hits_by_qacc(reader))
    LOG.info(
        "Found %i contigs with %i hits",
        len(hitgroups),
        sum(len(hitgroup) for hitgroup in hitgroups),
    )

    # Filter
    if in_coverage:
        assert isinstance(chain_tpl, CoverageHitChain)
        hitgroups = list(chain_tpl.filter_hitgroups(hitgroups, min_read_count))
    filtered_hits = prefilter_hits_taxonomy(hitgroups, taxfilter_pre)
    filtered_hits = prefilter_hits_score(filtered_hits)
    LOG.info(
        "After prefiltering, %i contigs with %i hits remain",
        len(filtered_hits),
        sum(len(hitgroup) for hitgroup in filtered_hits),
    )

    all_chains = chain_tpl.make_chains(
        hit for hitgroup in filtered_hits for hit in hitgroup
    )
    best_chains = greedy_select_chains(all_chains)

    for chains in best_chains:
        taxid = wordscorer.score_taxids(chains)
        if not taxfilter(taxid):
            taxname = taxonomy.get_name(taxid)
            words = wordscorer.score(chains)
            LOG.info("Excluding %s -- %s", taxname, words)
            continue

        if features is not None:
            saccs = [chain.sacc for chain in chains]
            fdata = features.get(saccs, [])
            # FIXME: finish this - merge features from multiple references where useful

        selected_chain = chains[0]

        row = selected_chain.to_dict()
        row.update(
            {
                "sample": sample,
                "taxid": taxid,
                "words": wordscorer.score(chains),
            }
        )

        if genome_sizes is not None:
            row["genome_size_source"], row["genome_size"] = genome_sizes.get(taxid)

        if not taxonomy.is_null():
            row.update(
                {
                    "lineage": "; ".join(taxonomy.get_lineage(taxid)),
                    "lineage_ranks": "; ".join(taxonomy.get_lineage_ranks(taxid)),
                    "taxname": taxonomy.get_name(taxid),
                    "species": taxonomy.get_rank(taxid, "species"),
                }
            )
        LOG.info("Found     %s -- %s", row.get("taxname", ""), row["words"])
        bin_writer.writerow(row)

        for row in selected_chain.hits_to_dict():
            row.update(
                {
                    "sample": sample,
                    "sacc": selected_chain.sacc,
                }
            )
            hits_writer.writerow(row)

        if features is not None:
            ftable = features.get(selected_chain.sacc, ["gene/gene", "CDS/product"])
            for row in ftable:
                feature_writer.writerow(row)

    return True
