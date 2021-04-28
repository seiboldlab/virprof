"""Merges and classifies BLAST HSPs

A BLAST search of contigs against a reference results in a large
number of HSPs. These match regions of the contigs ("query "sequences)
to regions of the reference ("subject") sequences. For the purpose of
further analysis, we instead desire a minimal set of subjects
explaining the set of contigs. This script attempts to find maximally
scoring chains of HSPs explaining the contigs with a minimal number of
subject sequences.
"""


import csv
import os
import sys
import logging
import gzip
import re

from typing import List, Tuple, Iterator, Iterable, Callable, Optional

import click
import tqdm  # type: ignore

from . import blast  # type: ignore
from .blastbin import BlastHit, HitChain, CoverageHitChain, greedy_select_chains
from .wordscore import WordScorer
from .taxonomy import load_taxonomy
from .fasta import filter_fasta, FastaFile
from .regionlist import RegionList
from .fasta import scaffold_contigs, Btop
from .entrez import FeatureTables, GenomeSizes, EntrezAPI

LOG = logging.getLogger(__name__)


# Increase CSV field size limit to 2GB
csv.field_size_limit(2 ** 31)


class TqdmHandler(logging.Handler):
    """Logging handler passing writes through tqdm

    This is used so progress bar and log messages co-habitate stdout without
    clobbering lines.

    """

    def __init__(self) -> None:
        logging.Handler.__init__(self)

    def emit(self, record: logging.LogRecord) -> None:
        try:
            msg = self.format(record)
            tqdm.tqdm.write(msg)
        except (KeyboardInterrupt, SystemExit, RecursionError):
            raise
        except Exception:  # pylint: disable=broad-except
            self.handleError(record)


def setup_logging() -> None:
    """Sets up python logging facility"""
    if getattr(setup_logging, "done", False):
        return
    setup_logging.done = True
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S",
        handlers=[TqdmHandler()],
    )


def setup_profiling() -> None:
    """Start yappi profiler

    Profiling ends and results are printed when the program ends.

    Requires "yappi" to have been installed.
    """
    # pylint: disable=import-outside-toplevel
    import yappi  # type: ignore
    import atexit

    def dump_profile() -> None:
        """Print the profile at exit"""
        stats = yappi.get_func_stats()
        stats.sort("ttot")
        stats.print_all(
            columns={
                0: ("name", 120),
                1: ("ncall", 10),
                2: ("tsub", 8),
                3: ("ttot", 8),
                4: ("tavg", 8),
            }
        )
        yappi.stop()

    atexit.register(dump_profile)
    yappi.start()


def setup_debug() -> None:
    import signal
    import traceback

    def handlesigusr1(sig, frame):
        print("------------------")
        print("Dumping stack (caught signal %i)" % sig)
        traceback.print_stack(f=frame, limit=5)
        print("------------------")

    signal.signal(signal.SIGUSR1, handlesigusr1)


def group_hits_by_qacc(hits: List[BlastHit]) -> Iterator[List[BlastHit]]:
    """Groups input hits into lists with same query accession

    Assumes that the hits are sorted by qacc (as is the case for BLAST)
    """
    qacc = None
    group: List[BlastHit] = []
    for hit in hits:
        if hit.qacc != qacc:
            if group:
                yield group
            group = []
            qacc = hit.qacc
        group.append(hit)
    if group:
        yield group


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


@click.group()
def cli() -> None:
    """Use any of the subcommands"""


@cli.command()
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
    "--in-fasta",
    "-f",
    type=click.File("r"),
    help="FASTA file with sequence data for contigs to be binned",
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
    help="Output CSV file containing reference feature annotation (one row per feature)",
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
    "--profile",
    is_flag=True,
    help="Enable performance profiling. Requires YAPPI to be installed.",
)
@click.option("--debug", is_flag=True, help="Dump stack on receipt of SIGUSR1")
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
def blastbin(
    in_blast7: click.utils.LazyFile,
    in_coverage: click.utils.LazyFile,
    in_fasta: click.utils.LazyFile,
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
    profile: bool = False,
    debug: bool = False,
    cache_path: str = "entrez_cache",
    ncbi_api_key: Optional[str] = None,
    genome_size: bool = True,
) -> bool:
    # pylint: disable=too-many-arguments
    """Merge and classify contigs based on BLAST search results"""
    LOG.info("Executing blastbin")
    if profile:
        setup_profiling()
    if debug:
        setup_debug()

    if out_features is not None:
        if not ncbi_api_key:
            # Make it so empty command line parameter equals no api key set
            ncbi_api_key = None
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
    LOG.info("Loading taxonomy from {}".format(ncbi_taxonomy))
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

    LOG.info("Excluding: {}".format(exclude))

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
        saccs = [chain.sacc for chain in chains]

        taxid = wordscorer.score_taxids(chains)
        if not taxfilter(taxid):
            taxname = taxonomy.get_name(taxid)
            words = wordscorer.score(chains)
            LOG.info("Excluding %s -- %s", taxname, words)
            continue

        if features is not None:
            fdata = features.get(saccs, [])

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


@cli.command()
@click.option("--out", "-o", type=click.File("w"), required=True, help="Output binary")
@click.option(
    "--library",
    "-b",
    type=click.Choice(["graph_tool", "networkx"]),
    default="graph_tool",
    help="Tree library to use",
)
@click.option(
    "--ncbi-taxonomy",
    type=click.Path(),
    required=True,
    help="Path to the NCBI taxonomy dump directory",
)
def index_tree(out: click.utils.LazyFile, library: str, ncbi_taxonomy: str) -> bool:
    """Parse NCBI taxonomy from dump files and write tree to binary"""
    taxonomy = load_taxonomy(ncbi_taxonomy, library=library)
    taxonomy.save_tree_binary(out.name)
    return True


@cli.command()
@click.option(
    "--in-blast7",
    "in_blast7",
    type=click.File("r"),
    required=True,
    help="input blast7 format file",
)
@click.option(
    "--in-fasta",
    "in_fasta",
    type=click.File("r"),
    required=True,
    help="input blast7 format file",
)
@click.option(
    "--out",
    "-o",
    "outfile",
    type=click.File("w"),
    required=True,
    help="output blast7 format file",
)
@click.option(
    "--min-unaligned-bp", type=int, help="minimum number of unaligned basepairs"
)
@click.option("--max-aligned-bp", type=int, help="maximum number of aligned basepairs")
def filter_blast(
    in_blast7: click.utils.LazyFile,
    in_fasta: click.utils.LazyFile,
    outfile: click.utils.LazyFile,
    min_unaligned_bp: int,
    max_aligned_bp: int,
) -> bool:
    """Filter sequences based on blast hits

    Reads fasta formatted file ``in-fasta`` and removes all sequences for
    which not at least ``min-unaligned-bp`` basepairs remain uncovered by
    blast hits.
    """
    toremove = set()
    LOG.info("Removing sequences having BLAST hits")
    if min_unaligned_bp is not None:
        LOG.info("  covering all but %i bp", min_unaligned_bp)
    if max_aligned_bp is not None:
        LOG.info("  covering no more than %i bp", max_aligned_bp)
    LOG.info("Input FASTA: %s", in_fasta.name)
    LOG.info("Input BLAST: %s", in_blast7.name)
    LOG.info("Outut FASTA: %s", outfile.name)

    LOG.info("Loading Blast7...")
    if in_blast7.name.endswith(".gz"):
        # unzip = read_from_command(["gunzip", "-dc", in_blast7.name])
        unzip = gzip.open(in_blast7.name, "rt")
        if unzip.read(1) == "":
            reader = []
        else:
            unzip.seek(0)
            reader = blast.reader(unzip)
    else:
        reader = blast.reader(in_blast7)
    hitgroups = list(group_hits_by_qacc(reader))
    LOG.info("%i distinct query sequences had matches", len(hitgroups))

    LOG.info("Making Chains...")
    for hitgroup in hitgroups:
        all_chains = HitChain().make_chains(hitgroup)
        best_chain_sets = list(greedy_select_chains(all_chains))
        best_chain_set = best_chain_sets[0]
        best_chain = best_chain_set[0]
        unaligned = best_chain.qlen - best_chain.slen
        if min_unaligned_bp is not None and unaligned < min_unaligned_bp:
            toremove.add(hitgroup[0].qacc)
        if max_aligned_bp is not None and best_chain.slen >= max_aligned_bp:
            toremove.add(hitgroup[0].qacc)

    LOG.info(
        "%i matched sequences had more than %i unaligned bp and will be kept",
        len(hitgroups) - len(toremove),
        min_unaligned_bp,
    )
    LOG.info("%i sequences will be removed", len(toremove))

    LOG.info("Filtering FASTA...")
    nin, nout = filter_fasta(in_fasta, outfile, toremove, remove=True)
    LOG.info("Read %i sequences", nin)
    LOG.info("Wrote %i sequences (%i fewer than read)", nout, nin - nout)
    return True


def as_file_name(name):
    return (
        name.replace(" ", "_")
        .replace("/", "_")
        .replace(".", "_")
        .replace("__", "_")
        .strip("_")
    )


@cli.command()
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
    help="Break scaffolds if connecting contigs would require inserting more basepairs than this number.",
)
def export_fasta(
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

    ## Handle filtering options
    if filter_lineage is not None:
        regex = re.compile(filter_lineage)
        LOG.info("Filtering bins")
        bins = {
            bin: rows for bin, rows in bins.items() if regex.match(rows[0]["lineage"])
        }
        LOG.info("  %i bins matched lineage", len(bins))

    ## Load FASTA
    LOG.info("Loading FASTA from '%s'...", in_fasta.name)
    contigs = FastaFile(in_fasta)
    LOG.info("  found %i sequences", len(contigs))

    ## Write FASTA
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


@cli.command()
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
def find_bins(in_call_files, in_fasta_files, filter_lineage, bin_by, out, out_bins):
    LOG.info("Input call files:   %s", in_call_files)
    LOG.info("Input fasta files:  %s", in_fasta_files)
    LOG.info("Filter lineage:     %s", filter_lineage)
    LOG.info("Out pattern:        %s", out)
    LOG.info("Binning by:         % s", bin_by)
    LOG.info("Bin list file:      %s", out_bins.name)

    with open(in_call_files) as filedes:
        call_fnames = [line.strip() for line in filedes.readlines()]
    LOG.info("Loading %i call files...", len(call_fnames))

    acc_to_bin = {}
    bin_to_file = {}
    for fname in call_fnames:
        with open(fname) as filedes:
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

    with open(in_fasta_files) as filedes:
        fasta_fnames = [line.strip() for line in filedes.readlines()]
    LOG.info("Processing %i fasta files...", len(fasta_fnames))

    for idx, fasta_fname in enumerate(fasta_fnames):
        LOG.info("  (%i/%i) %s", idx + 1, len(fasta_fnames), fasta_fname)
        with open(fasta_fname) as filedes:
            fasta = FastaFile(filedes, "r")
            for acc in fasta:
                _sample, _, sacc = acc.partition(".")
                sacc, _, _ = sacc.partition("_")
                if not sacc in acc_to_bin:
                    continue
                out = bin_to_file[acc_to_bin[sacc]]
                sequence = fasta.get(acc)
                bp = sum(sequence.count(base) for base in (b"A", b"G", b"C", b"T"))
                out.put(acc, fasta.get(acc), comment=f"[bp={bp}]")

    LOG.info("Closing output files...")
    for out in bin_to_file.values():
        out.close()
    LOG.info("FINISHED")

@cli.command()
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
)
@click.option("--ncbi-api-key", type=str, help="NCBI API Key")
@click.option("--out-accs", type=click.File("w"))
@click.option("--out-fasta", type=click.File("w"))
@click.option("--out-gb", type=click.File("w"))
def download_genomes(
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
    setup_logging()

    LOG.info("Loading taxonomy from %s", ncbi_taxonomy)
    taxonomy = load_taxonomy(ncbi_taxonomy)
    LOG.info("Searching for species '%s'...", species)
    taxid = taxonomy.get_taxid(species)
    if taxid is None:
        LOG.error("Species '%s' not found in taxonomy", species)
        return 1
    LOG.info("Found taxid %i", taxid)

    if ncbi_api_key is not None:
        defaults = {"api_key": ncbi_api_key}
    else:
        defaults = {}
    entrez = EntrezAPI(defaults=defaults)
    if auto_len:
        LOG.info("Determining genome size...")
        genome_sizes = GenomeSizes(entrez=entrez)
        genome_method, genome_size = genome_sizes.get(taxid)
        if not genome_size:
            LOG.info("Failed to determine genome size. Aborting.")
            return 1
        LOG.info("Found genome size %i using method %s", genome_size, genome_method)
        min_len = genome_size * 0.8
        max_len = genome_size * 1.2

    ids = []
    accs = []
    if outgroup in ("yes", "no"):
        LOG.info("Querying Entrez Nuccore for matching entries with len between %i and %i...", min_len, max_len)
        query = entrez.search(
            "nucleotide",
            f"(txid{taxid}[Organism:exp])"
            f"AND ({min_len}:{max_len}[SLEN])"
        )
        summaries = entrez.summary("nucleotide", query)
        LOG.info("Found %i matches", len(summaries))
        ids += [tag.text for tag in summaries.findall(".//Id")]
        accs += [tag.text for tag in summaries.findall(".//*[@Name='AccessionVersion']")]

    if outgroup in ("yes", "only"):
        taxids = taxonomy.get_siblings(taxid)
        LOG.info("Found %i sibling species", len(taxids))
        for taxid in taxids:
            LOG.info("  searching for reference for %i (%s)", taxid, taxonomy.get_name(taxid))
            query = entrez.search(
                "nucleotide",
                f"(txid{taxid}[Organism:exp])"
                f"AND (\"complete genome\") AND (refseq[filter])"
            )
            summaries = entrez.summary("nucleotide", query)
            if summaries.findall("ERROR"):
                LOG.info("    no results found")
                continue
            LOG.info("    found %i matches", len(summaries))
            ids += [next(iter(summaries.findall(".//Id"))).text]
            accs += [next(iter(summaries.findall(".//*[@Name='AccessionVersion']"))).text]
            LOG.info("    added first match (%s)", accs[-1])

    if out_accs:
        LOG.info("Writing accession.version to %s", out_accs.name)
        for acc in accs:
            print(acc, file=out_accs)

    if out_fasta:
        LOG.info("Writing FASTA sequences to %s", out_fasta.name)
        sequences = entrez.fetch(
            "nucleotide",
            ids=ids,
            rettype="fasta"
        )
        out_fasta.write(sequences)

    if out_gb:
        LOG.info("Writing genbank to %s", out_fasta.name)
        data = entrez.fetch(
            "nucleotide",
            ids=ids,
            rettype="gbwithparts"
        )
        out_gb.write(data)

if __name__ == "__main__":
    setup_logging()
    # pylint: disable=unexpected-keyword-arg, no-value-for-parameter
    cli(prog_name="python -m virprof")
