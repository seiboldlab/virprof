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
import logging
import gzip

from typing import List, Tuple, Dict, Type, Iterator, Iterable, Callable, Optional

import click
import tqdm  # type: ignore
import ymp.blast  # type: ignore

from .blastbin import BlastHit, HitChain, CoverageHitChain, greedy_select_chains
from .wordscore import WordScorer
from .taxonomy import load_taxonomy
from .fasta import get_accs_from_fasta, filter_fasta, read_from_command


LOG = logging.getLogger(__name__)


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
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S',
        handlers=[TqdmHandler()]
    )


def setup_profiling() -> None:
    """Start yappi profiler

    Profiling ends and results are printed when the program ends.

    Requires "yappi" to have been installed.
    """
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
                4: ("tavg", 8)
            })
        yappi.stop()

    atexit.register(dump_profile)
    yappi.start()


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


def prefilter_hits(hitgroups: Iterable[List[BlastHit]],
                   prefilter: Callable[[int], bool]) -> List[BlastHit]:
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
    n_queries = n_filtered = 0
    hits = []
    for hitgroup in hitgroups:
        n_queries += 1
        minscore = max(hit.bitscore for hit in hitgroup) * 0.9
        top_keep = [
            all(prefilter(taxid) for taxid in hit.staxids)
            for hit in hitgroup
            if hit.bitscore > minscore
        ]
        keep_group = len(top_keep)/2 < top_keep.count(True)
        if keep_group:
            hits += hitgroup
        else:
            n_filtered += 1
    LOG.info("  %i query sequences (contigs) had hits",
             n_queries)
    LOG.info("  %i matched excludes", n_filtered)
    LOG.info("  %i HSPs to process", len(hits))
    return hits


def load_coverage(chain_tpl: HitChain, cov_files: Dict[str, click.utils.LazyFile]) -> None:
    """Load coverage into chain template"""
    if not isinstance(chain_tpl, CoverageHitChain):
        return
    cov = {}
    for cov_sample, cov_fd in cov_files.items():
        LOG.info("Loading coverage file %s", cov_sample)
        cov_reader = csv.DictReader(cov_fd, delimiter="\t")
        cov_data = {row['#rname']: row for row in cov_reader}
        cov[cov_sample] = cov_data
        cov_fd.close()
    if cov:
        chain_tpl.set_coverage(cov)

@click.group()
def cli() -> None:
    """Use any of the subcommands"""


@cli.command()
@click.option('--in-blast7', '-b', type=click.File('r'), required=True,
              help="BLAST result file in format 7")
@click.option('--in-coverage', '-c', type=click.File('r'), multiple=True,
              help="Samtools coverage result file")
@click.option('--in-fasta', '-f', type=click.File('r'),
              help="FASTA file with contigs")
@click.option('--out', '-o', type=click.File('w'), default='-',
              help="Output CSV file")
@click.option('--ncbi-taxonomy', '-t', type=click.Path(),
              help="Path to NCBI taxonomy (tree or raw)")
@click.option('--include', '-i', multiple=True,
              help="Scientific name of NCBI taxonomy node identifying"
              " subtree to include in output")
@click.option('--exclude', '-e', multiple=True,
              help="Scientific name of NCBI taxonomy node identifying"
              " subtree to exclude from output")
@click.option('--no-standard-excludes', '-E', type=bool,
              help="Do not exclude Human, artificial and unclassified"
              " sequences by default", default=False)
@click.option('--chain-penalty', type=int, default=20,
              help="Cost to BLAST score for concatenating two hits")
@click.option('--num-words', type=int, default=4,
              help="Number of words to add to 'words' field")
@click.option('--profile', is_flag=True)
def blastbin(in_blast7: click.utils.LazyFile,
             in_coverage: click.utils.LazyFile,
             in_fasta: click.utils.LazyFile,
             out: click.utils.LazyFile,
             include: Tuple[str, ...],
             exclude: Tuple[str, ...],
             ncbi_taxonomy: Optional[str] = None,
             no_standard_excludes: bool = False,
             chain_penalty: int = 20,
             num_words: int = 4,
             profile: bool = False) -> bool:
    # pylint: disable=too-many-arguments
    """Merge and classify contigs based on BLAST search results

    """
    if profile:
        setup_profiling()

    if in_coverage:
        chain_tpl = CoverageHitChain(chain_penalty=chain_penalty)
    else:
        chain_tpl = HitChain(chain_penalty=chain_penalty)

    sample = os.path.basename(in_blast7.name)
    if in_blast7.name.endswith(".gz"):
        sample = sample.rstrip(".blast7.gz")
        #unzip = read_from_command(["gunzip", "-dc", in_blast7.name])
        unzip = gzip.open(in_blast7.name)
        reader = ymp.blast.reader(unzip)
    else:
        sample = sample.rstrip(".blast7")
        reader = ymp.blast.reader(in_blast7)

    wordscorer = WordScorer(keepwords=num_words)
    LOG.info("Loading taxonomy from {}".format(ncbi_taxonomy))
    taxonomy = load_taxonomy(ncbi_taxonomy)

    if not no_standard_excludes:
        exclude += (
            'Homo sapiens',
            'Mus musculus',
            'artificial sequences',
            'unclassified sequences',
        )
    LOG.info("Excluding: {}".format(exclude))

    prefilter = taxonomy.make_filter(exclude=exclude)
    taxfilter = taxonomy.make_filter(include, exclude)

    field_list = ['sample', 'words']
    if not taxonomy.is_null():
        field_list += ['taxname']
    field_list += chain_tpl.fields
    field_list += ['taxid', 'saccs']
    if not taxonomy.is_null():
        field_list += ['lineage']
    LOG.info("Output fields: %s", ' '.join(field_list))

    writer = csv.DictWriter(out, field_list)
    writer.writeheader()
    LOG.info("Writing to: %s", out.name)

    hitgroups = group_hits_by_qacc(reader)
    filtered_hits = prefilter_hits(hitgroups, prefilter)
    load_coverage(chain_tpl, in_coverage)
    all_chains = chain_tpl.make_chains(filtered_hits)
    best_chains = greedy_select_chains(all_chains)

    for chains in best_chains:
        taxid = wordscorer.score_taxids(chains)
        if not taxfilter(taxid):
            continue
        row = chains[0].to_dict()
        row['sample'] = sample
        row['taxid'] = taxid
        row['words'] = wordscorer.score(chains)
        row['saccs'] = " ".join(chain.sacc for chain in chains)
        if not taxonomy.is_null():
            row['lineage'] = taxonomy.get_lineage(taxid)
            row['taxname'] = taxonomy.get_name(taxid)
        writer.writerow(row)
    return True


@cli.command()
@click.option('--out', '-o', type=click.File('w'), required=True,
              help="Output binary")
@click.option('--library', '-b', type=click.Choice(['graph_tool', 'networkx']),
              default='graph_tool',
              help="Tree library to use")
@click.option('--ncbi-taxonomy', type=click.Path(), required=True,
              help="Path to the NCBI taxonomy dump directory")
def index_tree(out: click.utils.LazyFile,
               library: str,
               ncbi_taxonomy: str) -> bool:
    """Parse NCBI taxonomy from dump files and write tree to binary"""
    taxonomy = load_taxonomy(ncbi_taxonomy, library=library)
    taxonomy.save_tree_binary(out.name)
    return True


@cli.command()
@click.option("--in-blast7", "in_blast7", type=click.File('r'), required=True,
              help="input blast7 format file")
@click.option("--in-fasta", "in_fasta", type=click.File('r'), required=True,
              help="input blast7 format file")
@click.option("--out", "-o", "outfile", type=click.File('w'), required=True,
              help="output blast7 format file")
@click.option("--min-unaligned-bp", type=int, default=0,
              help="minimum number of unaligned basepairs")
def filter_blast(in_blast7: click.utils.LazyFile,
                 in_fasta: click.utils.LazyFile,
                 outfile: click.utils.LazyFile,
                 min_unaligned_bp: int) -> bool:
    toremove = set()

    LOG.info("Loading Blast7")
    if in_blast7.name.endswith(".gz"):
        #unzip = read_from_command(["gunzip", "-dc", in_blast7.name])
        unzip = gzip.open(in_blast7.name, "rt")
        reader = ymp.blast.reader(unzip)
    else:
        reader = ymp.blast.reader(in_blast7)
    hitgroups = list(group_hits_by_qacc(reader))
    LOG.info("%i distinct query sequences had matches", len(hitgroups))

    LOG.info("Making Chains")
    for hitgroup in hitgroups:
        all_chains = HitChain().make_chains(hitgroup)
        best_chains = HitChain().greedy_select_chains(all_chains)
        for chain, _ in best_chains:
            remain = chain.qlen - chain.alen
            if remain < min_unaligned_bp:
                toremove.add(hitgroup[0].qacc)
    LOG.info("%i query sequences had less than %i bp unaligned and will be removed",
             len(toremove), min_unaligned_bp)
    LOG.info("(i.e. %i query sequences were matched but are kept)",
             len(hitgroups) - len(toremove))

    LOG.info("Filtering FASTA")
    filter_fasta(in_fasta, outfile, toremove, remove=True)
    LOG.info("Done")


if __name__ == "__main__":
    setup_logging()
    # pylint: disable=unexpected-keyword-arg, no-value-for-parameter
    cli(prog_name="python -m virprof")
