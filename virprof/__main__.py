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

from typing import List, Tuple, Dict, Type, Iterator, Iterable, Callable, Optional

import click
import tqdm  # type: ignore
import ymp.blast  # type: ignore

from .blastbin import BlastHit, HitChain, CoverageHitChain
from .wordscore import WordScorer
from .taxonomy import load_taxonomy


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


def parse_file_args(files: List[click.utils.LazyFile]
                   ) -> Tuple[Type[HitChain],
                              Dict[str, Tuple[click.utils.LazyFile,
                                              Dict[str, click.utils.LazyFile]]]]:
    """Composes sample dictionary from list of files passed on cmdline

    - If we have only files ending in ``.blast7``, each is a sample of
      its own.
    - If we have an equal number of ``.blast7`` and ``.coverage``
      files, each pair is a sample and the basenames must match.
    - If we have one file ending in ``.blast7`` and a number of
      ``.coverage`` files, we have one sample with multiple coverage
      files (e.g. from a co-assembly over technical replicates or time
      series).

    Returns:
      A dictionary with sample names (file base names) as key and
      pairs of file-descriptor and coverage-file dictionary as
      value. The coverage-file dictionary itself has the coverage file
      basename as keys and the respective file-descriptor as values. E.g.

        ```
           {'sample1': (sample_1_blast7_file, {'coverage1': coverage1_fd})}
        ```
    """
    # Sort file arguments by extension
    cov_files = {}
    blast7_files = {}
    other_files = []
    for fdes in files:
        name = os.path.basename(fdes.name)
        if name.endswith('.coverage'):
            cov_files[name[:-len('.coverage')]] = fdes
        elif name.endswith('.blast7'):
            blast7_files[name[:-len('.blast7')]] = fdes
        else:
            other_files.append(name)
    if other_files:
        raise click.BadArgumentUsage(
            f"Unknown file type in arguments: {other_files}"
        )

    # Build result dictionary
    samples = {}
    if cov_files:
        chain_class: Type[HitChain] = CoverageHitChain
        if len(blast7_files) > 1:
            basenames = set(cov_files.keys()) | set(blast7_files.keys())
            missing_files = set(
                name for name in basenames
                if name not in cov_files
                or name not in blast7_files
            )
            if missing_files:
                raise click.BadArgumentUsage(
                    f"Missing either blast7 or coverage file for:"
                    f" {missing_files}"
                )
            for name in basenames:
                samples[name] = (blast7_files[name], {name: cov_files[name]})
            LOG.info("Processing %i samples with coverage", len(samples))
        else:
            name = next(iter(blast7_files.keys()))
            samples[name] = (blast7_files[name], cov_files)
            LOG.info("Processing single sample with %i coverage files",
                     len(cov_files))
    else:
        chain_class = HitChain
        for name in blast7_files:
            samples[name] = (blast7_files[name], {})
        LOG.info("Processing %i samples without coverage", len(samples))
    return chain_class, samples


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
    """Prefilter hits"""
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
@click.argument('files', type=click.File('r', lazy=True), nargs=-1)
@click.option('--out', '-o', type=click.File('w'), default='-',
              help="Output CSV file")
@click.option('--ncbi-taxonomy', '-t', type=click.Path(),
              help="Path to an unpacked NCBI taxonomy")
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
def blastbin(files: List[click.utils.LazyFile],
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
    setup_logging()
    if profile:
        setup_profiling()

    if not files:
        LOG.info("No files to process")
        return True
    chain_class, samples = parse_file_args(files)

    chain_tpl = chain_class(chain_penalty=chain_penalty)
    wordscorer = WordScorer(keepwords=num_words)
    taxonomy = load_taxonomy(ncbi_taxonomy)

    if not no_standard_excludes:
        exclude += (
            'Homo sapiens',
            'Mus musculus',
            'artificial sequences',
            'unclassified sequences',
        )
        # cellular organisms
        # Microviridae
        # Caudovirales
    prefilter = taxonomy.make_filter(exclude=exclude)
    taxfilter = taxonomy.make_filter(include, exclude)

    field_list = ['sample']
    field_list += ['words']
    if not taxonomy.is_null():
        field_list += ['taxname']
    field_list += chain_tpl.fields
    field_list += ['taxid']
    if not taxonomy.is_null():
        field_list += ['lineage']
    LOG.info("Output fields: %s", ' '.join(field_list))

    writer = csv.DictWriter(out, field_list)
    writer.writeheader()
    LOG.info("Writing to: %s", out.name)

    if len(samples) > 1:
        sample_iter = tqdm.tqdm(samples.items())
    else:
        sample_iter = samples.items()

    for sample, (sample_fd, cov_files) in sample_iter:
        LOG.info("Processing %s", sample)
        reader = ymp.blast.reader(sample_fd)
        hitgroups = group_hits_by_qacc(reader)
        filtered_hits = prefilter_hits(hitgroups, prefilter)
        load_coverage(chain_tpl, cov_files)
        all_chains = chain_tpl.make_chains(filtered_hits)
        sample_fd.close()
        best_chains = chain_tpl.greedy_select_chains(all_chains)

        for chain, altchains in best_chains:
            allchains = [chain] + (altchains or [])
            taxid = wordscorer.score_taxids(allchains)
            if not taxfilter(taxid):
                continue

            row = chain.to_dict()
            row['sample'] = sample
            row['taxid'] = taxid
            row['words'] = wordscorer.score(allchains)
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


if __name__ == "__main__":
    # pylint: disable=unexpected-keyword-arg, no-value-for-parameter
    cli(prog_name="python -m virprof")
