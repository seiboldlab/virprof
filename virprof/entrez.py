"""Fetching Feature Table from Entrez"""

import csv
import logging
import json
import os
import shutil
import tempfile
import time
import xml.etree.ElementTree as ET

from random import randint
from typing import Optional, Set, Dict, Iterator, Union, List, Any, Tuple

import requests
from requests import Session
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from requests.exceptions import RequestException


LOG = logging.getLogger(__name__)


class BaseAPI:
    """Base class for error-handled http API access"""
    #: Default codes assumed to be temporary server issues solveable
    #: by retrying the request.
    default_retry_codes = set(
        (
            500,  # Generic Server error
            502,  # Bad Gateway Error
            504,  # Gateway Timeout
            429,  # Too Many Requests
        )
    )

    #: API endpoint, must be filled by inheriting classes
    URL = ""

    def __init__(
        self,
        session: Optional[Session] = None,
        defaults: Optional[Dict[str, str]] = None,
        timeout: int = 120
    ) -> None:
        self._session = self.make_session() if session is None else session
        self._defaults = {} if defaults is None else defaults
        self.timeout = timeout
        self._old_loglevel: Optional[int] = None

    def _get(self, url_params: Dict[str, str], params: Dict[str, str]) -> str:
        url = self.URL.format(**url_params)
        all_params = self._defaults.copy()
        all_params.update(params)
        response = self._session.get(url, params=all_params, timeout=self.timeout)
        response.raise_for_status()
        return response.text.strip()

    @classmethod
    def make_session(
        cls,
        max_retries: int = 50,
        max_redirects: int = 20,
        max_timeouts: int = 50,
        backoff_factor: float = 0.5,
        retry_codes: Optional[Set[int]] = None,
    ) -> Session:
        """Makes a persistent requests session

        The session created has a custom HTTPAdapter attached with a
        Retry object that will handle retries on typical server errors
        and honor the retry-after value returned if we exceed rate limit.

        Args:
          max_retries: Maximum number of times the request is retried after
            connection errors, read errors or "other" errors.
          max_redirects: Maximum number of HTTP redirects
          max_timeouts: Maximum number of times we will retry after a 429
            timeout error
          backoff_factor: Scaling factor for exponential backoff
          retry_codes: Codes considered a timeout error

        Returns:
          Session with urllib3 Retry object in HTTPadapter mounted for http(s)
          protocol.
        """
        if retry_codes is None:
            retry_codes = cls.default_retry_codes
        session = Session()
        retries = Retry(
            total=max_retries,
            connect=max_retries,
            read=max_retries,
            redirect=max_redirects,
            status=max_timeouts,
            status_forcelist=retry_codes,
            backoff_factor=backoff_factor,
            raise_on_redirect=True,
            raise_on_status=True,
            respect_retry_after_header=True,
        )
        adapter = HTTPAdapter(max_retries=retries)
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        return session

    def enable_debug(self):
        """Enables debug logging from urllib3 library

        Requires logging to be setup (e.g. via logging.basicConfig())
        """
        if self._old_loglevel is not None:
            return
        LOG.error("Configuring 'urllib3' logging: ON")

        logger = logging.getLogger()
        self._old_loglevel = logger.getEffectiveLevel()
        logger.setLevel(logging.DEBUG)

        requests_log = logging.getLogger("urllib3")
        requests_log.setLevel(logging.DEBUG)
        requests_log.propagate = True

    def disable_debug(self):
        """Disables debug logging from urllib3 library"""
        if self._old_loglevel is None:
            return
        LOG.error("Configuring 'urllib3' logging: OFF")

        logging.getLogger().setLevel(self._old_loglevel)
        self._old_loglevel = None

        requests_log = logging.getLogger("urllib3")
        requests_log.propagate = False


class EntrezAPI(BaseAPI):
    """Caller for Entrez API methods

    Handles rate limiting by reacting to 429 error

    Args:
      session: Pre-prepared session. Not that the session object
         needs to be set up to handle retries internally.
      defaults: Additional parameters to add to the URL on each
         API call. This is useful to e.g. set `apikey`.
      timeout: Maximum time to wait for any ony server response.
      retries: Outer loop retries for batch requests failing on 400, 429 and 500
    """

    #: Entrez API base URL (with tool name as `{tool}`)
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{tool}.fcgi"

    #: Default codes assumed to be temporary server issues solveable
    #: by retrying the request. These have error 500 removed as we need
    #: to handle it separately.
    default_retry_codes = set(
        (
            502,  # Bad Gateway Error
            504,  # Gateway Timeout
            429,  # Too Many Requests
            # 500 may be caused by single entry in batch request, so
            # simple retry does not help.
        )
    )

    def __init__(
        self,
        session: Optional[Session] = None,
        defaults: Optional[Dict[str, str]] = None,
        timeout: int = 120,
        retries: int = 3,
    ) -> None:
        super().__init__(session=session, defaults=defaults, timeout=timeout)
        self._retries = retries

    def _get(self, tool: str, params: Dict[str, str]) -> str:
        """Executes remote API call"""
        text = super()._get({'tool': tool}, params)
        if text.startswith("Error:"):
            LOG.error(
                "Entrez request failed with '%s'. Returning empty string instead.", text
            )
            text = ""
        return text

    def _batch_get(
        self,
        tool: str,
        parms: Dict[str, str],
        batch_param: str,
        batch_data: List[str],
        batch_size: int,
    ) -> str:
        result = []
        done = 0
        retry = 0
        cur_batch_size = batch_size
        while done < len(batch_data):
            to_get = batch_data[done : done + cur_batch_size]
            LOG.info(
                "Calling Entrez %s (%i of %i done, fetching %i)",
                tool,
                done,
                len(batch_data),
                len(to_get),
            )
            try:
                parms[batch_param] = ",".join(to_get)
                data = self._get(tool, parms)
            except RequestException as exc:
                if len(to_get) > 1:
                    # Try again with smaller batch size
                    cur_batch_size = int(cur_batch_size / 2)
                    continue
                if exc.response is None:
                    # No http error?
                    LOG.error("Unknown error '%s'- re-raising", exc)
                    raise
                if exc.response.status_code in (400, 429, 500):
                    # Handle 400, 429 and 500 by retrying again after 60-90 seconds
                    if retry < self._retries:
                        retry += 1
                        waitfor = randint(60, 90)
                        LOG.error(
                            "Error fetching data for '%s'. Retrying %i/%i in %i seconds",
                            to_get[0],
                            retry,
                            self._retries,
                            waitfor,
                        )
                        time.sleep(waitfor)
                    else:
                        LOG.error(
                            "Skipping sequence '%s' due to recurring errors", to_get[0]
                        )
                        done += 1
                        cur_batch_size = batch_size
                        retry = 0
                else:
                    LOG.error("Unknown error - re-raising")
                    raise
            else:
                result.append(data)
                done += len(to_get)
                cur_batch_size = min(batch_size, cur_batch_size * 2)
                retry = 0

        LOG.info("Calling Entrez efetch: DONE")
        return "\n".join(result)

    def search(self, database: str, term: str):
        """Calls Entrez Search (esearch) API"""
        result = self._get(
            "esearch",
            {
                "db": database,
                "term": term,
                "usehistory": "y",
                "retmode": "json",
            },
        )
        parsed = json.loads(result)
        return parsed["esearchresult"]

    def summary(self, database: str, query):
        """Calls Entrez Search (esearch) API"""
        data = self._get(
            "esummary",
            {
                "db": database,
                "query_key": query["querykey"],
                "WebEnv": query["webenv"],
            },
        )
        xml = ET.fromstring(data)
        return xml

    def fetch(
        self,
        database: str,
        ids: Union[List[str], str],
        rettype: str,
        retmode: str = "",
        batch_size: int = 100,
    ) -> str:
        """Calls Entrez Fetch (efetch) API

        See Entrez Docs for detauls on parameters.

        Args:
           database: The database to query. E.g. "nucleotide"
           ids: Entrez ID or accession or list thereof
           rettype: Type of data to return
           retmode: Format of data to return
           batch_size: Number of ids to query in one call
        Returns:
           The raw text returned by the Entrez API.
        """
        if isinstance(ids, str):
            ids = [ids]
        return self._batch_get(
            "efetch",
            {
                "db": database,
                "rettype": rettype,
                "retmode": retmode,
            },
            batch_param="id",
            batch_data=ids,
            batch_size=batch_size,
        )


class FeatureTableParsingError(Exception):
    "failed to parse"

    def __init__(self, msg: str, lineno: int, line: str, lastok: str):
        super().__init__()
        self.msg = msg
        self.lineno = lineno
        self.line = line
        self.lastok = lastok

    def __str__(self):
        return (
            f"\n"
            f"  Parser error in line {self.lineno}:\n"
            f"    {self.msg}\n"
            f"  Line: '{self.line}'\n"
            f"  Last OK: '{self.lastok}'\n"
        )


class FeatureTableParser:
    """Parses Entrez feature table format"""

    def __init__(self) -> None:
        #: Line buffer
        self.lines: Iterator[str] = iter([])
        #: Current line
        self.cur_line: str = ""
        #: Current line number
        self.lineno: int = 0
        #: Last OK accession
        self.lastok: str = ""

    def parser_error(self, msg: str) -> FeatureTableParsingError:
        """Constructs parsing exception"""
        return FeatureTableParsingError(msg, self.lineno, self.cur_line, self.lastok)

    def next_line(self) -> str:
        """Fetches next line from line buffer"""
        lineno = self.lineno
        try:
            cur_line = next(self.lines)
            lineno += 1
            # skip empty lines
            while not cur_line:
                cur_line = next(self.lines)
                lineno += 1
        except StopIteration:
            cur_line = ""
            lineno += 1
        self.cur_line = cur_line
        self.lineno = lineno
        return cur_line

    def parse(self, txt: str):
        """Parse a set of feature tables"""
        self.lines = iter(txt.splitlines())
        self.lineno = 0
        self.lastok = ""
        self.next_line()
        result = {}
        while self.cur_line:
            acc = self.parse_header()
            features = self.parse_features()
            result[acc] = features
            self.lastok = acc
        return result

    def parse_header(self) -> str:
        """Parse header from current line"""
        fields = self.cur_line.split("|")
        if len(fields) != 3:
            raise self.parser_error(
                "Expected header in format '>Feature {src}|{acc}|{comment}'"
            )
        acc, _, _version = fields[1].partition(".")
        self.next_line()
        return acc

    def parse_features(self):
        """Parse individual feature table"""
        features = []
        while self.cur_line:
            regions = self.parse_feature_regions()
            if not regions:
                break
            annotation = self.parse_feature_annotation()
            for start, end, typ, meta in regions:
                meta.update(annotation)
                features.append((start, end, typ, meta))
        return features

    def parse_feature_regions(self):
        """Parse region line(s) for feature table entry"""
        regions = []
        typ = None
        while self.cur_line:
            fields = self.cur_line.split("\t")
            meta = {}
            start_open = end_open = False
            try:
                if fields[0][0] in "<>":
                    fields[0] = fields[0][1:]
                    start_open = True
                # Lines like "19^\t20\tintron" appear, not clear what this means,
                # perhaps removed sequence?
                if fields[0][-1] == "^":
                    fields[0] = fields[0][:-1]
                start = int(fields[0])
                if fields[1][0] in "<>":
                    fields[1] = fields[1][1:]
                    end_open = True
                if fields[1][-1] == "^":
                    fields[1] = fields[1][:-1]
                end = int(fields[1])
                if start > end:
                    start, end = end, start
                    start_open, end_open = end_open, start_open
                    meta["_reversed"] = True
                if start_open:
                    meta["_left_open"] = True
                if end_open:
                    meta["_right_open"] = True
                if typ is None:
                    # region line must have type
                    typ = fields[2]
                else:
                    # subsequent lines must not have type
                    if len(fields) > 2 and fields[2]:
                        break
                regions.append((start, end, typ, meta))
            except (ValueError, IndexError):
                break
            self.next_line()
        return regions

    def parse_feature_annotation(self):
        """Parse key value part of feature"""
        annotation = {}
        while self.cur_line:
            if self.cur_line[:3] != "\t\t\t":
                break
            fields = self.cur_line.split("\t")
            try:
                key = fields[3]
            except IndexError:
                raise self.parser_error("Expected key in 4th column") from None
            try:
                value = fields[4]
            except IndexError:
                value = ""
            annotation[key] = value
            self.next_line()
        return annotation


class Cache:
    "Disk Cache"

    #: Default location of cache
    DEFAULT_PATH = "~/.cache/virprof/entrez"

    def __init__(
        self,
        path: str = None,
    ) -> None:
        if path is None:
            path = self.DEFAULT_PATH
        path = os.path.expanduser(path)
        os.makedirs(path, exist_ok=True)
        self._path = path
        LOG.info("Using cache: '%s'", path)

    def _make_path(self, cache: str, entry: any):
        entry = str(entry)
        if len(entry) < 6:
            entry += "x" * (6 - len(entry))
        path = os.path.join(self._path, cache, entry[0:2], entry[2:5])
        os.makedirs(path, exist_ok=True)
        return os.path.join(path, entry[5:])

    def get(self, cache: str, ids: List[str]) -> Dict[str, Any]:
        """Retrieve from cache"""
        result = {}
        for entry in ids:
            path = self._make_path(cache, entry)
            if os.path.exists(path):
                with open(path, "r") as cachefd:
                    result[entry] = json.load(cachefd)
        return result

    def put(self, cache: str, data: Dict[str, str]) -> None:
        """Write to cache"""
        LOG.info("writing %i entries to cache", len(data))
        tmpdir = tempfile.mkdtemp(dir=self._path)
        for entry in data:
            path = self._make_path(cache, entry)
            tmp = os.path.join(tmpdir, str(entry))
            with open(tmp, "w") as cachefd:
                json.dump(data[entry], cachefd)
            os.rename(tmp, path)
        shutil.rmtree(tmpdir)


class FeatureTables:
    """Access Entrez Feature Tables"""

    def __init__(
        self,
        entrez: EntrezAPI = None,
        parser: FeatureTableParser = None,
        cache_path: str = None,
        api_key=None,
    ) -> None:
        if api_key is not None:
            defaults = {"api_key": api_key}
        else:
            defaults = None
        self.entrez = entrez if entrez is not None else EntrezAPI(defaults=defaults)
        self.parser = parser if parser is not None else FeatureTableParser()
        self.cache = Cache(cache_path)

    def get(
        self,
        accessions: Union[str, List[str]],
        fields: Optional[List[str]] = None,
        nocache: bool = False,
    ):
        if fields is None:
            fields = ["*/*"]
        if isinstance(accessions, str):
            accessions = [accessions]
        LOG.info("Fetching feature tables for %i accessions", len(accessions))
        result = {}
        if not nocache:
            result.update(self.cache.get("features", accessions))
            LOG.info("Retrieved %i accessions from disk cache", len(result))
            accessions = [acc for acc in accessions if acc not in result]
        if accessions:
            text = self.entrez.fetch(
                "nucleotide", accessions, rettype="ft", retmode="text"
            )
            try:
                parsed = self.parser.parse(text)
            except FeatureTableParsingError as exc:
                exc.msg += "\n  Accessions: "
                exc.msg += " ".join(accessions)
                raise
            self.cache.put("features", parsed)
            result.update(parsed)
        return self.to_table(result, fields)

    @staticmethod
    def to_table(parsed, fields):
        select = {}
        for field in fields:
            typ, _, key = field.partition("/")
            select.setdefault(typ, set()).add(key)

        result = []
        for acc, annotation in parsed.items():
            for start, end, typ, data in annotation:
                if not (typ in select or "*" in select):
                    continue
                select2 = select.get(typ, set()) | select.get("*", set())
                for key, value in data.items():
                    if not (key in select2 or "*" in select2):
                        continue
                    result.append(
                        {
                            "acc": acc,
                            "start": start,
                            "end": end,
                            "typ": typ,
                            "key": key,
                            "value": value,
                        }
                    )
        return result

    @property
    def fields(_self):
        return ["acc", "start", "end", "typ", "key", "value"]

    @staticmethod
    def write_table(table, out):
        if not table:
            return
        writer = csv.DictWriter(out, fieldnames=table[0].keys())
        for row in table:
            writer.writerow(row)


class NcbiGenomeAPI(BaseAPI):
    #: Base URL for NCBI genome size check API
    URL = "https://api.ncbi.nlm.nih.gov/genome/v0/{func}"

    def __init__(
            self,
            session: Optional[Session] = None,
            defaults: Optional[Dict[str, str]] = None,
            timeout: int = 120,
    ) -> None:
        super().__init__(session=session, defaults=defaults, timeout=timeout)

    def call_expected_genome_size(self, taxid: int) -> Dict[str, str]:
        """Retrieves expected genome size for given NCBI taxid"""
        try:
            text = self._get({"func": "expected_genome_size"}, params={"species_taxid": taxid})
        except RequestException as exc:
            return {}
        genome_size_response = ET.fromstring(text)
        return {node.tag: node.text for node in genome_size_response}

    def get_genome_size(self, taxid: int) -> Optional[int]:
        """Determines the genome size using NCBI genome size check API

        Returns None if the API did not yield a result or if the result was malformed.

        Note: This will return the total genome size, spanning
          multiple molecules (segments, chromosomes). Where this is
          the case, no single accession will ever cover the entire
          genome.
        """
        results = self.call_expected_genome_size(taxid)
        exp_len = results.get("expected_ungapped_length")
        if not exp_len:
            # No result, usually means that there are fewer than 4 accepted reference
            # genomes for the taxonomy ID in question
            return None
        try:
            return int(exp_len)
        except ValueError:
            # Result, but not an integer. Let's log this, but then
            # continue with alternate ways of finding the genome size.
            LOG.error(
                "NCBI genome size check API returned malformed value '%s'", exp_len
            )
            return None


class GenomeSizes:
    """Determine genome sizes for given NCBI taxonomy IDs"""

    def __init__(
            self,
            entrez: EntrezAPI = None,
            genome: NcbiGenomeAPI = None,
            cache_path: str = None,
            api_key=None,
    ) -> None:
        if api_key is not None:
            defaults = {"api_key": api_key}
        else:
            defaults = None
        self.entrez = entrez if entrez is not None else EntrezAPI(defaults=defaults)
        self.genome = genome if genome is not None else NcbiGenomeAPI()
        self.cache = Cache(cache_path)

    def entrez_genome_search(
        self,
        taxid: int,
        refseq: bool = False,
        complete: bool = False,
    ) -> int:
        """Determines the genome size using Entrez APIs

        The size is found by querying entrez for all sequences within
        the taxomony node and taking the average length. The result
        can be limited to only refseq, and/or only sequences with
        "complete genome" in their title.

        Returns None if no results were found

        Params:
          taxid: The NCBI taxid for which a size should be found
          refseq: Query only refseq sequences
          complete: Query only sequences marked as "complete genome"

        Note: This will return the average molecule size for taxa
          whith genomes spanning multiple molecules and therefore
          multiple accessions. This is different from the result
          returned by `genome_size_check`.

        """
        term = f"(txid{taxid}[Organism:exp])"
        if refseq:
            term += "AND (refseq[filter])"
        if complete:
            term += 'AND ("complete genome")'
        query = self.entrez.search("nucleotide", term)
        summaries = self.entrez.summary("nucleotide", query)
        lengths = [
            int(item.text)
            for doc in summaries
            for item in doc
            if item.attrib.get("Name") == "Length"
        ]
        if not lengths:
            return None
        return int(sum(lengths) / len(lengths))

    def get_one(self, taxid: int) -> Tuple[str, int]:
        """Determines genome size for a taxonomy ID using several methods. See `get`"""
        expected_length = self.genome.get_genome_size(taxid)
        if expected_length:
            return "genome_size_check", expected_length
        expected_length = self.entrez_genome_search(taxid, refseq=True, complete=True)
        if expected_length:
            return "avg_refseq", expected_length
        expected_length = self.entrez_genome_search(taxid, refseq=True)
        if expected_length:
            return "avg_refseq", expected_length
        expected_length = self.entrez_genome_search(taxid, complete=True)
        if expected_length:
            return "avg_genome", expected_length
        expected_length = self.entrez_genome_search(taxid)
        if expected_length:
            return "avg_sequence", expected_length
        return "failed", 0

    def get_many(self, taxids: List[int], nocache: bool = False) -> Dict[int, Tuple[str, int]]:
        """Fetches the genome sizes for the NCBI taxonomy IDs passed in ``taxids``.

        For each taxonomy ID, four different approaches are tried. The
        returned dictionary has for each taxonomic ID a tuple listing
        the source for the genome size and the genome size itself. The
        four approaches are:

        `genome_size_check`: Using the NCBI "genome size check" API (returns total
           accross molecules, e.g. segments or chromosomes)

        `avg_refseq`: Using an Entrez search for the taxomic ID, considering only
           sequences deposited in RefSeq. This will return an average
           value if the taxon features multi-molecule genomes.

        `avg_genome`: Using an Entrez search considering only sequences with the
           string "complete genome" in the title.

        `avg_sequence`: Using an Entrez search without constraints.
        """
        result = {}
        if not nocache:
            result = self.cache.get("genome_sizes", taxids)
            taxids = [taxid for taxid in taxids if taxid not in result]
        newresult = {taxid: self.get_one(taxid) for taxid in taxids}
        self.cache.put("genome_sizes", newresult)
        result.update(newresult)
        return result

    def get(self, taxid: int, nocache: bool = False) -> Tuple[str, int]:
        return self.get_many([taxid], nocache=nocache).get(taxid, (None, None))


def main():
    """test main"""
    logging.basicConfig()  # level=logging.DEBUG)
    import sys

    genomes = GenomeSizes(cache_path="cache_path")
    for taxid, exp_size in {
            290028: 29926,
            29832: 18844,
    }.items():
        method, size = genomes.get(taxid, nocache=True)
        if exp_size * 0.9 < size < exp_size * 1.1:
            print(f"{taxid}: OK with method '{method}' ({size} close to expected {exp_size})")
        else:
            print(f"{taxid}: FAIL with method '{method}' ({size} not close to expected {exp_size})")

    features = FeatureTables(cache_path="cache_path")
    features.entrez.enable_debug()
    accs = [
        "NC_045512",
        "X64011",
        "4A1D_1",
        "AJ309573",
        "3J62_AA",
    ]
    table = features.get(accs, ["gene/gene", "CDS/product"])
    features.write_table(table, sys.stdout)
    accs = [
        "3J62_AA",
    ]
    table = features.get(accs, ["gene/gene", "CDS/product"])
    features.write_table(table, sys.stdout)


if __name__ == "__main__":
    main()
