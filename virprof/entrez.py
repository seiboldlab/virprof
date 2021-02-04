"""Fetching Feature Table from Entrez"""

import csv
import logging
import json
import os
import shutil
import tempfile
import time

from typing import Optional, Set, Dict, Iterator, Union, List, Any

from requests import Session
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from requests.exceptions import RequestException


LOG = logging.getLogger(__name__)


class EntrezAPI:
    """Caller for Entrez API methods

    Handles rate limiting by reacting to 429 error

    Args:
      session: Pre-prepared session. Not that the session object
         needs to be set up to handle retries internally.
      defaults: Additional parameters to add to the URL on each
         API call. This is useful to e.g. set `apikey`.
      timeout: Maximum time to wait for any ony server response.
    """

    #: Entrez API base URL (with tool name as `{tool}`)
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{tool}.fcgi"

    #: Default codes assumed to be temporary server issues solveable
    #: by retrying the request.
    default_retry_codes = set(
        (
            500,  # Internal Server Error
            502,  # Bad Gateway Error
            504,  # Gateway Timeout
            429,  # Too Many Requests
        )
    )

    def __init__(
        self,
        session: Optional[Session] = None,
        defaults: Optional[Dict[str, str]] = None,
        timeout: int = 300,
        retry_unexpected = 10,
    ) -> None:
        self._session = self.make_session() if session is None else session
        self._defaults = {} if defaults is None else defaults
        self._retry_unexpected = retry_unexpected
        self.timeout = timeout
        self._old_loglevel: Optional[int] = None

    @classmethod
    def make_session(
        cls,
        max_retries: int = 100,
        max_redirects: int = 20,
        max_timeouts: int = 100,
        backoff_factor: float = 0.2,
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

    def _get(self, tool: str, params: Dict[str, str]) -> str:
        """Executes remote API call"""
        LOG.info("Calling Entrez API for %s...", tool)
        url = self.URL.format(tool=tool)
        all_params = self._defaults.copy()
        all_params.update(params)
        for n in range(self._retry_unexpected):
            try:
                request = self._session.get(url, params=all_params,
                                            timeout=self.timeout)
                text = request.text.strip()
                break
            except RequestException as exc:
                LOG.error("Unexpected exception %s from requests module: '%s'",
                          type(exc), str(exc))
                LOG.error("Retrying ... (attempt %i/%i)",
                          n+1, self._retry_unexpected)
                time.sleep(5)
                self.enable_debug()
        else:
            LOG.error("Entrez request failed after %i retries (see above).")
            raise exc
        self.disable_debug()
        if text.startswith("Error:"):
            LOG.info("Entrez request failed with '%s'. Returning empty string instead.", text)
            text = ""
        LOG.info("Calling Entrez API for %s: DONE.", tool)
        return text

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
        result = []
        for start in range(0, len(ids), batch_size):
            idparm = ids[start:start+batch_size]
            LOG.info("Calling Entrez efetch (%i..%i of %i) ...",
                     start+1, start+len(idparm), len(ids))
            result.append(self._get(
                "efetch",
                {
                    "db": database,
                    "id": ",".join(idparm),
                    "rettype": rettype,
                    "retmode": retmode,
                },
            ))
        return "\n".join(result)

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

    def _make_path(self, cache: str, entry: str):
        if len(entry) < 6:
            entry += "x" * (6 - len(entry))
        path = os.path.join(
            self._path, cache, entry[0:2], entry[2:5]
        )
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
            tmp = os.path.join(tmpdir, entry)
            with open(tmp, "w") as cachefd:
                json.dump(data[entry], cachefd)
            os.rename(tmp, path)
        shutil.rmtree(tmpdir)


class FeatureTables:
    """Access Entrez Feature Tables"""
    def __init__(
        self, entrez: EntrezAPI = None,
        parser: FeatureTableParser = None,
        cache_path: str = None,
        api_key = None
    ) -> None:
        if api_key is not None:
            defaults = {'api_key': api_key}
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
            text = self.entrez.fetch("nucleotide", accessions, rettype="ft", retmode="text")
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
            for start, stop, typ, data in annotation:
                if not (typ in select or '*' in select):
                    continue
                select2 = select.get(typ, set()) | select.get('*', set())
                for key, value in data.items():
                    if not (key in select2 or '*' in select2):
                        continue
                    result.append({
                        'acc': acc,
                        'start': start,
                        'stop': stop,
                        'typ': typ,
                        'key': key,
                        'value': value,
                    })
        return result

    @staticmethod
    def write_table(table, out):
        if not table:
            return
        writer = csv.DictWriter(out, fieldnames=table[0].keys())
        for row in table:
            writer.writerow(row)


def main():
    """test main"""
    logging.basicConfig()
    import sys
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
