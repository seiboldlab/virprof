"""Fetching Feature Table from Entrez"""

import logging

from typing import Optional, Set, Dict, Iterator

from requests import Session
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


LOG = logging.getLogger(__name__)


class EntrezAPI:
    """Caller for Entrez API methods

    Handles rate limiting by reacting to 429 error
    """

    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{tool}.fcgi"

    default_retry_codes = set(
        (
            500,  # Internal Server Error
            502,  # Bad Gateway Error
            504,  # Gateway Timeout
            429,  # Too Many Requests
        )
    )

    def __init__(self, session=None, defaults=None, timeout=300) -> None:
        if session is None:
            session = self.make_session()
        self._session = session

        if defaults is None:
            defaults = {}
        self._defaults = defaults

        self.timeout = timeout

    @classmethod
    def make_session(
        cls,
        max_retries: int = 10,
        max_redirects: int = 20,
        max_timeouts: int = 10,
        backoff_factor: float = 0.2,
        retry_codes: Optional[Set[int]] = None,
    ) -> Session:
        """Makes a persistent requests session

        The session created has a custom HTTPAdapter attached with a
        Retry object that will handle retries on typical server errors
        and honor the retry-after value returned if we exceed rate limit.
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

    def _get(self, tool: str, params: Dict[str, str]):
        url = self.URL.format(tool=tool)
        all_params = self._defaults.copy()
        all_params.update(params)
        return self._session.get(url, params=all_params, timeout=self.timeout)

    def fetch(self, database: str, ids, rettype: str, retmode: str = ""):
        """Calls Entrez Fetch API"""
        params = {"db": database, "id": ids, "rettype": rettype, "retmode": retmode}
        return self._get("efetch", params).text

    @staticmethod
    def enable_debug():
        """Enables debug logging from urllib3 library

        Requires logging to be setup (e.g. via logging.basicConfig())
        """
        logging.getLogger().setLevel(logging.DEBUG)
        requests_log = logging.getLogger("urllib3")
        requests_log.setLevel(logging.DEBUG)
        requests_log.propagate = True


class FeatureTableParsingError(Exception):
    "failed to parse"

    def __init__(self, msg, lineno, line):
        super().__init__(self)
        self.msg = msg
        self.lineno = lineno
        self.line = line


class FeatureTableParser:
    """Parses Entrez feature table format"""

    def __init__(self) -> None:
        #: Line buffer
        self.lines: Iterator[str] = iter([])
        #: Current line
        self.cur_line: str = ""
        #: Current line number
        self.lineno: int = 0

    def parser_error(self, msg: str) -> FeatureTableParsingError:
        """Constructs parsing exception"""
        return FeatureTableParsingError(msg, self.lineno, self.cur_line)

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
        self.next_line()
        result = {}
        while self.cur_line:
            acc = self.parse_header()
            features = self.parse_features()
            result[acc] = features
        return result

    def parse_header(self) -> str:
        """Parse header from current line"""
        fields = self.cur_line.split("|")
        if len(fields) != 3:
            raise self.parser_error(
                "Expected header in format '>Feature {src}|{acc}|{comment}'"
            ) from None
        self.next_line()
        return fields[1]

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
                start = int(fields[0])
                if fields[1][0] in "<>":
                    fields[1] = fields[1][1:]
                    end_open = True
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
