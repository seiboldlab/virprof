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
