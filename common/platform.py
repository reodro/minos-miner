"""Read-only client for the public Minos API. Submitting lives in main.submitter."""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from typing import Optional

from common.settings import SETTINGS


class Platform:
    def __init__(self, base_url: Optional[str] = None, timeout: float = 30.0):
        self.base = (base_url or SETTINGS.platform_url).rstrip("/")
        self.timeout = timeout

    def _get(self, path: str, **query):
        url = self.base + path + ("?" + urllib.parse.urlencode(query) if query else "")
        with urllib.request.urlopen(url, timeout=self.timeout) as r:
            return json.load(r)

    @staticmethod
    def _rid(round_id: str) -> str:
        return urllib.parse.quote(round_id, safe="")

    def parameter_ranges(self) -> dict:
        return self._get("/scoring/parameter-ranges")

    def rounds(self, limit: int = 100) -> list:
        """Most recent rounds, newest first (the API caps limit at 100)."""
        return self._get("/scoring/rounds", limit=min(limit, 100))["rounds"]

    def leaderboard(self, round_id: str) -> dict:
        return self._get(f"/scoring/rounds/{self._rid(round_id)}/leaderboard")

    def scores(self, round_id: str) -> dict:
        """Every validator score row with full metrics and submitted_at."""
        return self._get(f"/scoring/rounds/{self._rid(round_id)}/scores")

    def verification(self, round_id: str) -> dict:
        """Draw, file hashes and (after the window closes) the released files."""
        return self._get(f"/verification/round/{self._rid(round_id)}")
