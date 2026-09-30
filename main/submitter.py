"""Submitting configs to the platform. The last submission before the deadline is the one scored."""
from __future__ import annotations

import logging
from typing import List, Protocol, Tuple

from common.types import Config, Round

log = logging.getLogger(__name__)


class Submitter(Protocol):
    def submit(self, round: Round, config: Config) -> bool:
        """True if the platform accepted it."""


class DryRunSubmitter:
    """Records submissions instead of sending them. Used in backtests and tests."""

    def __init__(self):
        self.sent: List[Tuple[str, str]] = []

    def submit(self, round: Round, config: Config) -> bool:
        self.sent.append((round.round_id, config.key))
        log.info("dry-run submit %s %s", round.round_id, config.key)
        return True


# The live submitter is main.miner.PlatformSubmitter (over main.miner.UpstreamPlatform).
