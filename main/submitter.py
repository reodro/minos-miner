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


class PlatformSubmitter:
    """Real submission through the upstream miner's platform client. Not implemented yet.

    Must: sign with the hotkey, send config.to_submission(), never send a fee
    (resubmission_fee_enabled is false), and treat a 402/4xx as not accepted.
    """

    def submit(self, round: Round, config: Config) -> bool:
        raise NotImplementedError("main.PlatformSubmitter.submit")
