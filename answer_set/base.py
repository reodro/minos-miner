"""The interface main depends on. Owner: answer_set."""
from __future__ import annotations

from typing import Iterable, Protocol

from common.types import AnswerSet, CallSet, Round


class AnswerSetProvider(Protocol):
    def build(self, round: Round) -> AnswerSet:
        """Called once at round start. Must return within a few minutes; version=1."""

    def update(self, answer_set: AnswerSet, new_callsets: Iterable[CallSet]) -> AnswerSet:
        """Called repeatedly during the round with callsets evaluated so far.

        Returns the same object if nothing changed, or a new AnswerSet with version+1.
        Must use independent evidence (reads, other callers, priors) to label new
        sites - never a vote of the configs being ranked (legacy circularity trap).
        """
