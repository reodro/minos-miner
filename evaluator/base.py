"""The interface main depends on. Owner: evaluator."""
from __future__ import annotations

from typing import Protocol

from common.types import AnswerSet, CallSet, Config, Round, ScoreResult


class Evaluator(Protocol):
    def run(self, round: Round, config: Config) -> CallSet:
        """Expensive (~100 s full region): produce the callset exactly as a validator would.

        Deterministic: the same (round, config) must give the same callset, so results
        are cached by (round_id, config.key). A failed run returns CallSet(ok=False)
        and is never retried.
        """

    def score(self, callset: CallSet, answer_set: AnswerSet) -> ScoreResult:
        """Cheap (milliseconds): match callset vs answer set, then common.scoring.finalize.

        Called again for every stored callset whenever the answer set's version changes.
        With uncertain synthetic candidates, score is the expectation over their p.
        """
