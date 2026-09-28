"""Interfaces main depends on. Owner: evaluator.

Runner   - turns (round, config, fidelity) into a callset. Expensive, deterministic. Knows GATK.
Executor - decides where and when jobs run (local cores, fleet, simulation): queueing by priority,
           dedup of identical jobs, caching, retries. Never decides WHAT to run.
Scorer   - callset x answer set -> ScoreResult. Cheap, pure; called again for every stored callset
           whenever the answer set's version changes.
"""
from __future__ import annotations

from typing import Iterable, List, Protocol

from common.types import AnswerSet, CallSet, Config, EvalJob, Round, ScoreResult


class Runner(Protocol):
    def run(self, round: Round, config: Config, fidelity: str) -> CallSet:
        """Must match the validator's callset exactly. Raise or return CallSet(ok=False) on failure."""


class Executor(Protocol):
    def submit(self, jobs: Iterable[EvalJob]) -> None:
        """Queue jobs. Never blocks and never holds a job back to wait for others."""

    def poll(self, timeout: float) -> List[CallSet]:
        """Callsets completed since the last poll (job_id set). Blocks up to `timeout` if none are ready."""

    def cancel(self, job_ids: Iterable[str]) -> None:
        """Drop queued jobs. Jobs already running finish and are still returned (their result is data)."""

    def cost_hint(self, job: EvalJob) -> float:
        """Expected cost relative to one full evaluation: 0.0 if already done or running."""

    def idle_slots(self) -> int:
        """Free capacity right now (workers minus running minus queued)."""


class Scorer(Protocol):
    def score(self, callset: CallSet, answer_set: AnswerSet) -> ScoreResult:
        """Match callset vs answer set, then common.scoring.finalize. Must equal the validator's hap.py path."""
