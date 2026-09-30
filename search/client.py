"""EvalClient: lets search strategies read like direct calls while the executor runs asynchronously.

    results = await client.propose_all(configs)         # wait for a batch (e.g. a GA generation)
    results = await client.propose_all(configs, quorum=0.9)  # barrier that tolerates stragglers
    fut = client.propose(config); outcome = await fut   # single config

Many strategies run concurrently as asyncio tasks on one client, so while one waits at a barrier the others keep
the executor busy. Failed configs come back as Outcome(ok=False); they are never raised.
"""
from __future__ import annotations

import asyncio
import itertools
import logging
import math
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Tuple

from common.params import ParamRanges
from common.store import ResultStore
from common.types import FULL, AnswerSet, CallSet, Config, EvalJob, Round, ScoreResult
from evaluator.base import Executor, Scorer

log = logging.getLogger(__name__)


@dataclass
class Outcome:
    config: Config
    result: Optional[ScoreResult]
    error: Optional[str] = None
    tag: str = ""

    @property
    def ok(self) -> bool:
        return self.result is not None and self.result.score is not None

    @property
    def score(self) -> Optional[float]:
        return self.result.score if self.ok else None


class EvalClient:
    def __init__(self, round: Round, executor: Executor, scorer: Scorer, ranges: ParamRanges, store: ResultStore,
                 answer_set: AnswerSet, budget: Optional[int] = None):
        self.round, self.executor, self.scorer, self.ranges, self.store = round, executor, scorer, ranges, store
        self.answer_set = answer_set
        self.budget = budget
        self.spent = 0.0                                            # full-evaluation units sent to the executor
        self._ids = itertools.count()
        self._pending: Dict[str, Tuple[asyncio.Future, Config, str]] = {}
        self.callsets: Dict[Tuple[str, str], CallSet] = {}          # (config_key, fidelity) -> callset
        self.configs: Dict[str, Config] = {}
        self.latest: Dict[Tuple[str, str], ScoreResult] = {}         # scored against the current answer set
        self.on_result: List[Callable[[Outcome], None]] = []
        self.on_rescore: List[Callable[[], None]] = []

    # ---- strategy-facing API -------------------------------------------------------------------------------
    @property
    def exhausted(self) -> bool:
        return self.budget is not None and self.spent >= self.budget

    def propose(self, config: Config, *, fidelity: str = FULL, priority: int = 0, tag: str = "") -> asyncio.Future:
        return self.propose_many([config], fidelity=fidelity, priority=priority, tag=tag)[0]

    def propose_many(self, configs: Iterable[Config], *, fidelity: str = FULL, priority: int = 0,
                    tag: str = "") -> List[asyncio.Future]:
        """Several configs in ONE executor call, so related ones (a family: same config except cheap
        parameters) are queued together and can share one GATK run. Budget: each is charged the executor's
        cost for it as part of this group (a batch rider costs a fraction of a full run)."""
        loop = asyncio.get_running_loop()
        futs, jobs = [], []
        for config in configs:
            fut = loop.create_future()
            futs.append(fut)
            errs = self.ranges.errors(config)
            if errs:
                fut.set_result(Outcome(config, None, "invalid: " + "; ".join(errs), tag))
                continue
            cached = self.latest.get((config.key, fidelity))
            if cached is not None:
                fut.set_result(Outcome(config, cached, None, tag))
                continue
            jobs.append((fut, EvalJob(f"{self.round.tag}-j{next(self._ids)}", self.round, config, fidelity, priority, tag)))
        if not jobs:
            return futs
        hints = getattr(self.executor, "cost_hints", None)
        costs = hints([j for _, j in jobs]) if hints else [self.executor.cost_hint(j) for _, j in jobs]
        send = []
        for (fut, job), cost in zip(jobs, costs):
            if cost > 0:
                if self.exhausted:
                    fut.set_result(Outcome(job.config, None, "budget exhausted", tag))
                    continue
                self.spent += cost            # budget is in full-evaluation units
            self.configs[job.config.key] = job.config
            self._pending[job.job_id] = (fut, job.config, tag)
            send.append(job)
        if send:
            self.executor.submit(send)
        return futs

    def variant_cost(self, main: Config, variant: Config, fidelity: str = FULL) -> float:
        """Extra cost (fraction of a full run) of evaluating `variant` in the same run as `main`; 1.0 if they
        cannot share a run or the evaluator cannot tell."""
        runner = getattr(self.executor, "runner", None)
        f = getattr(runner, "variant_cost", None)
        c = f(main, variant, fidelity) if f else None
        return 1.0 if c is None else c

    async def propose_all(self, configs: Iterable[Config], *, quorum: float = 1.0, fidelity: str = FULL, priority: int = 0,
                  tag: str = "") -> List[Optional[Outcome]]:
        """Outcomes in input order. Returns once ceil(quorum * n) are done; the rest are None here but keep
        running, and their results still reach the store and every listener."""
        futs = self.propose_many(configs, fidelity=fidelity, priority=priority, tag=tag)
        need = math.ceil(quorum * len(futs))
        pending = set(futs)
        while len(futs) - len(pending) < need:
            _, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
        return [f.result() if f.done() else None for f in futs]

    def best(self, fidelity: str = FULL) -> Optional[ScoreResult]:
        ok = [r for (k, fid), r in self.latest.items()
              if fid == fidelity and r.score is not None and r.gate_passed is not False]
        return max(ok, key=lambda r: r.score) if ok else None

    # ---- driven by the round loop ----------------------------------------------------------------------------
    @property
    def in_flight(self) -> int:
        return len(self._pending)

    def deliver(self, cs: CallSet) -> None:
        """Score one callset now (unless its machine already did) and hand the outcome to the waiting strategy."""
        if cs.score is not None:
            self.complete(cs, cs.score)
            return
        result, error = (self.score_batch([cs])[0] if cs.ok else (None, None))
        self.complete(cs, result, error)

    def score_batch(self, callsets: List[CallSet]) -> List[Tuple[Optional[ScoreResult], Optional[str]]]:
        """(result, error) per callset against the current answer set. Safe to run in a worker thread; packs the
        batch into shared hap.py runs when the scorer supports it (score_many)."""
        many = getattr(self.scorer, "score_many", None)
        if many and len(callsets) > 1:
            try:
                return [(r, None) for r in many(callsets, self.answer_set)]
            except Exception as e:  # noqa: BLE001 - one bad callset must not sink the batch: fall back one by one
                log.warning("packed scoring failed (%s); scoring one by one", e)
        out = []
        for cs in callsets:
            try:
                out.append((self.scorer.score(cs, self.answer_set), None))
            except Exception as e:  # a scorer failure is that config's failure, not the round's
                log.warning("scoring failed for %s: %s", cs.config_key, e)
                out.append((None, f"scoring failed: {e}"))
        return out

    def complete(self, cs: CallSet, result: Optional[ScoreResult], error: Optional[str] = None) -> None:
        """Record a callset and its score (or failure), resolve the strategy's future, notify listeners."""
        if cs.round_id != self.round.round_id or cs.config_key not in self.configs:
            return                    # a late result of another round (a shared fleet executor outlives rounds)
        fut, config, tag = self._pending.pop(cs.job_id, (None, self.configs.get(cs.config_key), ""))
        outcome = Outcome(config, None, cs.error or error, tag)
        if cs.ok:
            self.callsets[(cs.config_key, cs.fidelity)] = cs
            if result is not None:
                self.latest[(cs.config_key, cs.fidelity)] = result
                self.store.add(result, config, cs.runtime_s, tag)
                outcome.result = result
        if fut is not None and not fut.done():
            fut.set_result(outcome)
        for cb in self.on_result:
            cb(outcome)

    def idle_slots(self) -> Optional[int]:
        """Executor capacity free right now (None if the executor cannot tell)."""
        f = getattr(self.executor, "idle_slots", None)
        return f() if f else None

    def rescore(self, answer_set: AnswerSet) -> None:
        """New answer-set version: re-score every stored callset that is on this machine (~0.6 s each; callsets
        scored on fleet nodes stay there and keep their score), then notify strategies. Failures keep the old score."""
        self.answer_set = answer_set
        for (key, fid), cs in self.callsets.items():
            if cs.score is not None:                     # scored on a fleet node; its VCF is not here
                continue
            try:
                r = self.scorer.score(cs, answer_set)
            except Exception as e:  # noqa: BLE001
                log.warning("rescoring %s failed: %s", key, e)
                continue
            self.latest[(key, fid)] = r
            self.store.add(r, self.configs[key], cs.runtime_s, "rescore")
        for cb in self.on_rescore:
            cb()
