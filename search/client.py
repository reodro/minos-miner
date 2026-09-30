"""EvalClient: lets search strategies read like direct calls while the executor runs asynchronously.

    results = await client.map(configs)                 # wait for a batch (e.g. a GA generation)
    results = await client.map(configs, quorum=0.9)     # barrier that tolerates stragglers
    fut = client.submit(config); outcome = await fut    # single job
    async for o in client.as_completed(futures): ...    # streaming (e.g. RF_UCB refill)

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

    def submit(self, config: Config, *, fidelity: str = FULL, priority: int = 0, tag: str = "") -> asyncio.Future:
        return self.submit_many([config], fidelity=fidelity, priority=priority, tag=tag)[0]

    def submit_many(self, configs: Iterable[Config], *, fidelity: str = FULL, priority: int = 0,
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
            jobs.append((fut, EvalJob(f"j{next(self._ids)}", self.round, config, fidelity, priority, tag)))
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

    async def map(self, configs: Iterable[Config], *, quorum: float = 1.0, fidelity: str = FULL, priority: int = 0,
                  tag: str = "") -> List[Optional[Outcome]]:
        """Outcomes in input order. Returns once ceil(quorum * n) are done; the rest are None here but keep
        running, and their results still reach the store and every listener."""
        futs = self.submit_many(configs, fidelity=fidelity, priority=priority, tag=tag)
        need = math.ceil(quorum * len(futs))
        pending = set(futs)
        while len(futs) - len(pending) < need:
            _, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
        return [f.result() if f.done() else None for f in futs]

    async def as_completed(self, futures: Iterable[asyncio.Future]):
        for f in asyncio.as_completed(list(futures)):
            yield await f

    def best(self, fidelity: str = FULL) -> Optional[ScoreResult]:
        ok = [r for (k, fid), r in self.latest.items()
              if fid == fidelity and r.score is not None and r.gate_passed is not False]
        return max(ok, key=lambda r: r.score) if ok else None

    # ---- driven by the round loop ----------------------------------------------------------------------------
    @property
    def in_flight(self) -> int:
        return len(self._pending)

    def deliver(self, cs: CallSet) -> None:
        fut, config, tag = self._pending.pop(cs.job_id, (None, self.configs.get(cs.config_key), ""))
        outcome = Outcome(config, None, cs.error, tag)
        if cs.ok:
            self.callsets[(cs.config_key, cs.fidelity)] = cs
            try:
                r = self.scorer.score(cs, self.answer_set)
                self.latest[(cs.config_key, cs.fidelity)] = r
                self.store.add(r, config, cs.runtime_s, tag)
                outcome.result = r
            except Exception as e:  # a scorer failure is that config's failure, not the round's
                log.warning("scoring failed for %s: %s", cs.config_key, e)
                outcome.error = f"scoring failed: {e}"
        if fut is not None and not fut.done():
            fut.set_result(outcome)
        for cb in self.on_result:
            cb(outcome)

    def rescore(self, answer_set: AnswerSet) -> None:
        """New answer-set version: re-score every stored callset (milliseconds each), then notify strategies."""
        self.answer_set = answer_set
        for (key, fid), cs in self.callsets.items():
            r = self.scorer.score(cs, answer_set)
            self.latest[(key, fid)] = r
            self.store.add(r, self.configs[key], cs.runtime_s, "rescore")
        for cb in self.on_rescore:
            cb()
