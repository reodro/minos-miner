"""One live (or replayed) round: build truth, run strategies concurrently, keep the best config submitted."""
from __future__ import annotations

import asyncio
import logging
import time
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple

from answer_set.base import AnswerSetProvider
from common.params import ParamRanges
from common.store import ResultStore
from common.types import FULL, CallSet, Config, Round, ScoreResult
from evaluator.base import Executor, Scorer
from search.client import EvalClient, Outcome
from search.base import Strategy
from main.submitter import Submitter

log = logging.getLogger(__name__)


@dataclass
class LoopSettings:
    update_every: int = 50              # new callsets between answer_set.update() calls
    resubmit_margin: float = 0.0        # score points a new best must gain before resubmitting
    stop_before_deadline_s: float = 120.0
    max_evals: Optional[int] = None     # budget by count (backtests must match config count, not wall time)
    poll_timeout: float = 1.0


@dataclass
class RoundOutcome:
    answer_set_version: int = 1
    evaluated: int = 0
    failed: int = 0
    best: Optional[ScoreResult] = None
    submissions: List[Tuple[float, str, float]] = field(default_factory=list)   # (t, config_key, score)
    utilization: Optional[float] = None                                        # SimExecutor only


async def run_round_async(round: Round, provider: AnswerSetProvider, executor: Executor, scorer: Scorer,
                          strategies: List[Strategy], submitter: Submitter, ranges: ParamRanges,
                          settings: LoopSettings, baselines: List[Config], store: Optional[ResultStore] = None,
                          deadline_ts: Optional[float] = None,
                          clock: Optional[Callable[[], float]] = None) -> RoundOutcome:
    virtual = hasattr(executor, "now")
    clock = clock or ((lambda: executor.now) if virtual else time.time)
    client = EvalClient(round, executor, scorer, ranges, store or ResultStore(), provider.build(round),
                        budget=settings.max_evals)
    out = RoundOutcome(client.answer_set.version)
    submitted: List[Optional[str]] = [None]
    since_update: List[CallSet] = []

    def maybe_submit():
        best = client.best(FULL)
        if best is None:
            return
        out.best = best
        key = submitted[0]
        if best.config_key == key:
            return
        current = client.latest[(key, FULL)].score if key else None
        if current is None or best.score > current + settings.resubmit_margin:
            if submitter.submit(round, client.configs[best.config_key]):
                submitted[0] = best.config_key
                out.submissions.append((clock(), best.config_key, best.score))

    def on_result(o: Outcome):
        if o.ok:
            out.evaluated += 1
            since_update.append(client.callsets[(o.config.key, o.result.fidelity)])
        else:
            out.failed += 1
        maybe_submit()

    client.on_result.append(on_result)
    client.on_rescore.append(maybe_submit)

    def time_left() -> bool:
        return deadline_ts is None or clock() < deadline_ts - settings.stop_before_deadline_s

    async def search():
        # Baselines first, in order: the first one that runs is submitted, so the round never ends empty.
        for c in baselines:
            await client.submit(c, priority=10, tag="baseline")
            if submitted[0] is not None:
                break
        await asyncio.gather(*(s(client) for s in strategies))

    task = asyncio.create_task(search())
    idle_spins = 0
    try:
        while time_left():
            for _ in range(20):                  # let every strategy react before time moves on
                await asyncio.sleep(0)
            if task.done() and client.in_flight == 0:
                break
            batch = executor.poll(0) if virtual else await asyncio.to_thread(executor.poll, settings.poll_timeout)
            for cs in batch:
                client.deliver(cs)
            if len(since_update) >= settings.update_every:
                new = provider.update(client.answer_set, list(since_update))
                since_update.clear()
                if new.version != client.answer_set.version:
                    client.rescore(new)
                    out.answer_set_version = new.version
            idle_spins = idle_spins + 1 if (virtual and not batch and client.in_flight == 0) else 0
            if idle_spins > 1000:                 # strategies alive but submitting nothing: a strategy bug
                log.warning("strategies stalled with nothing in flight; ending round")
                break
    finally:
        task.cancel()
        try:
            await task
        except (asyncio.CancelledError, Exception) as e:
            if not isinstance(e, asyncio.CancelledError):
                log.warning("search task ended with %r", e)
    if hasattr(executor, "utilization"):
        out.utilization = executor.utilization()
    return out


def run_round(*args, **kwargs) -> RoundOutcome:
    return asyncio.run(run_round_async(*args, **kwargs))
