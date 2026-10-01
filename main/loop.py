"""One live (or replayed) round: build truth, run strategies concurrently, keep the best config submitted."""
from __future__ import annotations

import asyncio
import logging
import time
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

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
    resubmit_interval_s: float = 0.0    # min seconds between submissions (chain commitments are rate-limited);
                                        # a better config found meanwhile is submitted when the interval ends
                                        # and, at the latest, just before the loop stops
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
    callsets: Dict[str, CallSet] = field(default_factory=dict)                 # config_key -> full-fidelity callset


async def run_round_async(round: Round, provider: AnswerSetProvider, executor: Executor, scorer: Scorer,
                          strategies: List[Strategy], submitter: Submitter, ranges: ParamRanges,
                          settings: LoopSettings, baselines: List[Config], store: Optional[ResultStore] = None,
                          deadline_ts: Optional[float] = None,
                          clock: Optional[Callable[[], float]] = None,
                          on_tick: Optional[Callable] = None) -> RoundOutcome:
    virtual = hasattr(executor, "now")
    clock = clock or ((lambda: executor.now) if virtual else time.time)
    client = EvalClient(round, executor, scorer, ranges, store or ResultStore(), provider.build(round),
                        budget=settings.max_evals)
    out = RoundOutcome(client.answer_set.version)
    submitted: List[Optional[str]] = [None]
    last_submit = [float("-inf")]
    last_refused = [float("-inf")]                # a refused submission (platform down, fee) is retried after 30 s
    since_update: List[CallSet] = []

    def maybe_submit(force: bool = False):
        best = client.best(FULL)
        if best is None:
            return
        out.best = best
        key = submitted[0]
        if best.config_key == key:
            return
        if key is not None and not force and clock() - last_submit[0] < settings.resubmit_interval_s:
            return                            # the first submission is never delayed; later ones wait
        if not force and clock() - last_refused[0] < 30.0:
            return
        current = client.latest[(key, FULL)].score if key else None
        if current is None or best.score > current + settings.resubmit_margin:
            if submitter.submit(round, client.configs[best.config_key]):
                submitted[0] = best.config_key
                last_submit[0] = clock()
                out.submissions.append((clock(), best.config_key, best.score))
            else:
                last_refused[0] = clock()

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
            await client.propose(c, priority=10, tag="baseline")
            if submitted[0] is not None:
                break
        await asyncio.gather(*(s(client) for s in strategies))

    task = asyncio.create_task(search())
    # Live executors: scoring (one hap.py run per callset, seconds each) runs in a worker thread, so the loop keeps
    # polling GATK and strategies keep proposing meanwhile. The virtual executor scores inline (deterministic).
    to_score: List[CallSet] = []
    scoring: List[Optional[asyncio.Task]] = [None]
    stopping = [False]

    def start_scoring():
        if scoring[0] is None and to_score and not stopping[0]:
            batch = list(to_score)
            to_score.clear()

            async def run():
                results = await asyncio.to_thread(client.score_batch, batch)
                for cs, (res, err) in zip(batch, results):
                    client.complete(cs, res, err)
                scoring[0] = None
                start_scoring()
            scoring[0] = asyncio.create_task(run())

    idle_spins = 0
    try:
        while time_left():
            for _ in range(20):                  # let every strategy react before time moves on
                await asyncio.sleep(0)
            if task.done() and client.in_flight == 0:
                break
            batch = executor.poll(0) if virtual else await asyncio.to_thread(executor.poll, settings.poll_timeout)
            for cs in batch:
                if virtual or not cs.ok or cs.score is not None:     # already scored on its node: no work here
                    client.deliver(cs)
                else:
                    to_score.append(cs)
            if not virtual:
                start_scoring()
            if len(since_update) >= settings.update_every and scoring[0] is None:   # never rescore mid-batch
                new = provider.update(client.answer_set, list(since_update))
                since_update.clear()
                if new.version != client.answer_set.version:
                    client.rescore(new)
                    out.answer_set_version = new.version
            maybe_submit()                        # an improvement held back by the interval goes out when it ends
            if on_tick:
                on_tick(client, out, executor)    # the dashboard snapshot (throttled by the callee)
            idle_spins = idle_spins + 1 if (virtual and not batch and client.in_flight == 0) else 0
            if idle_spins > 1000:                 # strategies alive but proposing nothing: a strategy bug
                log.warning("strategies stalled with nothing in flight; ending round")
                break
    finally:
        stopping[0] = True                        # finish the batch in flight, start no new one
        if scoring[0] is not None:                # results already computed are worth submitting on
            try:
                await scoring[0]
            except Exception:  # noqa: BLE001
                pass
        maybe_submit(force=True)                  # last chance: the best found so far, whatever the interval
        end = getattr(executor, "end_round", None)  # nobody will read this round's results any more: free every
        if end:                                   # worker (queued jobs dropped, running containers killed)
            end(round.round_id)
        elif client.in_flight:
            executor.cancel(list(client._pending))
        task.cancel()
        try:
            await task
        except (asyncio.CancelledError, Exception) as e:
            if not isinstance(e, asyncio.CancelledError):
                log.warning("search task ended with %r", e)
    if hasattr(executor, "utilization"):
        out.utilization = executor.utilization()
    out.callsets = {k: cs for (k, fid), cs in client.callsets.items() if fid == FULL}
    return out


def run_round(*args, **kwargs) -> RoundOutcome:
    return asyncio.run(run_round_async(*args, **kwargs))
