"""One live (or replayed) round: start the search at once, build truth meanwhile, keep the best config submitted."""
from __future__ import annotations

import asyncio
import concurrent.futures
import logging
import time
from dataclasses import dataclass, field
from typing import Awaitable, Callable, Dict, List, Optional, Tuple

from answer_set.base import AnswerSetProvider
from common.params import ParamRanges
from common.store import ResultStore
from common.types import FULL, AnswerSet, CallSet, Config, Round, ScoreResult
from evaluator.base import Executor, Scorer
from search.client import EvalClient, Outcome
from main.submitter import Submitter

log = logging.getLogger(__name__)


class AnswerSetUnavailable(Exception):
    """The round has no usable answer set (the builder failed, or the caller refused what it built): nothing can be
    scored, so the round falls back to a baseline. The executor's work for the round has been dropped."""


@dataclass
class LoopSettings:
    update_every: int = 50              # new callsets between answer_set.update() calls
    resubmit_margin: float = 0.0        # score points a new best must gain before resubmitting
    resubmit_interval_s: float = 0.0    # min seconds between submissions (chain commitments are rate-limited);
                                        # a better config found meanwhile is submitted when the interval ends
                                        # and, at the latest, just before the loop stops
    stop_before_deadline_s: float = 120.0
    build_timeout_s: float = 1500.0     # a build not done by then (a hung 1000G path) ends the round on the baseline
    score_margin_s: float = 45.0        # scoring after GATK: work must end this long before the final submission
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
    failures: Dict[str, int] = field(default_factory=dict)                     # failure reason -> count
    failed_configs: Dict[str, dict] = field(default_factory=dict)              # config_key -> {"params", "reason"}:
                                                                               # what to exclude from the space later


async def run_round_async(round: Round, provider: AnswerSetProvider, executor: Executor, scorer: Scorer,
                          searcher: Callable[[EvalClient], Awaitable[None]], submitter: Submitter, ranges: ParamRanges,
                          settings: LoopSettings, baselines: List[Config], store: Optional[ResultStore] = None,
                          deadline_ts: Optional[float] = None,
                          clock: Optional[Callable[[], float]] = None,
                          on_tick: Optional[Callable] = None,
                          accept: Optional[Callable[[AnswerSet], Optional[str]]] = None,
                          scoring=None) -> RoundOutcome:
    """The search starts at once and the answer set is built meanwhile (minutes on a live round): the executor runs
    the baselines and the opening wave, and their results are held until the truth exists, then scored. `accept`
    sees every answer set the provider builds and may refuse it with a reason: the round then ends with
    AnswerSetUnavailable (the caller submits a baseline). `scoring`: an evaluator.ScoringPool for the callsets this
    machine scores (processes, as many at once as it has); without one they are scored one unit at a time."""
    virtual = hasattr(executor, "now")
    clock = clock or ((lambda: executor.now) if virtual else time.time)
    client = EvalClient(round, executor, scorer, ranges, store or ResultStore(), None,
                        budget=settings.max_evals, clock=clock,
                        deadline_ts=None if deadline_ts is None else deadline_ts - settings.stop_before_deadline_s)
    out = RoundOutcome()
    submitted: List[Optional[str]] = [None]
    last_submit = [float("-inf")]
    last_refused = [float("-inf")]                # a refused submission (platform down, fee) is retried after 30 s
    since_update: List[CallSet] = []
    held: List[CallSet] = []                      # GATK results that arrived before the answer set
    started = clock()
    defaults = ranges.defaults(round_tool(baselines)).params if baselines else {}

    def changed(config: Config) -> dict:
        return {k: v for k, v in config.params.items() if v != defaults.get(k)}

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
            reason = _failure_reason(o.error)
            out.failures[reason] = out.failures.get(reason, 0) + 1
            if o.config is not None and len(out.failed_configs) < 500:
                out.failed_configs[o.config.key] = {"params": changed(o.config), "reason": reason}
            if out.failures[reason] in (1, 10, 100, 1000):
                log.warning("config %s failed (%d so far with this reason): %s | differs from the defaults in %s",
                            o.config.key if o.config else "?", out.failures[reason], reason,
                            changed(o.config) if o.config else "?")
        maybe_submit()

    client.on_result.append(on_result)
    client.on_rescore.append(maybe_submit)

    def time_left() -> bool:
        return deadline_ts is None or clock() < deadline_ts - settings.stop_before_deadline_s

    def publish_answer_set(aset: AnswerSet) -> None:
        publish = getattr(executor, "publish_answer_set", None)   # fleet nodes score with it from now on
        if publish:
            publish(round.round_id, aset)

    async def build():
        """The answer set, while the executor already runs the first wave (the virtual executor builds inline: its
        clock only moves in poll, so the result is the same every run)."""
        # its own thread, not the loop's default executor: asyncio.run waits for the default executor at exit, and a
        # build that hangs must not hold the round's end (the thread ends on its own, unheard)
        pool = concurrent.futures.ThreadPoolExecutor(1, thread_name_prefix="answer-set")
        try:
            aset = provider.build(round) if virtual else await asyncio.get_running_loop().run_in_executor(pool, provider.build, round)
        except Exception as e:  # noqa: BLE001 - the round still gets a baseline
            log.exception("answer set failed for %s", round.round_id)
            raise AnswerSetUnavailable(f"answer set failed: {e}") from e
        finally:
            pool.shutdown(wait=False)
        why = accept(aset) if accept else None
        if why:
            raise AnswerSetUnavailable(why)
        client.set_answer_set(aset)
        out.answer_set_version = aset.version
        publish_answer_set(aset)
        log.info("%s: answer set v%d ready at t+%.0fs; %d results waited for it", round.round_id, aset.version,
                 clock() - started, len(held))

    async def search():
        # Baselines (highest priority, so they run first) and the search start together: the fleet fills at once
        # instead of idling while the baseline runs. Baselines are tried in order until one has been submitted (or a
        # better config has), so the round never ends empty.
        async def baselines_in_order():
            for c in baselines:
                await client.propose(c, priority=10, tag="baseline")
                if submitted[0] is not None:
                    break

        await asyncio.gather(baselines_in_order(), searcher(client))

    # The executor starts no GATK work that could not be scored before the final submission (it would only be killed).
    finish_by = getattr(executor, "set_finish_by", None)
    if finish_by and deadline_ts is not None:
        finish_by(deadline_ts - settings.stop_before_deadline_s - settings.score_margin_s)
    building = asyncio.create_task(build())
    task = asyncio.create_task(search())
    # Live executors: scoring (hap.py, seconds per callset) runs off the loop - in the pool's processes, as many
    # units at once as it has - so the loop keeps polling GATK and the search keeps proposing meanwhile. The
    # virtual executor scores inline (deterministic).
    to_score: List[CallSet] = []
    scoring_tasks: set = set()
    stopping = [False]
    unit_size = 8                                 # callsets per scoring unit (about one GATK batch's worth)
    parallel = scoring.parallel if scoring is not None else 1

    async def score_unit(unit: List[CallSet]):
        try:
            if scoring is not None:
                pairs = await asyncio.to_thread(scoring.score_unit, unit, client.answer_set)
                for cs, r in pairs:
                    client.complete(cs, r, None if r is not None else (cs.error or "scoring failed"))
            else:
                results = await asyncio.to_thread(client.score_batch, unit)
                for cs, (res, err) in zip(unit, results):
                    client.complete(cs, res, err)
        except Exception as e:  # noqa: BLE001 - a broken scorer fails these callsets, never the round
            log.warning("scoring %d callsets failed: %s", len(unit), e)
            for cs in unit:
                client.complete(cs, None, f"scoring failed: {e}")
        finally:
            scoring_tasks.discard(asyncio.current_task())
            start_scoring()

    def start_scoring():
        while to_score and len(scoring_tasks) < parallel and not stopping[0]:
            unit = to_score[:unit_size]
            del to_score[:unit_size]
            scoring_tasks.add(asyncio.create_task(score_unit(unit)))

    idle_spins = 0
    try:
        while time_left():
            for _ in range(20):                  # let the search react before time moves on
                await asyncio.sleep(0)
            if building.done() and building.exception() is not None:
                raise building.exception()       # no truth for this round: the caller falls back to a baseline
            if not building.done() and clock() - started > settings.build_timeout_s:
                raise AnswerSetUnavailable(f"answer set not ready after {settings.build_timeout_s:.0f} s")
            if task.done() and building.done() and client.in_flight == 0:
                break
            batch = executor.poll(0) if virtual else await asyncio.to_thread(executor.poll, settings.poll_timeout)
            if client.answer_set is None:        # GATK ran ahead of the truth: successes wait, failures need none
                held += [cs for cs in batch if cs.ok]
                batch = [cs for cs in batch if not cs.ok]
            elif held:
                batch, held = held + batch, []
            for cs in batch:
                if virtual or not cs.ok or cs.score is not None:     # already scored on its node: no work here
                    client.deliver(cs)
                else:
                    to_score.append(cs)
            if not virtual:
                start_scoring()
            if (client.answer_set is not None and len(since_update) >= settings.update_every
                    and not scoring_tasks):                           # never rescore mid-unit
                new = provider.update(client.answer_set, list(since_update))
                since_update.clear()
                if new.version != client.answer_set.version:
                    client.rescore(new)
                    out.answer_set_version = new.version
                    publish_answer_set(new)
            maybe_submit()                        # an improvement held back by the interval goes out when it ends
            if on_tick:
                on_tick(client, out, executor)    # the dashboard snapshot (throttled by the callee)
            idle_spins = idle_spins + 1 if (virtual and not batch and client.in_flight == 0) else 0
            if idle_spins > 1000:                 # the search alive but proposing nothing: a search bug
                log.warning("search stalled with nothing in flight; ending round")
                break
    finally:
        stopping[0] = True                        # finish the units in flight, start no new one
        if scoring_tasks:                         # results already computed are worth submitting on
            await asyncio.gather(*scoring_tasks, return_exceptions=True)
        maybe_submit(force=True)                  # last chance: the best found so far, whatever the interval
        end = getattr(executor, "end_round", None)  # nobody will read this round's results any more: free every
        if end:                                   # worker (queued jobs dropped, running containers killed)
            end(round.round_id)
        if finish_by:                             # a fleet executor serves the next round too
            finish_by(None)
        elif client.in_flight:
            executor.cancel(list(client._pending))
        for t in (task, building):
            if not t.done():
                t.cancel()
                if t is building:                 # a build in a thread cannot be cancelled: it ends on its own, unheard
                    t.add_done_callback(lambda b: b.cancelled() or b.exception())
                    continue
            try:
                await t
            except (asyncio.CancelledError, Exception) as e:
                if not isinstance(e, (asyncio.CancelledError, AnswerSetUnavailable)):
                    log.warning("%s task ended with %r", "search" if t is task else "answer set", e)
    if hasattr(executor, "utilization"):
        out.utilization = executor.utilization()
    out.callsets = {k: cs for (k, fid), cs in client.callsets.items() if fid == FULL}
    return out


def round_tool(baselines: List[Config]) -> str:
    return baselines[0].tool if baselines else "gatk"


def _failure_reason(error: Optional[str]) -> str:
    """One short line per kind of failure (the reason text's first line, without run-specific ids)."""
    import re
    text = (error or "unknown").strip().splitlines()[0]
    text = re.sub(r"0x[0-9a-f]+|[0-9a-f]{12,}|/tmp/\S+", "*", text)
    return text[:120]


def run_round(*args, **kwargs) -> RoundOutcome:
    return asyncio.run(run_round_async(*args, **kwargs))
