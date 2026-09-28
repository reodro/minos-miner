"""One live (or replayed) round: build truth, search, keep the best config submitted."""
from __future__ import annotations

import logging
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

from answer_set.base import AnswerSetProvider
from common.params import ParamRanges
from common.types import AnswerSet, CallSet, Config, Round, ScoreResult
from evaluator.base import Evaluator
from main.search import Searcher
from main.submitter import Submitter

log = logging.getLogger(__name__)


@dataclass
class LoopSettings:
    workers: int = 4
    update_every: int = 50              # new callsets between answer_set.update() calls
    resubmit_margin: float = 0.0        # score points a new best must gain before resubmitting
    stop_before_deadline_s: float = 120.0
    max_evals: Optional[int] = None     # budget by count (backtests must match config count, not wall time)


@dataclass
class RoundOutcome:
    answer_set_version: int
    evaluated: int
    failed: int
    best: Optional[ScoreResult]
    submissions: List[Tuple[float, str, float]] = field(default_factory=list)   # (t, config_key, score)


def run_round(round: Round, provider: AnswerSetProvider, evaluator: Evaluator, searcher: Searcher,
              submitter: Submitter, ranges: ParamRanges, settings: LoopSettings,
              baselines: List[Config], deadline_ts: Optional[float] = None,
              clock: Callable[[], float] = time.time) -> RoundOutcome:
    answer_set: AnswerSet = provider.build(round)
    configs: Dict[str, Config] = {}
    callsets: Dict[str, CallSet] = {}
    results: Dict[str, ScoreResult] = {}
    outcome = RoundOutcome(answer_set.version, 0, 0, None)
    submitted_key: Optional[str] = None
    since_update: List[CallSet] = []

    def time_left() -> bool:
        return deadline_ts is None or clock() < deadline_ts - settings.stop_before_deadline_s

    def budget_left() -> bool:
        return settings.max_evals is None or outcome.evaluated + outcome.failed < settings.max_evals

    def maybe_submit():
        nonlocal submitted_key
        ok = [r for r in results.values() if r.score is not None and r.gate_passed]
        if not ok:
            return
        best = max(ok, key=lambda r: r.score)
        outcome.best = best
        if best.config_key == submitted_key:
            return
        # Compare against the submitted config's CURRENT score: an answer-set upgrade rescored it.
        current = results[submitted_key].score if submitted_key else None
        if current is None or best.score > current + settings.resubmit_margin:
            if submitter.submit(round, configs[best.config_key]):
                submitted_key = best.config_key
                outcome.submissions.append((clock(), best.config_key, best.score))

    def accept(cs: CallSet):
        nonlocal answer_set
        callsets[cs.config_key] = cs
        if not cs.ok:
            outcome.failed += 1
            return
        outcome.evaluated += 1
        r = evaluator.score(cs, answer_set)
        results[cs.config_key] = r
        searcher.observe([r])
        since_update.append(cs)
        if len(since_update) >= settings.update_every:
            new = provider.update(answer_set, list(since_update))
            since_update.clear()
            if new.version != answer_set.version:
                answer_set = new
                outcome.answer_set_version = new.version
                for k, c in callsets.items():
                    if c.ok:
                        results[k] = evaluator.score(c, answer_set)
                searcher.rescored(list(results.values()))
        maybe_submit()

    def valid_new(cands: List[Config]) -> List[Config]:
        out = []
        for c in cands:
            errs = ranges.errors(c)
            if errs:
                log.warning("skipping invalid config %s: %s", c.key, errs)
            elif c.key not in configs:
                configs[c.key] = c
                out.append(c)
        return out

    def run_safe(c: Config) -> CallSet:
        try:
            return evaluator.run(round, c)
        except Exception as e:  # a crashed run is a failed config, never a crashed round
            log.warning("run failed for %s: %s", c.key, e)
            return CallSet(round.round_id, c.key, "", "", 0.0, ok=False, error=str(e))

    # Baselines first and in order: the first one that runs is submitted immediately,
    # so the round never ends without a submission.
    for c in valid_new(baselines):
        accept(run_safe(c))
        if submitted_key is not None or not time_left():
            break

    with ThreadPoolExecutor(max_workers=settings.workers) as pool:
        pending = {}
        while time_left() and budget_left():
            room = settings.workers - len(pending)
            if settings.max_evals is not None:
                room = min(room, settings.max_evals - outcome.evaluated - outcome.failed - len(pending))
            for c in valid_new(searcher.propose(room)) if room > 0 else []:
                pending[pool.submit(run_safe, c)] = c
            if not pending:
                break
            done, _ = wait(pending, timeout=5.0, return_when=FIRST_COMPLETED)
            for f in done:
                pending.pop(f)
                accept(f.result())
        for f in pending:
            f.cancel()
    return outcome
