"""Scoring in worker processes on this machine. Owner: evaluator.

hap.py's Python half (upstream's parsers over its output) holds the GIL for seconds per result, so one interpreter
scores ~6 results a minute while 14 GATK workers deliver ~17: a round scored by one thread is scoring-bound. A fleet
node therefore scores in one process per GATK worker; ScoringPool is the same thing for main's own executor. Every
process shares the machine's per-round decision cache (a cluster scored by any of them is never sent to hap.py
again) and scores one GATK batch's callsets at a time.
"""
from __future__ import annotations

import logging
import os
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from typing import List, Optional, Tuple

from common.types import AnswerSet, CallSet, ScoreResult

_PROC_SCORER = None


def _proc_init(work_dir: str) -> None:
    """A scorer process: its own interpreter (own GIL), scoring with upstream's hap.py path. It exits as soon as its
    parent (the node agent, or any other caller) is gone - even after a hard kill - so no scorer runs on its own."""
    global _PROC_SCORER
    parent = os.getppid()

    def watchdog():
        while os.getppid() == parent:
            time.sleep(2)
        os._exit(1)

    threading.Thread(target=watchdog, daemon=True).start()
    logging.basicConfig(level=logging.INFO, stream=sys.stderr, format="%(asctime)s %(levelname)s scorer: %(message)s")
    from evaluator.sandbox import install_scoring_shim
    install_scoring_shim()                     # this process imports scoring itself; the parent's shim does not carry over
    from evaluator.decision_cache import DecisionCacheScorer
    _PROC_SCORER = DecisionCacheScorer(work_dir=work_dir)    # every scorer process shares work_dir's per-round store


def _proc_score(cs: CallSet, aset: AnswerSet) -> tuple:
    """Score one callset in a scorer process: (callset, ScoreResult), or (failed callset, None) carrying the error."""
    try:
        return cs, _PROC_SCORER.score(cs, aset)
    except Exception as e:  # noqa: BLE001 - a callset hap.py cannot score is a failed result, never a 0.0
        return replace(cs, ok=False, error=f"scoring failed on node: {e}"), None


def _proc_score_unit(css: list, aset: AnswerSet) -> list:
    """Score one GATK batch's callsets in a scorer process, one after another: their clusters not yet cached on this
    node are packed into as few hap.py runs as possible (score_many). If that fails, each is scored on its own, so one
    bad callset cannot fail the others."""
    if len(css) == 1:
        return [_proc_score(css[0], aset)]
    try:
        return list(zip(css, _PROC_SCORER.score_many(css, aset)))
    except Exception:  # noqa: BLE001
        return [_proc_score(c, aset) for c in css]



def scorer_processes(workers: int) -> int:
    """Scorer processes for a machine that runs this many GATK workers.

    One hap.py process scores about six callsets a minute. A worker finishes less than one walk a
    minute. One scorer per six workers matches that, with at least two on any machine that has more
    than one worker, and never more scorers than workers. One per worker puts a second runnable on
    every core.
    """
    n = max(1, int(workers))
    if n == 1:
        return 1
    return min(n, max(2, n // 6))


class ScoringPool:
    def __init__(self, parallel: int, score_dir: Optional[str] = None, scorer=None):
        """`score_dir` set: `parallel` scorer processes (upstream's hap.py path, the decision cache in score_dir).
        `scorer` set (tests): that one scorer, called from up to `parallel` threads."""
        self.parallel = max(1, parallel)
        self._scorer = scorer
        self._pool = None
        if scorer is None:
            import multiprocessing
            self._pool = ProcessPoolExecutor(self.parallel, mp_context=multiprocessing.get_context("spawn"),
                                             initializer=_proc_init, initargs=(score_dir,))

    def score_unit(self, callsets: List[CallSet], answer_set: AnswerSet) -> List[Tuple[CallSet, Optional[ScoreResult]]]:
        """(callset, score) per callset; one that could not be scored comes back ok=False with the error and no score.
        Blocks until done: call it from a thread, up to `parallel` at once."""
        if self._scorer is not None:
            out = []
            for cs in callsets:
                try:
                    out.append((cs, self._scorer.score(cs, answer_set)))
                except Exception as e:  # noqa: BLE001 - a callset hap.py cannot score is a failed result, never a 0.0
                    out.append((replace(cs, ok=False, error=f"scoring failed: {e}"), None))
            return out
        return self._pool.submit(_proc_score_unit, callsets, answer_set).result()

    def close(self) -> None:
        if self._pool is not None:
            self._pool.shutdown(wait=False, cancel_futures=True)
