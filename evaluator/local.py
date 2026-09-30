"""Executor on this machine's cores: priority queue, dedup of identical jobs, result cache, opportunistic batching."""
from __future__ import annotations

import heapq
import itertools
import threading
from collections import defaultdict
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import replace
from typing import Dict, Iterable, List, Optional, Tuple

from common.settings import default_workers
from common.types import CallSet, EvalJob
from evaluator.base import Runner

DIVERGED = "minos:diverged"   # same marker as evaluator.gatk.DIVERGED (kept here to avoid a GATK import)


def job_key(job: EvalJob) -> Tuple[str, str, str]:
    """Two jobs with the same key produce the same callset (runs are deterministic)."""
    return job.round.round_id, job.config.key, job.fidelity


def failed(job: EvalJob, error: str) -> CallSet:
    return CallSet(job.round.round_id, job.config.key, "", "", 0.0, ok=False, error=error,
                   job_id=job.job_id, fidelity=job.fidelity)


class LocalExecutor:
    """Runs jobs on local worker slots.

    Batching: if the runner offers `branch_key`/`run_batch`, a free slot takes the highest-priority job plus every
    queued job sharing its branch key (up to `max_batch`) and evaluates them in one run. Nothing ever waits to be
    batched - only jobs already queued ride along. A variant that comes back DIVERGED is re-queued to run alone.
    """

    def __init__(self, runner: Runner, workers: Optional[int] = None, max_batch: int = 32):
        self.runner = runner
        workers = workers or default_workers()
        self.workers = workers
        self.max_batch = max_batch
        self._batching = hasattr(runner, "run_batch") and hasattr(runner, "branch_key")
        self._pool = ThreadPoolExecutor(workers)
        self._queue: list = []                     # (-priority, seq, job)
        self._seq = itertools.count()
        self._running: Dict = {}                   # future -> [jobs]
        self._inflight: Dict = {}                  # key -> job actually queued/running
        self._waiters = defaultdict(list)          # key -> duplicate jobs waiting on the same run
        self._done: Dict = {}                      # key -> CallSet
        self._ready: List[CallSet] = []
        self._cancelled = set()
        self._solo = set()                         # job ids that must not be batched again
        self._lock = threading.Lock()
        self.batches_run = 0
        self.batched_jobs = 0

    # ---- running ---------------------------------------------------------------------------------------------
    def _branch(self, job: EvalJob) -> Optional[str]:
        if not self._batching or job.job_id in self._solo:
            return None
        return self.runner.branch_key(job.config, job.fidelity)

    def _run(self, jobs: List[EvalJob]) -> List[CallSet]:
        try:
            if len(jobs) == 1:
                out = [self.runner.run(jobs[0].round, jobs[0].config, jobs[0].fidelity)]
            else:
                out = self.runner.run_batch(jobs[0].round, [j.config for j in jobs], jobs[0].fidelity)
        except Exception as e:  # a crashed config is a result, not a crashed executor
            return [failed(j, str(e)) for j in jobs]
        return [replace(cs, job_id=j.job_id, fidelity=j.fidelity) for cs, j in zip(out, jobs)]

    def _drop_cancelled(self, job: EvalJob) -> None:
        k = job_key(job)
        waiters = self._waiters.pop(k, [])
        if waiters:            # a duplicate still wants this result: it takes the cancelled job's place
            first, rest = waiters[0], waiters[1:]
            self._inflight[k] = first
            self._waiters[k] = rest
            heapq.heappush(self._queue, (-first.priority, next(self._seq), first))
        else:
            self._inflight.pop(k, None)

    def _fill(self):
        while len(self._running) < self.workers and self._queue:
            _, _, job = heapq.heappop(self._queue)
            if job.job_id in self._cancelled:
                self._drop_cancelled(job)
                continue
            batch = [job]
            bk = self._branch(job)
            if bk is not None:
                keep = []
                while self._queue:
                    item = heapq.heappop(self._queue)
                    other = item[2]
                    if (len(batch) < self.max_batch and other.job_id not in self._cancelled
                            and other.round.round_id == job.round.round_id and self._branch(other) == bk):
                        batch.append(other)
                    else:
                        keep.append(item)
                for item in keep:
                    heapq.heappush(self._queue, item)
            if len(batch) > 1:
                self.batches_run += 1
                self.batched_jobs += len(batch)
            self._running[self._pool.submit(self._run, batch)] = batch

    def submit(self, jobs: Iterable[EvalJob]) -> None:
        with self._lock:
            for job in jobs:
                k = job_key(job)
                if k in self._done:
                    self._ready.append(replace(self._done[k], job_id=job.job_id))
                elif k in self._inflight:
                    self._waiters[k].append(job)
                else:
                    self._inflight[k] = job
                    heapq.heappush(self._queue, (-job.priority, next(self._seq), job))
            self._fill()

    def poll(self, timeout: float) -> List[CallSet]:
        with self._lock:
            running = list(self._running) if not self._ready and self._running else []
        if running:
            wait(running, timeout=timeout, return_when=FIRST_COMPLETED)
        with self._lock:
            for f in [f for f in self._running if f.done()]:
                jobs = self._running.pop(f)
                for job, cs in zip(jobs, f.result()):
                    k = job_key(job)
                    if not cs.ok and cs.error == DIVERGED:
                        # its batch output was not exact: run it alone (still counts as in flight)
                        self._solo.add(job.job_id)
                        heapq.heappush(self._queue, (-job.priority, next(self._seq), job))
                        continue
                    self._inflight.pop(k, None)
                    if cs.ok:
                        self._done[k] = cs
                    self._ready.append(cs)
                    self._ready += [replace(cs, job_id=w.job_id) for w in self._waiters.pop(k, [])]
            self._fill()
            out, self._ready = self._ready, []
            return out

    # ---- queries ---------------------------------------------------------------------------------------------
    def cancel(self, job_ids: Iterable[str]) -> None:
        with self._lock:
            self._cancelled.update(job_ids)

    def cost_hint(self, job: EvalJob) -> float:
        """0 if done or in flight; ~0.01 if a queued job would carry it in the same batch; else 1."""
        with self._lock:
            k = job_key(job)
            if k in self._done or k in self._inflight:
                return 0.0
            bk = self._branch(job)
            if bk is not None and any(self._branch(q[2]) == bk and q[2].round.round_id == job.round.round_id
                                      for q in self._queue):
                return 0.01
            return 1.0

    def idle_slots(self) -> int:
        with self._lock:
            return max(0, self.workers - len(self._running) - len(self._queue))
