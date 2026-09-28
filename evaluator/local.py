"""Executor on this machine's cores: priority queue, dedup of identical jobs, in-memory result cache."""
from __future__ import annotations

import heapq
import itertools
import threading
from collections import defaultdict
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import replace
from typing import Dict, Iterable, List, Tuple

from common.types import CallSet, EvalJob
from evaluator.base import Runner


def job_key(job: EvalJob) -> Tuple[str, str, str]:
    """Two jobs with the same key produce the same callset (runs are deterministic)."""
    return job.round.round_id, job.config.key, job.fidelity


def failed(job: EvalJob, error: str) -> CallSet:
    return CallSet(job.round.round_id, job.config.key, "", "", 0.0, ok=False, error=error,
                   job_id=job.job_id, fidelity=job.fidelity)


class LocalExecutor:
    def __init__(self, runner: Runner, workers: int = 4):
        self.runner = runner
        self.workers = workers
        self._pool = ThreadPoolExecutor(workers)
        self._queue: list = []                     # (-priority, seq, job)
        self._seq = itertools.count()
        self._running: Dict = {}                   # future -> job
        self._inflight: Dict = {}                  # key -> job actually running/queued
        self._waiters = defaultdict(list)          # key -> duplicate jobs waiting on the same run
        self._done: Dict = {}                      # key -> CallSet
        self._ready: List[CallSet] = []
        self._cancelled = set()
        self._lock = threading.Lock()

    def _run(self, job: EvalJob) -> CallSet:
        try:
            cs = self.runner.run(job.round, job.config, job.fidelity)
        except Exception as e:  # a crashed config is a result, not a crashed executor
            return failed(job, str(e))
        return replace(cs, job_id=job.job_id, fidelity=job.fidelity)

    def _fill(self):
        while len(self._running) < self.workers and self._queue:
            _, _, job = heapq.heappop(self._queue)
            if job.job_id in self._cancelled:
                k = job_key(job)
                waiters = self._waiters.pop(k, [])
                if waiters:            # a duplicate still wants this result: it takes the cancelled job's place
                    first, rest = waiters[0], waiters[1:]
                    self._inflight[k] = first
                    self._waiters[k] = rest
                    heapq.heappush(self._queue, (-first.priority, next(self._seq), first))
                else:
                    self._inflight.pop(k, None)
                continue
            self._running[self._pool.submit(self._run, job)] = job

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
            if not self._ready and self._running:
                running = list(self._running)
            else:
                running = []
        if running:
            wait(running, timeout=timeout, return_when=FIRST_COMPLETED)
        with self._lock:
            for f in [f for f in self._running if f.done()]:
                job = self._running.pop(f)
                cs = f.result()
                k = job_key(job)
                self._inflight.pop(k, None)
                if cs.ok:
                    self._done[k] = cs
                self._ready.append(cs)
                self._ready += [replace(cs, job_id=w.job_id) for w in self._waiters.pop(k, [])]
            self._fill()
            out, self._ready = self._ready, []
            return out

    def cancel(self, job_ids: Iterable[str]) -> None:
        with self._lock:
            self._cancelled.update(job_ids)

    def cost_hint(self, job: EvalJob) -> float:
        k = job_key(job)
        return 0.0 if k in self._done or k in self._inflight else 1.0

    def idle_slots(self) -> int:
        with self._lock:
            return max(0, self.workers - len(self._running) - len(self._queue))
