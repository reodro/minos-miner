"""Executor on a virtual clock: simulated runtimes, crashes and stragglers, no waiting.

Lets main test multi-strategy search at thousands of evals per second, with a measured core utilisation.
Same semantics as LocalExecutor (priority queue, dedup, cache).
"""
from __future__ import annotations

import heapq
import itertools
from collections import defaultdict
from dataclasses import replace
from typing import Callable, Dict, Iterable, List

from common.types import CallSet, EvalJob
from evaluator.base import Runner
from evaluator.local import failed, job_key


class SimExecutor:
    def __init__(self, runner: Runner, workers: int = 16,
                 duration: Callable[[EvalJob], float] = lambda job: 60.0,
                 crashes: Callable[[EvalJob], bool] = lambda job: False):
        self.runner, self.workers, self.duration, self.crashes = runner, workers, duration, crashes
        self.now = 0.0
        self.busy = 0.0
        self._queue: list = []
        self._seq = itertools.count()
        self._running: list = []                  # (finish_time, seq, job)
        self._inflight: Dict = {}
        self._waiters = defaultdict(list)
        self._done: Dict = {}
        self._ready: List[CallSet] = []
        self._cancelled = set()

    def _fill(self):
        while len(self._running) < self.workers and self._queue:
            _, _, job = heapq.heappop(self._queue)
            if job.job_id in self._cancelled:
                k = job_key(job)
                waiters = self._waiters.pop(k, [])
                if waiters:
                    self._inflight[k], self._waiters[k] = waiters[0], waiters[1:]
                    heapq.heappush(self._queue, (-waiters[0].priority, next(self._seq), waiters[0]))
                else:
                    self._inflight.pop(k, None)
                continue
            d = self.duration(job)
            self.busy += d
            heapq.heappush(self._running, (self.now + d, next(self._seq), job))

    def submit(self, jobs: Iterable[EvalJob]) -> None:
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

    def poll(self, timeout: float = 0.0) -> List[CallSet]:
        if not self._ready and self._running:
            t = self._running[0][0]
            self.now = t
            while self._running and self._running[0][0] == t:
                _, _, job = heapq.heappop(self._running)
                k = job_key(job)
                self._inflight.pop(k, None)
                if self.crashes(job):
                    cs = failed(job, "simulated crash")
                else:
                    cs = replace(self.runner.run(job.round, job.config, job.fidelity), job_id=job.job_id,
                                 fidelity=job.fidelity)
                    self._done[k] = cs
                self._ready.append(cs)
                self._ready += [replace(cs, job_id=w.job_id) for w in self._waiters.pop(k, [])]
            self._fill()
        out, self._ready = self._ready, []
        return out

    def cancel(self, job_ids: Iterable[str]) -> None:
        self._cancelled.update(job_ids)

    def cost_hint(self, job: EvalJob) -> float:
        k = job_key(job)
        return 0.0 if k in self._done or k in self._inflight else 1.0

    def idle_slots(self) -> int:
        return max(0, self.workers - len(self._running) - len(self._queue))

    def utilization(self) -> float:
        """Fraction of worker-seconds spent running jobs so far."""
        if self.now == 0:
            return 0.0
        not_yet_run = sum(f - self.now for f, _, _ in self._running)
        return (self.busy - not_yet_run) / (self.workers * self.now)
