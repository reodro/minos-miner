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
                 crashes: Callable[[EvalJob], bool] = lambda job: False, batching: bool = True):
        """batching: pack queued configs that can share a run, as the real executors do (False: one job per run,
        the older strategies' tests were written against that)."""
        self.runner, self.workers, self.duration, self.crashes = runner, workers, duration, crashes
        self.batching = batching
        self.max_batch = 32
        self.batches = 0                          # runs started (a batch is one run)
        self.now = 0.0
        self.busy = 0.0
        self._queue: list = []
        self._seq = itertools.count()
        self._running: list = []                  # (finish_time, seq, batch: list of jobs)
        self._inflight: Dict = {}
        self._waiters = defaultdict(list)
        self._done: Dict = {}
        self._ready: List[CallSet] = []
        self._cancelled = set()

    def _branch(self, job: EvalJob):
        f = getattr(self.runner, "branch_key", None)
        return f(job.config, job.fidelity) if f else None

    def _drop_cancelled(self, job: EvalJob) -> None:
        k = job_key(job)
        waiters = self._waiters.pop(k, [])
        if waiters:
            self._inflight[k], self._waiters[k] = waiters[0], waiters[1:]
            heapq.heappush(self._queue, (-waiters[0].priority, next(self._seq), waiters[0]))
        else:
            self._inflight.pop(k, None)

    def _fill(self):
        """Like LocalExecutor: a free slot takes the top job plus queued jobs that can share its run (same branch
        key, up to max_batch); the batch takes the head's duration plus each rider's cost fraction of it, and all
        its results arrive together. A latency-lane head takes only cheap riders; a latency job never rides."""
        while len(self._running) < self.workers and self._queue:
            _, _, job = heapq.heappop(self._queue)
            if job.job_id in self._cancelled:
                self._drop_cancelled(job)
                continue
            batch, d = [job], self.duration(job)
            bk = self._branch(job) if self.batching else None
            cost = getattr(self.runner, "variant_cost", None)
            if bk is not None and cost is not None:
                cap = 0.02 if job.lane == "latency" else 1.0
                keep = []
                while self._queue and len(batch) < self.max_batch:
                    item = heapq.heappop(self._queue)
                    other = item[2]
                    frac = None
                    if (other.job_id not in self._cancelled and other.lane != "latency"
                            and other.round.round_id == job.round.round_id and self._branch(other) == bk):
                        frac = cost(job.config, other.config, job.fidelity)
                    if frac is not None and frac <= cap:
                        batch.append(other)
                        d += frac * self.duration(job)
                    else:
                        keep.append(item)
                for item in keep:
                    heapq.heappush(self._queue, item)
            self.busy += d
            self.batches += 1
            heapq.heappush(self._running, (self.now + d, next(self._seq), batch))

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
                _, _, batch = heapq.heappop(self._running)
                for job in batch:
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
        """Like LocalExecutor: 0 if done or in flight; a queued job that would carry it prices it at the runner's
        variant cost; else 1 (a solo run)."""
        k = job_key(job)
        if k in self._done or k in self._inflight:
            return 0.0
        bk = self._branch(job) if self.batching else None
        cost = getattr(self.runner, "variant_cost", None)
        if bk is None or cost is None:
            return 1.0
        costs = [cost(q[2].config, job.config, job.fidelity) for q in self._queue
                 if q[2].round.round_id == job.round.round_id and self._branch(q[2]) == bk]
        return min((c for c in costs if c is not None), default=1.0)

    def cost_hints(self, jobs: List[EvalJob]) -> List[float]:
        """cost_hint for jobs submitted together: a job that can ride in the batch of one listed before it costs
        only its variant cost."""
        out = []
        cost = getattr(self.runner, "variant_cost", None)
        for i, job in enumerate(jobs):
            c = self.cost_hint(job)
            bk = self._branch(job) if self.batching else None
            if c > 0 and bk is not None and cost is not None:
                for prev in jobs[:i]:
                    if prev.round.round_id == job.round.round_id and self._branch(prev) == bk:
                        v = cost(prev.config, job.config, job.fidelity)
                        if v is not None:
                            c = min(c, v)
            out.append(c)
        return out

    def capacity(self) -> int:
        return self.workers

    def idle_slots(self) -> int:
        queued_runs = (len({self._branch(j) or j.job_id for _, _, j in self._queue}) if self.batching
                       else len(self._queue)) if self._queue else 0
        return max(0, self.workers - len(self._running) - queued_runs)

    def utilization(self) -> float:
        """Fraction of worker-seconds spent running jobs so far."""
        if self.now == 0:
            return 0.0
        not_yet_run = sum(f - self.now for f, _, _ in self._running)
        return (self.busy - not_yet_run) / (self.workers * self.now)
