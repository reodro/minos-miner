"""Executor on this machine's cores: priority queue, dedup of identical jobs, result cache, opportunistic batching."""
from __future__ import annotations

import heapq
import itertools
import threading
import time
from collections import defaultdict
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import replace
from typing import Dict, Iterable, List, Optional, Tuple

from common.settings import default_workers
from common.types import CallSet, EvalJob
from evaluator.base import Runner

DIVERGED = "minos:diverged"   # same marker as evaluator.gatk.DIVERGED (kept here to avoid a GATK import)
SOLO_S = 132.0                # one full GATK run at production load
REPAIR_S = 34.0               # a diverged variant's repair run
WARM_S = SOLO_S * 0.011       # genotyping one variant in a JVM that already walked the family
MAX_BATCH_S = 4.5 * 60        # split a family that would run longer than this and leave the rest queued


def extra_seconds(runner, main: EvalJob, variant: EvalJob) -> Optional[float]:
    """Seconds `variant` adds inside `main`'s batch. None = it cannot join that batch.

    A runner without variant_cost adds nothing: the old behavior, batch up to max_batch.
    """
    cost = getattr(runner, "variant_cost", None)
    if cost is None:
        return 0.0
    frac = cost(main.config, variant.config, main.fidelity)
    if frac is None:
        return None
    return SOLO_S * float(frac)


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
    batched - only jobs already queued ride along. A variant that comes back DIVERGED is re-queued alone and fixed by
    runner.repair (or simply run again when the runner cannot repair).
    """

    def __init__(self, runner: Runner, workers: Optional[int] = None, max_batch: int = 32, now=None):
        self.runner = runner
        workers = workers or default_workers()
        self.workers = workers
        self.max_batch = max_batch
        self._now = now or time.time
        self._finish_by: Optional[float] = None   # do not start a batch that would still be running then
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
        self._repair: Dict[str, CallSet] = {}      # job id -> its diverged batch output, for runner.repair
        self._ended = set()                        # round ids whose results are no longer wanted (end_round)
        self._lock = threading.Lock()
        self._wake = threading.Event()             # set by submit(): an idle poll() waits on it instead of spinning
        self.batches_run = 0
        self.repairs_run = 0
        self.gatk_runs = 0                              # GATK runs actually started (solo, batch or repair)
        self.gatk_runs_done = 0                         # ... and finished
        self.done = 0                                   # callsets delivered (dashboard)

    # ---- running ---------------------------------------------------------------------------------------------
    def _served_warm(self, job: EvalJob) -> bool:
        """True when a live JVM can genotype this config without another traversal."""
        ready = getattr(self.runner, "warm_ready", None)
        if ready is None or job.job_id in self._repair or job.job_id in self._solo:
            return False
        return bool(ready(job.round, job.config, job.fidelity))

    def _evaluate(self, jobs: List[EvalJob]) -> List[CallSet]:
        """Genotype on a live JVM when one already walked this family; otherwise one traversal."""
        warm = getattr(self.runner, "warm_submit", None)
        if warm is not None and not any(j.job_id in self._solo for j in jobs):
            hits = [warm(j.round, j.config, j.fidelity) for j in jobs]
            if any(h is not None for h in hits):
                missed = [j for j, h in zip(jobs, hits) if h is None]
                by_id = {}
                if len(missed) == 1:
                    j = missed[0]
                    by_id[j.job_id] = self.runner.run(j.round, j.config, j.fidelity)
                elif missed:
                    batched = self.runner.run_batch(missed[0].round, [j.config for j in missed], missed[0].fidelity)
                    for j, cs in zip(missed, batched):
                        by_id[j.job_id] = cs
                return [h if h is not None else by_id[j.job_id] for j, h in zip(jobs, hits)]
        if len(jobs) == 1:
            j = jobs[0]
            return [self.runner.run(j.round, j.config, j.fidelity)]
        return self.runner.run_batch(jobs[0].round, [j.config for j in jobs], jobs[0].fidelity)

    def _branch(self, job: EvalJob) -> Optional[str]:
        if not self._batching or job.job_id in self._solo:
            return None
        return self.runner.branch_key(job.config, job.fidelity)

    def _run(self, jobs: List[EvalJob]) -> List[CallSet]:
        self.gatk_runs += 1
        try:
            if len(jobs) == 1 and jobs[0].job_id in self._repair:
                j = jobs[0]
                self.repairs_run += 1
                out = [self.runner.repair(j.round, j.config, self._repair.pop(j.job_id), j.fidelity)]
            else:
                out = self._evaluate(jobs)
        except Exception as e:  # a crashed config is a result, not a crashed executor
            return [failed(j, str(e)) for j in jobs]
        finally:
            self.gatk_runs_done += 1
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

    def set_finish_by(self, ts: Optional[float]) -> None:
        """Do not start a batch that would still be running at `ts` (the final submission). None clears it."""
        with self._lock:
            self._finish_by = ts

    def _room(self) -> float:
        """Seconds a batch started now may take: the family cap, and the time left before the final submission."""
        room = MAX_BATCH_S
        if self._finish_by is not None:
            room = min(room, self._finish_by - self._now())
        return room

    def _fill(self):
        skipped = []
        while len(self._running) < self.workers and self._queue:
            _, _, job = heapq.heappop(self._queue)
            if job.job_id in self._cancelled:
                self._drop_cancelled(job)
                continue
            room = self._room()
            served = self._served_warm(job)
            main_s = WARM_S if served else (REPAIR_S if job.job_id in self._repair else SOLO_S)
            if main_s > room:
                # it would be killed at the deadline; leave it queued and try a shorter job behind it
                skipped.append(job)
                continue
            batch, seconds = [job], main_s
            bk = self._branch(job)
            if bk is not None:
                keep = []
                while self._queue:
                    item = heapq.heappop(self._queue)
                    other = item[2]
                    extra = None
                    if (len(batch) < self.max_batch and other.job_id not in self._cancelled
                            and other.round.round_id == job.round.round_id and self._branch(other) == bk
                            and (not served or self._served_warm(other))):
                        extra = extra_seconds(self.runner, job, other)
                    if extra is not None and seconds + extra <= room:
                        batch.append(other)
                        seconds += extra
                    else:
                        keep.append(item)
                for item in keep:
                    heapq.heappush(self._queue, item)
            if len(batch) > 1 and not served:
                self.batches_run += 1
            self._running[self._pool.submit(self._run, batch)] = batch
        for job in skipped:
            heapq.heappush(self._queue, (-job.priority, next(self._seq), job))

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
            self._wake.set()

    def poll(self, timeout: float) -> List[CallSet]:
        with self._lock:
            running = list(self._running) if not self._ready and self._running else []
            idle = not self._ready and not self._running
            if idle:
                self._wake.clear()
        if running:
            wait(running, timeout=timeout, return_when=FIRST_COMPLETED)
        elif idle:
            # nothing to wait for: block until a submit (or the timeout) instead of returning at once - callers poll
            # in loops, and a spinning loop starves every other thread of the GIL (seen: node scoring 20x slower)
            self._wake.wait(timeout)
        with self._lock:
            for f in [f for f in self._running if f.done()]:
                jobs = self._running.pop(f)
                for job, cs in zip(jobs, f.result()):
                    k = job_key(job)
                    if job.round.round_id in self._ended:        # its round is over: nobody wants this result
                        self._inflight.pop(k, None)
                        self._waiters.pop(k, None)
                        self._repair.pop(job.job_id, None)
                        continue
                    if not cs.ok and cs.error == DIVERGED:
                        # its batch output was not exact: repair it if the runner can, else run it alone
                        # (still counts as in flight)
                        self._solo.add(job.job_id)
                        if hasattr(self.runner, "repair") and cs.vcf:
                            self._repair[job.job_id] = cs
                        heapq.heappush(self._queue, (-job.priority, next(self._seq), job))
                        continue
                    self._inflight.pop(k, None)
                    if cs.ok:
                        self._done[k] = cs
                    self._ready.append(cs)
                    self._ready += [replace(cs, job_id=w.job_id) for w in self._waiters.pop(k, [])]
            self._fill()
            out, self._ready = self._ready, []
            self.done += len(out)
            return out

    # ---- queries ---------------------------------------------------------------------------------------------
    def cancel(self, job_ids: Iterable[str]) -> None:
        with self._lock:
            self._cancelled.update(job_ids)

    def stats(self) -> dict:
        with self._lock:
            return {"kind": "local", "workers": self.workers, "running": len(self._running), "queued": len(self._queue),
                    "batches_run": self.batches_run, "repairs_run": self.repairs_run, "done": self.done,
                    "gatk_runs": self.gatk_runs, "gatk_runs_done": self.gatk_runs_done}

    def end_round(self, round_id: str) -> None:
        """The round is over: drop its queued jobs, kill its running GATK containers, deliver nothing more for it.
        Frees every worker for the next round at once instead of letting stale batches run on for minutes."""
        with self._lock:
            self._ended.add(round_id)
            keep = [item for item in self._queue if item[2].round.round_id != round_id]
            for item in self._queue:
                if item[2].round.round_id == round_id:
                    k = job_key(item[2])
                    self._inflight.pop(k, None)
                    self._waiters.pop(k, None)
                    self._repair.pop(item[2].job_id, None)
            self._queue = keep
            heapq.heapify(self._queue)
        kill = getattr(self.runner, "kill_round", None)
        if kill:
            kill(round_id)

    def cost_hint(self, job: EvalJob) -> float:
        """0 if done or in flight; if a queued job would carry it in its batch, the runner's variant cost
        (a fraction of a solo run, by the most expensive GATK stage it changes); else 1."""
        with self._lock:
            k = job_key(job)
            if k in self._done or k in self._inflight:
                return 0.0
            bk = self._branch(job)
            if bk is None:
                return 1.0
            costs = [self.runner.variant_cost(q[2].config, job.config, job.fidelity) for q in self._queue
                     if q[2].round.round_id == job.round.round_id and self._branch(q[2]) == bk]
            costs = [c for c in costs if c is not None]
            return min(costs, default=1.0)

    def cost_hints(self, jobs: List[EvalJob]) -> List[float]:
        """cost_hint for jobs submitted together: a job that can ride in the batch of one listed before it
        costs only its variant cost, as if that one were already queued."""
        out = []
        for i, job in enumerate(jobs):
            c = self.cost_hint(job)
            bk = self._branch(job)
            if c > 0 and bk is not None and hasattr(self.runner, "variant_cost"):
                for prev in jobs[:i]:
                    if prev.round.round_id == job.round.round_id and self._branch(prev) == bk:
                        v = self.runner.variant_cost(prev.config, job.config, job.fidelity)
                        if v is not None:
                            c = min(c, v)
            out.append(c)
        return out

    def capacity(self) -> int:
        return self.workers

    def idle_slots(self) -> int:
        with self._lock:
            return max(0, self.workers - len(self._running) - len(self._queue))
