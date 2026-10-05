"""Executor on this machine's cores: priority queue, dedup of identical jobs, result cache, opportunistic batching,
and provably identical configs run once (evaluator/equivalence.py, the patch's haplotype-cap report)."""
from __future__ import annotations

import heapq
import itertools
import logging
import math
import os
import threading
import time
from collections import Counter, defaultdict, deque
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import replace
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from common.settings import default_workers
from common.types import CallSet, Config, EvalJob, Round
from evaluator.base import Runner
from evaluator.equivalence import RoundFacts, canonical, scan

DIVERGED = "minos:diverged"   # same marker as evaluator.gatk.DIVERGED (kept here to avoid a GATK import)
SOLO_S = 132.0                # one full GATK run at production load
REPAIR_S = 34.0               # a diverged variant's repair run
WARM_S = SOLO_S * 0.011       # genotyping one variant in a JVM that already walked the family
SMALL_BATCH = 1.3             # a small batch holds at most this much work, in solo runs: its head and cheap siblings
BIG_BATCH = 2.3               # a big batch at most this much: packed, yet never a long serial tail on one worker
MAX_BATCH_S = BIG_BATCH * SOLO_S   # no batch runs longer: a family beyond it is split and the rest stays queued
MNH = "max_num_haplotypes_in_population"
CLOCK_WINDOW = 30             # recent runs the deadline estimate is taken from
CLOCK_QUANTILE = 0.9          # ... and how pessimistic it is (a run this slow or slower is 1 in 10)
TIMELINE_S = 60.0             # one executor sample per minute (stats()["timeline"] and the log)

log = logging.getLogger(__name__)


def quantile_scale(ratios, q: float = CLOCK_QUANTILE) -> float:
    """Real seconds per table second: the q-quantile of the measured ratios, 1.0 (the table) before any."""
    if not ratios:
        return 1.0
    s = sorted(ratios)
    return s[max(0, math.ceil(q * len(s)) - 1)]


def _slow_tier(frac: float) -> Optional[str]:
    """Assembly or activity, and only when `frac` is that stage's table fraction. PairHMM, genotyping, and a
    fake cost stay on the 132 s table so a measured solo cannot split those families."""
    from evaluator.gatk import VARIANT_COST
    for name in ("activity", "assembly"):
        if abs(float(frac) - VARIANT_COST[name]) < 1e-6:
            return name
    return None


def extra_seconds(runner, main: EvalJob, variant: EvalJob, learned: Optional[Dict[str, float]] = None) -> Optional[float]:
    """Seconds `variant` adds inside `main`'s batch. None = it cannot join that batch.

    A runner without variant_cost adds nothing: the old behavior, batch up to max_batch.
    Genotyping and PairHMM stay at SOLO_S times their fraction. A learned assembly or activity extra replaces
    that table price on this machine; callers that pass no `learned` map (the fleet coordinator) keep the table.
    """
    cost = getattr(runner, "variant_cost", None)
    if cost is None:
        return 0.0
    frac = cost(main.config, variant.config, main.fidelity)
    if frac is None:
        return None
    tier = _slow_tier(float(frac))
    if learned and tier and tier in learned:
        return float(learned[tier])
    return SOLO_S * float(frac)


def _link_as(vcf: str, ran_key: str, asked_key: str) -> None:
    """The VCF of a config that ran as another (provably identical) one, also under its own name
    (<key>_<fidelity>.vcf.gz, + .tbi), as hard links: whatever finds callsets by name (main.review) still does."""
    src = Path(vcf) if vcf else None
    if src is None or not src.name.startswith(ran_key) or not src.exists():
        return
    for suffix in ("", ".tbi"):
        a, b = Path(str(src) + suffix), src.with_name(asked_key + src.name[len(ran_key):] + suffix)
        if a.exists() and not b.exists():
            try:
                os.link(a, b)
            except OSError:
                pass


def _without_cap(config: Config) -> str:
    """The config's key with max_num_haplotypes_in_population left out: configs differing only in that cap."""
    return Config(config.tool, {k: v for k, v in config.params.items() if k != MNH}).key


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
        self._dropped = set()                      # job ids whose results are no longer wanted (cancel_running)
        self._retried = set()                      # job ids whose batch died in GATK's own code (head rerun alone)
        self._lock = threading.Lock()
        self._wake = threading.Event()             # set by submit(): an idle poll() waits on it instead of spinning
        self.batches_run = 0
        self.repairs_run = 0
        self.retries = 0                                # jobs re-queued after an rc=3 crash
        self.gatk_runs = 0                              # GATK runs actually started (solo, batch or repair)
        self.gatk_runs_done = 0                         # ... and finished
        self.done = 0                                   # callsets delivered (dashboard)
        self._ratios: deque = deque(maxlen=CLOCK_WINDOW)   # real wall / table seconds of recent cold runs
        self._learned_extra: Dict[str, float] = {}     # assembly / activity seconds measured on this machine
        self.warm_runs = 0                              # runs that genotyped on a live JVM (no traversal)
        self._held = 0                                  # queued jobs the last fill left waiting for the deadline
        self._submitted = 0                             # jobs submitted since the last timeline sample
        self._last_submit: Optional[float] = None
        self._timeline: deque = deque(maxlen=150)       # one sample per minute: where the workers' time went
        self._next_sample = 0.0
        # Provably identical configs run once (evaluator/equivalence.py and the patch's haplotype-cap report): a job
        # runs as its canonical config and its result is delivered under the config it asked for.
        self._facts: Dict[str, Optional[RoundFacts]] = {}   # round id -> its BAM's facts (None while scanning)
        self._as_key: Dict[str, str] = {}                # job id -> the config key it asked for, when it runs as another
        self._caps: Dict[tuple, List[Tuple[int, Config]]] = defaultdict(list)  # (round, fidelity, config w/o cap) ->
        #                                                  (cap needed, the config that ran): any cap >= needed is that run
        self.canonicalized = 0                           # jobs that ran (or reused a result) as a canonical config

    # ---- running ---------------------------------------------------------------------------------------------
    def _served_warm(self, job: EvalJob) -> bool:
        """True when a live JVM can genotype this config without another traversal."""
        ready = getattr(self.runner, "warm_ready", None)
        if ready is None or job.job_id in self._repair or job.job_id in self._solo:
            return False
        return bool(ready(job.round, job.config, job.fidelity))

    def _evaluate(self, jobs: List[EvalJob]) -> Tuple[List[CallSet], bool]:
        """Genotype on a live JVM when one already walked this family; otherwise one traversal.
        Also says whether a live JVM served any of them (then the wall is not a traversal's)."""
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
                return [h if h is not None else by_id[j.job_id] for j, h in zip(jobs, hits)], True
        if len(jobs) == 1:
            j = jobs[0]
            return [self.runner.run(j.round, j.config, j.fidelity)], False
        return self.runner.run_batch(jobs[0].round, [j.config for j in jobs], jobs[0].fidelity), False

    def _branch(self, job: EvalJob) -> Optional[str]:
        if not self._batching or job.job_id in self._solo:
            return None
        return self.runner.branch_key(job.config, job.fidelity)

    def _run(self, jobs: List[EvalJob], table_s: Optional[float] = None) -> List[CallSet]:
        """One GATK run. `table_s` is what the fill priced it at (table seconds); None for a run not to learn from."""
        self.gatk_runs += 1
        t0 = self._now()
        repair = len(jobs) == 1 and jobs[0].job_id in self._repair
        warm = False
        try:
            if repair:
                j = jobs[0]
                self.repairs_run += 1
                out = [self.runner.repair(j.round, j.config, self._repair.pop(j.job_id), j.fidelity)]
            else:
                out, warm = self._evaluate(jobs)
        except Exception as e:  # a crashed config is a result, not a crashed executor
            return [failed(j, str(e)) for j in jobs]
        finally:
            self.gatk_runs_done += 1
        if warm:
            self.warm_runs += 1
        elif not repair:
            self._learn(jobs, out, self._now() - t0, table_s)
        return [replace(cs, job_id=j.job_id, fidelity=j.fidelity) for cs, j in zip(out, jobs)]

    def _scale(self) -> float:
        """Real seconds per table second on this machine (recent runs, pessimistic quantile; 1.0 before any).
        Caller holds the lock."""
        return quantile_scale(self._ratios)

    def _deadline_s(self) -> float:
        """Seconds a cold solo walk takes here, by the same estimate the deadline uses (the fleet node reports it)."""
        return SOLO_S * self._scale()

    def _learn(self, jobs: List[EvalJob], results: List[CallSet], wall: float, table_s: Optional[float]) -> None:
        """Learn this machine's speed from a finished cold run, and an assembly or activity extra that ran longer
        than the table.

        Speed: wall / the table seconds the fill priced the run at, kept for the last CLOCK_WINDOW runs. The deadline
        uses a high quantile of it, so a run (solo or family) is started only when it will most likely finish.
        Packing itself stays in table seconds: scaling the extras would push a six-copy PairHMM family over the
        family cap and split it. A repair, a failure, a live-JVM run, or a run under half a table solo (a fake or
        cached run) is not a lesson.
        """
        if not results or len(results) != len(jobs) or any(not cs.ok for cs in results):
            return
        with self._lock:
            if table_s and table_s > 0 and wall >= SOLO_S * 0.5:
                self._ratios.append(wall / table_s)
            if len(jobs) == 1:
                return
            cost = getattr(self.runner, "variant_cost", None)
            if cost is None:
                return
            fracs = []
            for other in jobs[1:]:
                frac = cost(jobs[0].config, other.config, jobs[0].fidelity)
                if frac is None:
                    return
                fracs.append(float(frac))
            top = max(fracs)
            tier = _slow_tier(top)
            if tier is None or tier in self._learned_extra:
                return
            implied = (wall - SOLO_S) / (len(jobs) - 1)
            if implied > SOLO_S * top:
                self._learned_extra[tier] = implied

    def _route(self, job: EvalJob) -> None:
        """A (canonical) job: reuse a finished result, wait on the same run, or queue it. Caller holds the lock."""
        k = job_key(job)
        if k in self._done:
            self._ready.append(self._restore(job.job_id, self._done[k]))
        elif k in self._inflight:
            self._waiters[k].append(job)
        else:
            self._inflight[k] = job
            heapq.heappush(self._queue, (-job.priority, next(self._seq), job))

    # ---- provably identical configs ----------------------------------------------------------------------------
    def _as_canonical(self, job: EvalJob) -> EvalJob:
        """The job as the config it will run: its canonical config on this round (BAM facts), or the config of a
        finished run its haplotype cap provably reproduces. Remembers the config it asked for. Caller holds the lock."""
        config = self._canonical_config(job)
        if config.key == job.config.key:
            return job
        self._as_key.setdefault(job.job_id, job.config.key)
        self.canonicalized += 1
        return replace(job, config=config)

    def _canonical_config(self, job: EvalJob) -> Config:
        c = canonical(job.config, self._facts_for(job.round))
        cap = c.params.get(MNH)
        if isinstance(cap, int) and not isinstance(cap, bool):
            for needed, ran in self._caps.get((job.round.round_id, job.fidelity, _without_cap(c)), ()):
                if cap >= needed:
                    return ran
        return c

    def _facts_for(self, round: Round) -> Optional[RoundFacts]:
        """This round's BAM facts; the first call starts the scan (a few seconds) in the background and returns None
        until it is done - jobs submitted meanwhile simply run as they are. Caller holds the lock."""
        rid = round.round_id
        if rid not in self._facts:
            self._facts[rid] = None
            if round.bam and os.path.exists(round.bam):
                threading.Thread(target=self._scan, args=(rid, round.bam), daemon=True).start()
        return self._facts[rid]

    def _scan(self, round_id: str, bam: str) -> None:
        try:
            facts = scan(bam)
        except Exception as e:  # noqa: BLE001 - no facts: configs run as they are
            log.warning("equivalence: cannot scan %s: %s", bam, e)
            return
        with self._lock:
            if round_id in self._facts:
                self._facts[round_id] = facts
        log.info("equivalence %s: %d reads, max %d per start, %d base qualities", round_id, facts.reads,
                 facts.max_reads_per_start, len(facts.quals))

    def _learn_cap(self, job: EvalJob, cs: CallSet) -> None:
        """After an exact run: the haplotype cap that reproduces it (the patch's report). Every config identical but
        for a cap at least that large is this run - including jobs already queued. Caller holds the lock."""
        needed = getattr(self.runner, "kbest_needed", None)
        cap = job.config.params.get(MNH)
        if needed is None or not isinstance(cap, int) or isinstance(cap, bool):
            return
        n = needed(job.round.round_id, job.config.key, job.fidelity)
        if n is None or n > cap:
            return
        cls = (job.round.round_id, job.fidelity, _without_cap(job.config))
        if any(ran.key == job.config.key for _, ran in self._caps[cls]):
            return
        self._caps[cls].append((n, job.config))
        moved, keep = [], []
        for item in self._queue:
            q = item[2]
            qcap = q.config.params.get(MNH)
            if (q.round.round_id, q.fidelity) == cls[:2] and isinstance(qcap, int) and qcap >= n \
                    and q.config.key != job.config.key and _without_cap(q.config) == cls[2]:
                moved.append(q)
            else:
                keep.append(item)
        if not moved:
            return
        self._queue = keep
        heapq.heapify(self._queue)
        for q in moved:                                  # it and every duplicate waiting on it become this run
            k = job_key(q)
            self._inflight.pop(k, None)
            for w in [q] + self._waiters.pop(k, []):
                self._as_key.setdefault(w.job_id, w.config.key)
                self.canonicalized += 1
                self._route(replace(w, config=job.config))
        log.info("equivalence: %d queued configs reuse %s (any max_num_haplotypes >= %d)", len(moved),
                 job.config.key, n)

    def _restore(self, job_id: str, cs: CallSet) -> CallSet:
        """A result for this job, under the config key it asked for. Caller holds the lock."""
        key = self._as_key.pop(job_id, None)
        cs = replace(cs, job_id=job_id)
        if key is not None and key != cs.config_key:
            _link_as(cs.vcf, cs.config_key, key)
            cs = replace(cs, config_key=key, score=None if cs.score is None else replace(cs.score, config_key=key))
        return cs

    def _drop_cancelled(self, job: EvalJob) -> None:
        self._as_key.pop(job.job_id, None)
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
        """Table seconds a batch started now may take: the family cap, and the time left before the final
        submission divided by this machine's measured speed (a 250 s walk on a 132 s table is 1.9x).
        Caller holds the lock."""
        room = MAX_BATCH_S
        if self._finish_by is not None:
            room = min(room, (self._finish_by - self._now()) / self._scale())
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
            repairing = job.job_id in self._repair
            main_s = WARM_S if served else (REPAIR_S if repairing else SOLO_S)
            if main_s > room:
                # it would be killed at the deadline; leave it queued and try a shorter job behind it
                skipped.append(job)
                continue
            batch, seconds = [job], main_s
            bk = self._branch(job)
            limit_s = (SMALL_BATCH if job.batch == "small" else BIG_BATCH) * SOLO_S
            if bk is not None:
                keep = []
                while self._queue:
                    item = heapq.heappop(self._queue)
                    other = item[2]
                    extra = None
                    if (len(batch) < self.max_batch and other.job_id not in self._cancelled
                            and other.batch != "small"                      # a small batch is always its own head
                            and other.round.round_id == job.round.round_id and self._branch(other) == bk
                            and (not served or self._served_warm(other))):
                        extra = extra_seconds(self.runner, job, other, self._learned_extra)
                        if extra is not None and seconds + extra > limit_s:
                            extra = None
                    if extra is not None and seconds + extra <= room:
                        batch.append(other)
                        seconds += extra
                    else:
                        keep.append(item)
                for item in keep:
                    heapq.heappush(self._queue, item)
            if len(batch) > 1 and not served:
                self.batches_run += 1
            self._running[self._pool.submit(self._run, batch, None if served else seconds)] = batch
        self._held = len(skipped)
        for job in skipped:
            heapq.heappush(self._queue, (-job.priority, next(self._seq), job))

    def submit(self, jobs: Iterable[EvalJob]) -> None:
        with self._lock:
            for job in jobs:
                self._route(self._as_canonical(job))
                self._submitted += 1
                self._last_submit = self._now()
            self._fill()
            self._sample()
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
                    if job.round.round_id in self._ended or job.job_id in self._dropped:   # nobody wants it
                        self._forget(job)
                        continue
                    if (not cs.ok and (cs.error or "").startswith("rc=3") and len(jobs) > 1
                            and job.job_id not in self._retried):
                        # GATK died in its own code inside a batch (the patch's variant path, for some heads on some
                        # regions; the same config passes on its own). The head runs again alone, still in flight;
                        # its riders fail here (riders only exist inside a batch).
                        self._retried.add(job.job_id)
                        if job is jobs[0]:
                            self._solo.add(job.job_id)
                            self.retries += 1
                            heapq.heappush(self._queue, (-job.priority, next(self._seq), job))
                            continue
                        cs = replace(cs, error=(cs.error or "") + " | batch crashed; its head runs alone")
                    if not cs.ok and cs.error == DIVERGED:
                        # its batch output was not exact: repair it if the runner can, else run it alone
                        # (still counts as in flight)
                        self._solo.add(job.job_id)
                        if hasattr(self.runner, "repair") and cs.vcf:
                            self._repair[job.job_id] = cs
                        heapq.heappush(self._queue, (-job.priority, next(self._seq), job))
                        continue
                    if cs.ok:
                        self._learn_cap(job, cs)
                    self._finish(job, cs)
            self._fill()
            out, self._ready = self._ready, []
            self.done += len(out)
            self._sample()
            return out

    def _forget(self, job: EvalJob) -> None:
        """A result of an ended round: nobody wants it. Caller holds the lock."""
        k = job_key(job)
        self._inflight.pop(k, None)
        for w in self._waiters.pop(k, []):
            self._as_key.pop(w.job_id, None)
        self._as_key.pop(job.job_id, None)
        self._repair.pop(job.job_id, None)

    def _finish(self, job: EvalJob, cs: CallSet) -> None:
        """Deliver a result to its job and every duplicate waiting on it. Caller holds the lock."""
        k = job_key(job)
        self._inflight.pop(k, None)
        if cs.ok:
            self._done[k] = cs                         # under the config that ran
        self._ready.append(self._restore(job.job_id, cs))
        self._ready += [self._restore(w.job_id, cs) for w in self._waiters.pop(k, [])]

    def _sample(self) -> None:
        """At most once a minute, record where the workers' time is going (caller holds the lock).

        busy / idle: worker slots running a GATK run or not. queued / families: jobs waiting and the GATK runs they
        would make. held: queued jobs the deadline keeps from starting. submitted / since_submit: what search sent in
        the last minute and how long ago it last sent anything. Idle workers with an empty queue and a long
        since_submit mean search is not proposing (or waiting for scores); idle workers with held jobs mean the
        deadline.
        """
        now = self._now()
        if now < self._next_sample:
            return
        self._next_sample = now + TIMELINE_S
        families = {(it[2].round.round_id, self._branch(it[2]) or it[2].job_id) for it in self._queue}
        tags = Counter(j.tag for batch in self._running.values() for j in batch)
        queued_tags = Counter(it[2].tag for it in self._queue)
        busy = len(self._running)
        sample = {"t": round(now, 1), "workers": self.workers, "busy": busy, "idle": max(0, self.workers - busy),
                  "running_configs": sum(len(b) for b in self._running.values()), "queued": len(self._queue),
                  "families": len(families), "held": self._held, "submitted": self._submitted,
                  "since_submit": None if self._last_submit is None else round(now - self._last_submit, 1),
                  "canonicalized": self.canonicalized,
                  "done": self.done, "gatk_runs": self.gatk_runs, "warm_runs": self.warm_runs,
                  "scale": round(self._scale(), 3),
                  "time_left": None if self._finish_by is None else round(self._finish_by - now, 1),
                  "by_tag": {t: [tags.get(t, 0), queued_tags.get(t, 0)] for t in sorted(set(tags) | set(queued_tags))}}
        self._submitted = 0
        self._timeline.append(sample)
        log.info("executor: busy %d/%d queued %d (%d families, %d held) submitted %d since_submit %s done %d "
                 "scale %.2f by_tag %s", busy, self.workers, sample["queued"], sample["families"], sample["held"],
                 sample["submitted"], sample["since_submit"], self.done, sample["scale"], sample["by_tag"])

    # ---- queries ---------------------------------------------------------------------------------------------
    def cancel(self, job_ids: Iterable[str]) -> None:
        with self._lock:
            self._cancelled.update(job_ids)

    def cancel_running(self, job_ids: Iterable[str]) -> None:
        """Jobs whose result is no longer wanted even if they are running (a duplicate run elsewhere delivered
        first): queued ones are dropped, a running batch made only of them is killed when the runner can
        (`kill_configs`), and nothing is delivered for them."""
        ids = set(job_ids)
        targets: Dict[str, List[str]] = {}
        with self._lock:
            self._cancelled.update(ids)
            self._dropped.update(ids)
            for jobs in self._running.values():
                if all(j.job_id in ids for j in jobs):      # a batch with wanted results runs on
                    targets.setdefault(jobs[0].round.round_id, []).extend(j.config.key for j in jobs)
        kill = getattr(self.runner, "kill_configs", None)
        if kill:
            for rid, keys in targets.items():
                kill(rid, keys)

    def stats(self) -> dict:
        with self._lock:
            return {"kind": "local", "workers": self.workers, "running": len(self._running), "queued": len(self._queue),
                    "batches_run": self.batches_run, "repairs_run": self.repairs_run, "retries": self.retries, "done": self.done,
                    "gatk_runs": self.gatk_runs, "gatk_runs_done": self.gatk_runs_done, "warm_runs": self.warm_runs,
                    "solo_s": round(self._deadline_s()), "scale": round(self._scale(), 3),
                    "extra_s": dict(self._learned_extra), "timeline": list(self._timeline),
                    "canonicalized": self.canonicalized, "cap_classes": sum(len(v) for v in self._caps.values())}

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
            self._facts.pop(round_id, None)
            for cls in [c for c in self._caps if c[0] == round_id]:
                del self._caps[cls]
        forget = getattr(self.runner, "forget_kbest", None)
        if forget:
            forget(round_id)
        kill = getattr(self.runner, "kill_round", None)
        if kill:
            kill(round_id)

    def cost_hint(self, job: EvalJob) -> float:
        """0 if done or in flight; if a queued job would carry it in its batch, the runner's variant cost
        (a fraction of a solo run, by the most expensive GATK stage it changes); else 1."""
        with self._lock:
            job = replace(job, config=self._canonical_config(job))
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
