"""FleetExecutor: the Executor on main when GATK and scoring run on sub-VPSs, controlled over SSH.

Same protocol as LocalExecutor (submit / poll / cancel / cost_hint(s) / idle_slots), so search and main do not
change. For every enabled host in the inventory (fleet.json, repo folder) main keeps one SSH session that runs the node
agent (evaluator/fleet/node.py --stdio); jobs go down and scored callsets come back as JSON lines on that session.

- joining: main pushes the self-test bundle, starts the agent, and accepts the node only if its self-test passed and
  its code digest equals main's
- rounds: publish_round() ships the round itself (the BAM with rsync in replays; on live rounds nodes fetch it from
  the presigned URL) and nodes start GATK at once; publish_answer_set() ships the truth files when the answer set
  exists (and again on every new version) - results made before it are held on the node and scored then
- scheduling: one queue with dedup and a result cache; a family (configs that can share a GATK run) is always sent as
  one group; a node gets groups while it has fewer open groups than workers. A small batch (a result the search
  waits for) goes to the fastest free node that has its round, or to a slow one when no fast node has a free
  worker; one running late is duplicated on a free faster node while the fleet has idle workers: the first result
  wins and the other run is aborted
- failure: when a session ends, that node's unfinished jobs go straight back to the queue; main reconnects later
- the inventory is re-read at every round start (publish_round): new hosts join, disabled ones are closed
"""
from __future__ import annotations

import heapq
import itertools
import json
import logging
import queue
import shlex
import subprocess
import threading
import time
from collections import defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from common.settings import SETTINGS
from common.types import AnswerSet, CallSet, EvalJob, Round
from evaluator.fleet import codec, selftest
from evaluator.fleet.inventory import INVENTORY, Host, code_digest, ensure_agent, load
from common.tuning import tune
from evaluator.local import SOLO_S, batch_limit, failed, job_key, max_batch_s

log = logging.getLogger(__name__)
_UNSET = object()
LOG_DIR = Path.home() / ".minos" / "fleet" / "logs"
SPEEDS = Path.home() / ".minos" / "fleet" / "speeds.json"   # host -> measured solo time, so routing is right from the start
SOLO_EMA = 0.3                # weight of each new solo-time sample in a node's running estimate
KEEP_ROUNDS = 3               # rounds kept on main and nodes (older ones are pruned)
# Tunable (common/tuning.py): evaluator.overdue (a small batch running this many times its expected time is a
# straggler and may be copied), evaluator.slow_factor (a node slower than this times the fastest is "slow": small
# batches prefer the fast ones), evaluator.small_batch / evaluator.big_batch (evaluator.local.batch_limit).


@dataclass
class _Published:
    round: Round
    answer_set: Optional[AnswerSet] = None                   # None until main has built it
    bam_url: Optional[str] = None
    bai_url: Optional[str] = None


@dataclass
class _Group:
    """One GATK run on a node: a head and its riders."""
    jobs: Dict[str, EvalJob]                                 # job_id -> job, still wanted from this run
    ids: Tuple[str, ...]                                     # every job sent in it
    units: float                                             # predicted cost in solo runs
    t0: float
    batch: str                                               # small | big (its head's)
    due: int                                                 # results this node still owes for it
    duplicate: bool = False                                  # a speculative copy of a run elsewhere
    twin: Optional[Tuple["_Session", int]] = field(default=None, repr=False)   # the other copy (session, group id)


@dataclass
class _Session:
    host: Host
    proc: Optional[subprocess.Popen] = None
    workers: int = 0
    ready: bool = False
    rounds_ready: set = field(default_factory=set)
    rounds_sent: set = field(default_factory=set)            # round_ids shipped to the node
    asets_sent: set = field(default_factory=set)             # (round_id, answer-set version) shipped to the node
    rounds_failed: set = field(default_factory=set)          # round_ids the node could not prepare or score
    groups: Dict[int, _Group] = field(default_factory=dict)  # open GATK runs by group id
    jobs: Dict[str, EvalJob] = field(default_factory=dict)   # outstanding job_id -> job
    job_group: Dict[str, int] = field(default_factory=dict)  # job_id -> group id
    group_seq: int = 0
    done: int = 0                                            # results delivered (dashboard)
    gatk_runs: int = 0                                       # groups handed to it: one GATK run each (dashboard)
    gatk_runs_done: int = 0                                  # ... whose every result came back from it
    outq: "queue.Queue[Optional[str]]" = field(default_factory=queue.Queue)
    closed: bool = False
    finished: bool = False                                   # its session thread has ended
    rejected: str = ""
    logf: Optional[object] = None                            # the node's stderr, kept per host
    solo_s: float = SOLO_S                                   # measured: seconds one solo GATK run takes THERE (EMA)
    solo_samples: int = 0

    @property
    def name(self) -> str:
        return self.host.name if self.host is not None else "?"

    def open_groups(self) -> int:
        return len(self.groups)

    def note_group_done(self, units: float, seconds: float) -> None:
        """A whole group came back after `seconds`: one solo-run-equivalent took seconds / units here. Learned
        from small groups only (at most a small batch's work): the cost table is least wrong there, and big
        batches' errors made healthy nodes look 2-3x slower in the v4 fleet run."""
        if units > batch_limit("small") + 1e-6:
            return
        sample = min(3600.0, max(0.05, seconds / max(units, 0.05)))
        self.solo_s = sample if self.solo_samples == 0 else (1 - SOLO_EMA) * self.solo_s + SOLO_EMA * sample
        self.solo_samples += 1

    def send(self, msg: dict) -> bool:
        """Queue a message for the node. Writing happens on the session's writer thread, so a slow link never
        blocks main's scheduler (which calls this while holding the executor lock)."""
        if self.closed or self.proc is None:
            return False
        self.outq.put(json.dumps(msg) + "\n")
        return True

    def writer(self, proc: subprocess.Popen) -> None:
        while True:
            line = self.outq.get()
            if line is None:
                break
            try:
                proc.stdin.write(line)
                proc.stdin.flush()
            except (BrokenPipeError, ValueError, OSError):
                break                                          # the reader sees the session end and re-queues
        try:
            proc.stdin.close()
        except OSError:
            pass

    def hang_up(self) -> None:
        """Close the node's stdin once the queued messages are written (the agent then exits)."""
        self.outq.put(None)


class FleetExecutor:
    def __init__(self, planner, inventory: Path = INVENTORY, max_batch: int = 32, reconnect_s: float = 30.0):
        """planner: an object with branch_key/variant_cost (a GatkRunner with patch_jar set) used only to group and
        price families - nothing runs on main."""
        self.runner = planner
        self.inventory, self.max_batch, self.reconnect_s = inventory, max_batch, reconnect_s
        self.code = code_digest()
        self._lock = threading.RLock()
        self._cv = threading.Condition(self._lock)
        self._seq = itertools.count()
        self._queue: list = []                                   # heap of (-priority, seq, job); taken items stay until popped
        self._queued: set = set()                                # job ids actually queued (taken riders leave the heap lazily)
        self._by_branch: Dict[Tuple[str, Optional[str]], Dict[str, EvalJob]] = defaultdict(dict)   # (round, branch) -> queued jobs
        self._branch_of: Dict[str, Optional[str]] = {}           # job id -> branch key, computed once
        self._inflight: Dict[Tuple, EvalJob] = {}
        self._waiters = defaultdict(list)
        self._done: Dict[Tuple, CallSet] = {}
        self._ready: List[CallSet] = []
        self._cancelled = set()
        self._rounds: Dict[str, _Published] = {}
        self._sessions: Dict[str, _Session] = {}
        self._ended = set()                                      # round ids ended by end_round()
        self._stopping = False
        self._now = time.time
        self._finish_by: Optional[float] = None                  # do not send a batch that would still be running then
        self.speculated = 0                                      # latency groups duplicated on a faster node
        self._speeds: Dict[str, float] = self._load_speeds()
        self.speculation_won = 0                                 # ... whose duplicate delivered a result first
        LOG_DIR.mkdir(parents=True, exist_ok=True)

    # ---- hosts --------------------------------------------------------------------------------------------------------
    @staticmethod
    def _load_speeds() -> Dict[str, float]:
        try:
            return {k: float(v) for k, v in json.loads(SPEEDS.read_text()).items()}
        except (OSError, ValueError):
            return {}

    def _save_speeds(self) -> None:
        """Every node's measured solo time, by host: the next run routes right from its first minute."""
        for s in self._sessions.values():
            if s.solo_samples:
                self._speeds[s.name] = round(s.solo_s, 1)
        try:
            SPEEDS.parent.mkdir(parents=True, exist_ok=True)
            tmp = SPEEDS.with_suffix(".tmp")
            tmp.write_text(json.dumps(self._speeds, indent=1))
            tmp.replace(SPEEDS)
        except OSError:
            pass

    def start(self) -> "FleetExecutor":
        self.refresh_hosts()
        return self

    def refresh_hosts(self) -> None:
        """Re-read the inventory: start sessions for new enabled hosts, close the ones no longer enabled."""
        hosts = {h.name: h for h in load(self.inventory)}
        ensure_agent(list(hosts.values()))                    # unlock passphrase-protected keys once
        with self._lock:
            for name, s in list(self._sessions.items()):
                if name not in hosts:
                    self._drop(s, "disabled in the inventory")
                    s.closed = True
                    del self._sessions[name]
            for name, h in hosts.items():
                old = self._sessions.get(name)
                if old is not None and old.rejected and old.finished:
                    del self._sessions[name]                    # rejected earlier (old code?): try again now
                if name not in self._sessions:
                    s = _Session(h)
                    self._sessions[name] = s
                    threading.Thread(target=self._run_session, args=(s,), daemon=True, name=f"fleet-{name}").start()

    def _run_session(self, s: _Session) -> None:
        try:
            while not s.closed and not self._stopping:
                try:
                    self._connect(s)
                    self._read(s)
                except Exception as e:  # noqa: BLE001 - one host failing must not stop the fleet
                    log.warning("%s: session error: %s", s.name, e)
                with self._lock:
                    self._drop(s, "session ended")
                if s.rejected or s.closed or self._stopping:
                    return
                time.sleep(self.reconnect_s)
        finally:
            s.finished = True

    def _connect(self, s: _Session) -> None:
        h = s.host
        bundle = selftest.make_bundle(SETTINGS.cache_dir, SETTINGS.rounds_dir)
        push(h, selftest.bundle_files(bundle), f"{h.repo}/cache/selftest/")   # slice BAM + truth: the score check
        args = h.extra.get("node_args", "") + (f" --workers {int(h.workers)}" if h.workers else "")
        shared = {"MINOS_NODE_SHARED_DATA": "1"} if h.local else {}          # main's rounds are not the node's to delete
        env = " ".join(f"{k}={shlex.quote(str(v))}" for k, v in {"DATA_DIR": h.data_dir, **shared, **h.extra.get("env", {})}.items())
        repo = shlex.quote(h.repo)
        cmd = h.ssh(f"cd {repo} && {env} {repo}/.venv/bin/python -m evaluator.fleet.node --stdio {args}")
        s.logf = open(LOG_DIR / f"{h.name.replace(':', '_')}.log", "a")
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=s.logf, text=True, bufsize=1)
        s.outq = queue.Queue()
        threading.Thread(target=s.writer, args=(proc,), daemon=True, name=f"fleet-w-{h.name}").start()
        with self._lock:
            s.proc = proc
            if s.closed or self._stopping:                       # disabled/closed while we were connecting
                s.hang_up()
        log.info("%s: session started", h.name)

    def _read(self, s: _Session) -> None:
        for line in s.proc.stdout:
            try:
                msg = json.loads(line)
            except json.JSONDecodeError:
                continue
            with self._cv:
                self._handle(s, msg)
                self._cv.notify_all()
        s.proc.wait()
        if s.logf:
            s.logf.close()

    def _handle(self, s: _Session, msg: dict) -> None:
        t = msg["type"]
        if s.closed:
            return
        if t == "ready":
            if not msg.get("ok"):
                s.rejected = selftest.reject_reason(msg)
            elif msg.get("code") != self.code:
                s.rejected = f"code differs from main ({msg.get('code')} != {self.code}): run control sync"
            if s.rejected:
                log.error("%s: rejected - %s", s.name, s.rejected)
                s.hang_up()
                return
            s.workers, s.ready = int(msg["workers"]), True
            known = self._speeds.get(s.name)
            if known and s.solo_samples == 0:                     # last time's measurement, until this run has its own
                s.solo_s = float(known)
            log.info("%s: ready, %d workers%s", s.name, s.workers, f", {s.solo_s:.0f} s per run last time" if known else "")
            threading.Thread(target=self._send_rounds, args=(s,), daemon=True).start()
        elif t == "round_ready":
            s.rounds_ready.add(msg["round_id"])
            self._pump()
        elif t == "round_failed":
            log.error("%s: round %s unusable there: %s", s.name, msg["round_id"], msg.get("error"))
            self._unserve(s, msg["round_id"])
        elif t == "result":
            self._result(s, codec.callset_from(msg["callset"]))

    def _result(self, s: _Session, cs: CallSet) -> None:
        """A scored callset from node `s`. The first result for a job wins (a speculative duplicate elsewhere is then
        aborted); a group whose every result came from this node times the node. Lock held."""
        if cs.round_id in self._ended:
            return                                                # the round is over; its jobs were already freed
        job = s.jobs.pop(cs.job_id, None)
        gid = s.job_group.pop(cs.job_id, None)
        g = s.groups.get(gid) if gid is not None else None
        if g is not None:
            g.jobs.pop(cs.job_id, None)
            g.due -= 1
            if g.due <= 0:                                        # every result of the run came back from here
                if not g.duplicate:                               # a copy's short time would teach the node it is faster than it is
                    s.note_group_done(g.units, self._now() - g.t0)
                s.gatk_runs_done += 1
                del s.groups[gid]
                if s.solo_samples % 10 == 1:
                    self._save_speeds()
            elif not g.jobs:                                      # the rest was delivered by its twin
                del s.groups[gid]
        if job is not None:
            s.done += 1
            if g is not None and g.duplicate and job_key(job) not in self._done:
                self.speculation_won += 1
            self._finish(cs, job)
            if g is not None and g.twin is not None:
                self._twin_served(g, cs.job_id)
        self._pump()

    def _twin_served(self, g: _Group, job_id: str) -> None:
        """`job_id` came back from `g`'s node first: its copy in the twin run is not needed; once nothing of that run
        is needed any more, its node aborts it (the GATK container is killed, the worker freed)."""
        other, ogid = g.twin
        og = other.groups.get(ogid)
        if og is None:
            g.twin = None
            return
        og.jobs.pop(job_id, None)
        other.jobs.pop(job_id, None)
        other.job_group.pop(job_id, None)
        if not og.jobs:
            del other.groups[ogid]
            other.send({"type": "abort_jobs", "job_ids": list(og.ids)})
            g.twin = None

    def _drop(self, s: _Session, why: str) -> None:
        """Session over: its unfinished jobs go back to the queue at once - except those its twin run elsewhere is
        still computing. Caller holds the lock."""
        if s.jobs:
            log.warning("%s: %s - %d jobs re-queued", s.name, why, len(s.jobs))
        for g in s.groups.values():
            running_elsewhere: set = set()
            if g.twin is not None:
                other, ogid = g.twin
                og = other.groups.get(ogid)
                if og is not None:
                    og.twin = None                                # the other copy goes on alone
                    running_elsewhere = set(og.jobs)
            for jid, job in g.jobs.items():
                if jid not in running_elsewhere and job_key(job) not in self._done:
                    self._push(job)
        s.jobs, s.groups, s.job_group, s.ready = {}, {}, {}, False
        s.rounds_ready, s.rounds_sent, s.asets_sent, s.rounds_failed = set(), set(), set(), set()
        if s.proc and s.proc.poll() is None:
            s.hang_up()
        self._cv.notify_all()
        self._pump()

    def _unserve(self, s: _Session, rid: str) -> None:
        """Node `s` cannot serve round `rid` (it could not prepare it, or its truth never arrived there): nothing
        more of it goes there, its unfinished work of the round goes back to the queue, and if no ready node can
        serve the round its queued jobs fail. Lock held."""
        s.rounds_failed.add(rid)
        s.rounds_ready.discard(rid)
        mine = []
        for gid, g in list(s.groups.items()):
            if any(j.round.round_id == rid for j in g.jobs.values()):
                mine += list(g.jobs)
                for jid, job in g.jobs.items():
                    if job_key(job) not in self._done:
                        self._push(job)
                if g.twin is not None and g.twin[0].groups.get(g.twin[1]) is not None:
                    g.twin[0].groups[g.twin[1]].twin = None
                del s.groups[gid]
        for jid in mine:
            s.jobs.pop(jid, None)
            s.job_group.pop(jid, None)
        if mine:
            s.send({"type": "abort", "round_id": rid})
        self._fail_if_unservable(rid)
        self._pump()

    def _fail_if_unservable(self, rid: str) -> None:
        """No ready node can prepare this round: its queued jobs fail now instead of waiting forever. Lock held."""
        ready = [s for s in self._sessions.values() if s.ready]
        if not ready or not all(rid in s.rounds_failed for s in ready):
            return
        keep, n = [], 0
        for item in self._queue:
            job = item[2]
            if job.round.round_id == rid:
                self._inflight.pop(job_key(job), None)
                self._ready.append(failed(job, "no fleet node could prepare this round"))
                self._ready += [failed(w, "no fleet node could prepare this round") for w in self._waiters.pop(job_key(job), [])]
                n += 1
            else:
                keep.append(item)
        self._queue = keep
        heapq.heapify(self._queue)
        self._reindex()
        if n:
            log.error("round %s: no node can prepare it - %d queued jobs failed", rid, n)
            self._cv.notify_all()

    # ---- rounds -------------------------------------------------------------------------------------------------------
    def publish_round(self, round: Round, bam_url: Optional[str] = None, bai_url: Optional[str] = None) -> None:
        """The round itself - its BAM (pushed from main in replays, fetched from the presigned URL by each node on
        live rounds) and region - before its answer set exists: nodes start GATK at once and hold the results."""
        self.refresh_hosts()
        with self._lock:
            self._rounds.pop(round.round_id, None)
            self._rounds[round.round_id] = _Published(round, None, bam_url, bai_url)
            while len(self._rounds) > KEEP_ROUNDS:              # oldest first (dicts keep insertion order)
                self._rounds.pop(next(iter(self._rounds)))
            sessions = [s for s in self._sessions.values() if s.ready]
        for s in sessions:
            threading.Thread(target=self._send_rounds, args=(s,), daemon=True).start()

    def publish_answer_set(self, round_id: str, answer_set: AnswerSet) -> None:
        """The round's truth files (a first version, or a newer one): every node scores with it from now on, starting
        with the results it held."""
        with self._lock:
            p = self._rounds.get(round_id)
            if p is None:
                log.warning("answer set for an unpublished round %s ignored", round_id)
                return
            p.answer_set = answer_set
            sessions = [s for s in self._sessions.values() if s.ready]
        for s in sessions:
            threading.Thread(target=self._send_rounds, args=(s,), daemon=True).start()

    def _send_rounds(self, s: _Session) -> None:
        """Ship to node `s` whatever it has not got yet: each kept round (once), then its answer set (once per
        version). Files go by rsync before the message; a failed shipment makes the round unservable there."""
        with self._lock:
            rounds = list(self._rounds.items())
        keep = [rid for rid, _ in rounds]
        for rid, p in rounds:
            with self._lock:
                if not s.ready:
                    return
                send_round = rid not in s.rounds_sent and rid not in s.rounds_failed
                if send_round:
                    s.rounds_sent.add(rid)
            if send_round:
                try:
                    s.send({**round_message(s.host, p.round, p.bam_url, p.bai_url), "keep": keep})
                except Exception as e:  # noqa: BLE001 - e.g. rsync failed: this node cannot serve the round
                    log.error("%s: could not ship round %s: %s", s.name, rid, e)
                    with self._lock:
                        s.rounds_sent.discard(rid)
                        self._unserve(s, rid)
                    continue
            with self._lock:
                a = p.answer_set
                key = (rid, a.version) if a is not None else None
                send_aset = key is not None and key not in s.asets_sent and rid not in s.rounds_failed
                if send_aset:
                    s.asets_sent.add(key)
            if send_aset:
                try:
                    s.send(answer_set_message(s.host, p.round, a))
                except Exception as e:  # noqa: BLE001
                    log.error("%s: could not ship the answer set of %s: %s", s.name, rid, e)
                    with self._lock:
                        s.asets_sent.discard(key)
                        self._unserve(s, rid)

    # ---- Executor protocol ----------------------------------------------------------------------------------------
    def set_finish_by(self, ts: Optional[float]) -> None:
        """Do not send a batch that would still be running at `ts` (the final submission). None clears it."""
        with self._lock:
            self._finish_by = ts

    def _room(self, s: Optional["_Session"] = None) -> float:
        """Seconds a batch started now on node `s` may take: the family cap (max_batch_s() for a node at the planning
        solo time, scaled by the node's own measured solo time) and the time left before the final submission."""
        room = max_batch_s() * (s.solo_s / SOLO_S if s is not None else 1.0)
        if self._finish_by is not None:
            room = min(room, self._finish_by - self._now())
        return room

    def _branch(self, job: EvalJob) -> Optional[str]:
        bk = self._branch_of.get(job.job_id, _UNSET)
        if bk is _UNSET:
            f = getattr(self.runner, "branch_key", None)
            bk = f(job.config, job.fidelity) if f else None
            self._branch_of[job.job_id] = bk
        return bk

    # ---- the queue: a priority heap plus an index by (round, branch) so a head finds its riders without a scan ------
    def _push(self, job: EvalJob) -> None:
        heapq.heappush(self._queue, (-job.priority, next(self._seq), job))
        self._queued.add(job.job_id)
        self._by_branch[(job.round.round_id, self._branch(job))][job.job_id] = job

    def _taken(self, job: EvalJob) -> None:
        """`job` leaves the queue (dispatched as a head or a rider, or dropped); its heap item is skipped later."""
        self._queued.discard(job.job_id)
        self._by_branch[(job.round.round_id, self._branch(job))].pop(job.job_id, None)

    def _reindex(self) -> None:
        """After the heap was rebuilt from a `keep` list: the index follows."""
        self._queued = {it[2].job_id for it in self._queue}
        self._by_branch = defaultdict(dict)
        for it in self._queue:
            self._by_branch[(it[2].round.round_id, self._branch(it[2]))][it[2].job_id] = it[2]

    def submit(self, jobs: Iterable[EvalJob]) -> None:
        with self._cv:
            for job in jobs:
                k = job_key(job)
                if k in self._done:
                    self._ready.append(replace(self._done[k], job_id=job.job_id))
                elif k in self._inflight:
                    self._waiters[k].append(job)
                else:
                    self._inflight[k] = job
                    self._push(job)
            self._pump()
            self._cv.notify_all()

    def poll(self, timeout: float) -> List[CallSet]:
        deadline = time.time() + timeout
        with self._cv:
            self._speculate()                                     # main polls every second: the straggler check
            while not self._ready and time.time() < deadline:
                self._cv.wait(max(0.0, deadline - time.time()))
            out, self._ready = self._ready, []
            return out

    def cancel(self, job_ids: Iterable[str]) -> None:
        with self._lock:
            self._cancelled.update(job_ids)

    def end_round(self, round_id: str) -> None:
        """The round is over: drop its queued jobs, tell every node holding its work to abort (the node kills those
        GATK containers), and free those workers at once for the next round. Late results for it are ignored."""
        with self._cv:
            self._ended.add(round_id)
            keep = []
            for item in self._queue:
                job = item[2]
                if job.round.round_id == round_id:
                    self._inflight.pop(job_key(job), None)
                    self._waiters.pop(job_key(job), None)
                else:
                    keep.append(item)
            self._queue = keep
            heapq.heapify(self._queue)
            self._reindex()
            for s in self._sessions.values():
                mine = []
                for gid, g in list(s.groups.items()):
                    if any(j.round.round_id == round_id for j in g.jobs.values()):
                        mine += list(g.jobs)
                        del s.groups[gid]
                if not mine:
                    continue
                for jid in mine:
                    self._inflight.pop(job_key(s.jobs.pop(jid)), None)
                    s.job_group.pop(jid, None)
                s.send({"type": "abort", "round_id": round_id})
            self._pump()

    def cost_hint(self, job: EvalJob) -> float:
        return self.cost_hints([job])[0]

    def cost_hints(self, jobs: List[EvalJob]) -> List[float]:
        with self._lock:
            out = []
            for i, job in enumerate(jobs):
                k = job_key(job)
                if k in self._done or k in self._inflight:
                    out.append(0.0)
                    continue
                c, bk = 1.0, self._branch(job)
                if bk is not None and hasattr(self.runner, "variant_cost"):
                    same = list(self._by_branch.get((job.round.round_id, bk), {}).values())
                    same += [p for p in jobs[:i] if self._branch(p) == bk]
                    for p in same:
                        v = self.runner.variant_cost(p.config, job.config, job.fidelity)
                        if v is not None:
                            c = min(c, v)
                out.append(c)
            return out

    def capacity(self) -> int:
        """GATK workers on every ready node."""
        with self._lock:
            return sum(s.workers for s in self._sessions.values() if s.ready)

    def idle_slots(self) -> int:
        with self._lock:
            free = sum(max(0, s.workers - s.open_groups()) for s in self._sessions.values() if s.ready)
            return max(0, free - len(self._queued))

    def stats(self) -> dict:
        """For the public dashboard: nodes are numbered, never named (no IPs)."""
        with self._lock:
            ready = [s for s in self._sessions.values() if s.ready]
            fastest = min((s.solo_s for s in ready), default=SOLO_S)
            slow, sick = tune("evaluator.slow_factor") * fastest, self._sick_limit(ready)
            nodes = [dict(node=f"node {i}", ready=s.ready, workers=s.workers, busy=s.open_groups(), jobs=len(s.jobs),
                          small=sum(1 for g in s.groups.values() if g.batch == "small"),
                          done=s.done, gatk_runs=s.gatk_runs, gatk_runs_done=s.gatk_runs_done,
                          rejected=bool(s.rejected), solo_s=round(s.solo_s), slow=s.ready and s.solo_s > slow,
                          sick=s.ready and s.solo_s > sick)
                     for i, s in enumerate(self._sessions.values(), 1)]
            cores: dict = {}
            for s in self._sessions.values():
                if not s.ready:
                    continue
                for g in s.groups.values():                    # one open group is one core
                    tag = next((j.tag for j in g.jobs.values() if j.tag), "")
                    cores[tag] = cores.get(tag, 0) + 1
            return {"kind": "fleet", "queued": len(self._queued), "nodes": nodes, "cores": cores,
                    "done": sum(n["done"] for n in nodes),
                    "gatk_runs": sum(n["gatk_runs"] for n in nodes),
                    "gatk_runs_done": sum(n["gatk_runs_done"] for n in nodes),
                    "workers": sum(n["workers"] for n in nodes if n["ready"]),
                    "running": sum(n["busy"] for n in nodes if n["ready"]),
                    "speculated": self.speculated, "speculation_won": self.speculation_won}

    def nodes(self) -> List[dict]:
        with self._lock:
            return [dict(host=n, ready=s.ready, workers=s.workers, open_groups=s.open_groups(), jobs=len(s.jobs),
                         rounds=len(s.rounds_ready), rejected=s.rejected, solo_s=round(s.solo_s, 1))
                    for n, s in self._sessions.items()]

    def close(self) -> None:
        """Stop every node (closing its stdin ends the agent, which kills its GATK containers)."""
        self._stopping = True
        with self._lock:
            for s in self._sessions.values():
                s.closed = True
                if s.proc and s.proc.poll() is None:
                    s.hang_up()

    # ---- scheduling -------------------------------------------------------------------------------------------------
    def _peek(self) -> Optional[EvalJob]:
        """The top live queued job (cancelled or done ones are dropped on the way), or None."""
        while self._queue:
            job = self._queue[0][2]
            if job.job_id not in self._queued or self._skip(job):
                heapq.heappop(self._queue)
                self._taken(job)
                continue
            return job
        return None

    def _pump(self) -> None:
        """Hand queued groups to nodes. The top job picks the node: a small batch goes to the fastest node with a
        free worker that has its round (its result is on the search's critical path), or to a slow node when no fast
        one has a free worker - a slow node never sits idle while small batches wait; a big batch goes to the node
        with the largest free share, so work spreads over the fleet instead of filling the first node.
        Caller holds the lock."""
        stuck = set()                                          # nodes that could serve nothing this pass
        ready = [s for s in self._sessions.values() if s.ready]
        fastest = min((s.solo_s for s in ready), default=SOLO_S)
        sick = self._sick_limit(ready)
        while self._queue:
            job = self._peek()
            if job is None:
                break
            free = [s for s in self._sessions.values()
                    if s.ready and id(s) not in stuck and s.open_groups() < s.workers]
            able = [s for s in free if job.round.round_id in s.rounds_ready]
            big_only = False
            if able and job.batch == "small":
                well = [s for s in able if s.solo_s <= sick]
                fast = [s for s in well if s.solo_s <= tune("evaluator.slow_factor") * fastest]
                if well:
                    s = min(fast or well, key=lambda x: x.solo_s)   # no fast node free: the fastest free slow one
                else:                                               # only sick nodes free: they take big work
                    s, big_only = max(able, key=lambda x: (x.workers - x.open_groups()) / max(1, x.workers)), True
            elif able:
                s = max(able, key=lambda x: (x.workers - x.open_groups()) / max(1, x.workers))
            elif free:
                s = free[0]                                    # the top job's round is not ready there: serve what is
            else:
                break
            if not self._dispatch(s, big_only):
                stuck.add(id(s))

    def _sick_limit(self, ready: List["_Session"]) -> float:
        """Solo time above which a node is sick (evaluator.sick_factor x the ready nodes' median)."""
        times = sorted(s.solo_s for s in ready)
        if not times:
            return float("inf")
        return tune("evaluator.sick_factor") * times[len(times) // 2]

    def _dispatch(self, s: _Session, big_only: bool = False) -> bool:
        """One group to node `s` (big batches only when `big_only`); False if it had nothing it could take (or the
        session is gone)."""
        groups = self._take_groups(s, big_only)
        if not groups:
            return False
        g, units = groups[0]
        self._send_group(s, g, units)
        return not s.closed and s.proc is not None

    def _send_group(self, s: _Session, jobs: List[EvalJob], units: float, duplicate: bool = False,
                    twin: Optional[Tuple[_Session, int]] = None) -> int:
        """Send `jobs` as one GATK run to node `s`; returns its group id. Lock held."""
        gid = s.group_seq
        s.group_seq += 1
        s.groups[gid] = _Group({j.job_id: j for j in jobs}, tuple(j.job_id for j in jobs), units, self._now(),
                               jobs[0].batch, len(jobs), duplicate, twin)
        for j in jobs:
            s.jobs[j.job_id] = j
            s.job_group[j.job_id] = gid
        s.gatk_runs += 1
        s.send({"type": "jobs", "groups": [[codec.job_to(j) for j in jobs]],
                "limits": {"small": batch_limit("small"), "big": batch_limit("big")}})   # the node packs as main does
        return gid

    def _speculate(self) -> None:
        """A small batch that has run past the time this node itself takes is duplicated on a free faster node,
        while the fleet has a worker to spare (nothing queued could use it: _pump ran first). The first result
        wins; the other run is aborted. A copy is not started early just because a faster node is free.
        Lock held."""
        now = self._now()
        for s in list(self._sessions.values()):
            if not s.ready:
                continue
            for gid, g in list(s.groups.items()):
                if g.batch != "small" or g.twin is not None or g.duplicate or not g.jobs:
                    continue
                rid = next(iter(g.jobs.values())).round.round_id
                others = [o for o in self._sessions.values()
                          if o is not s and o.ready and rid in o.rounds_ready and o.open_groups() < o.workers]
                if not others:
                    continue
                f = min(others, key=lambda o: o.solo_s)
                expected, elapsed, fresh = s.solo_s * g.units, now - g.t0, f.solo_s * g.units
                if fresh > self._room(f) or elapsed <= tune("evaluator.overdue") * expected:
                    continue
                jobs = list(g.jobs.values())
                gid2 = self._send_group(f, jobs, g.units, duplicate=True, twin=(s, gid))
                g.twin = (f, gid2)
                self.speculated += 1
                log.info("%s: small batch of %d jobs at %.0f s of %.0f expected - duplicated on %s (%.0f s there)",
                         s.name, len(jobs), elapsed, expected, f.name, fresh)

    def _variant_frac(self, main: EvalJob, variant: EvalJob) -> Optional[float]:
        """Extra cost of `variant` inside `main`'s run, as a fraction of a solo run; None = cannot share it."""
        cost = getattr(self.runner, "variant_cost", None)
        if cost is None:
            return 0.0                                         # a planner without costs: the old batch-anything rule
        return cost(main.config, variant.config, main.fidelity)

    def _take_groups(self, s: _Session, big_only: bool = False) -> List[Tuple[List[EvalJob], float]]:
        """The top queued job whose round the node has and that can finish in time there, plus queued jobs that can
        share its GATK run, up to batch_limit() solo runs of work in total (and the time left); a small
        job never rides in another job's batch. Times are this node's own (its measured solo time scales every cost).
        Returns [(group, predicted cost in solo runs)]."""
        skipped, too_long, group = [], [], None
        room = self._room(s)
        solo = s.solo_s
        units = 1.0
        while self._queue and group is None:
            item = heapq.heappop(self._queue)
            job = item[2]
            if job.job_id not in self._queued:                  # taken as a rider earlier: a stale heap item
                continue
            if self._skip(job):
                self._taken(job)
                continue
            if job.round.round_id not in s.rounds_ready or (big_only and job.batch == "small"):
                skipped.append(item)
                continue
            if solo > room:
                too_long.append(job)
                self._taken(job)
                continue
            self._taken(job)
            group, seconds, bk = [job], solo, self._branch(job)
            limit = batch_limit(job.batch) * solo
            if bk is not None:
                for other in list(self._by_branch.get((job.round.round_id, bk), {}).values()):   # the riders: same traversal
                    if len(group) >= self.max_batch:
                        break
                    if other.job_id not in self._queued or other.batch == "small":
                        continue
                    if self._skip(other):
                        self._taken(other)
                        continue
                    frac = self._variant_frac(job, other)
                    if frac is None:
                        continue
                    extra = frac * solo
                    if seconds + extra > limit:
                        continue
                    if seconds + extra <= room:
                        self._taken(other)
                        group.append(other)
                        seconds += extra
                        units += extra / solo
        for it in skipped:
            heapq.heappush(self._queue, it)
        for job in too_long:
            self._push(job)
        return [(group, units)] if group else []

    def _skip(self, job: EvalJob) -> bool:
        k = job_key(job)
        if job.job_id in self._cancelled:
            waiters = self._waiters.pop(k, [])
            if waiters:
                self._inflight[k] = waiters[0]
                self._waiters[k] = waiters[1:]
                self._push(waiters[0])
            else:
                self._inflight.pop(k, None)
            return True
        return k in self._done

    def _finish(self, cs: CallSet, job: EvalJob) -> None:
        k = job_key(job)
        if k in self._done:
            return
        cs = replace(cs, job_id=job.job_id)
        self._inflight.pop(k, None)
        if cs.ok:
            self._done[k] = cs
        self._ready.append(cs)
        self._ready += [replace(cs, job_id=w.job_id) for w in self._waiters.pop(k, [])]


def push(h: Host, files: List[str], dest_dir: str) -> None:
    """Copy local files into dest_dir on the host (rsync over ssh; only what changed is sent)."""
    if h.local and all(str(Path(f).parent) == dest_dir.rstrip("/") for f in files):
        return                                                  # already there (this machine, same path)
    mk = subprocess.run(h.ssh(f"mkdir -p {shlex.quote(dest_dir)}"), capture_output=True, text=True, timeout=120)
    if mk.returncode != 0:
        raise RuntimeError(f"mkdir on {h.name}: {mk.stderr[-200:]}")
    p = subprocess.run(h.rsync(files, dest_dir), capture_output=True, text=True, timeout=3600)
    if p.returncode != 0:
        raise RuntimeError(f"rsync to {h.name}: {p.stderr[-300:]}")


def _node_paths(h: Host, r: Round) -> Tuple[str, str]:
    """(the round's directory, the reference FASTA) as the node sees them."""
    return f"{h.data_dir}/rounds/{r.tag}", f"{h.data_dir}/reference/{r.chrom}/{r.chrom}.fa"


def round_message(h: Host, r: Round, bam_url: Optional[str], bai_url: Optional[str]) -> dict:
    """Push this round's BAM to the host when main has it (a replay) and describe the round with the host's paths."""
    rdir, ref = _node_paths(h, r)
    if not (bam_url and bai_url):                              # replay: main's own copy of the BAM
        push(h, [r.bam, r.bam + ".bai"], rdir + "/")
    rnd = replace(r, bam=f"{rdir}/{Path(r.bam).name}", reference=ref, released_dir=None)
    return {"type": "round", "round": codec.round_to(rnd), "bam_url": bam_url, "bai_url": bai_url}


def answer_set_message(h: Host, r: Round, a: AnswerSet) -> dict:
    """Push this version's truth files to the host's own cache and describe the answer set with the host's paths."""
    _, ref = _node_paths(h, r)
    adir = f"{h.repo}/cache/fleet/{r.tag}/v{a.version}"
    truth = [a.notes["truth_vcf"], a.notes["mutations_vcf"]]
    push(h, truth + [f + ".tbi" for f in truth if Path(f + ".tbi").exists()], adir + "/")
    notes = {"truth_vcf": f"{adir}/{Path(a.notes['truth_vcf']).name}",
             "mutations_vcf": f"{adir}/{Path(a.notes['mutations_vcf']).name}", "region": r.region, "reference": ref}
    return {"type": "answer_set", "round_id": r.round_id,
            "answer_set": {"version": a.version, "n_targets": a.n_targets, "notes": notes}}
