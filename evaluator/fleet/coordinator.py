"""FleetExecutor: the Executor on main when GATK and scoring run on sub-VPSs, controlled over SSH.

Same protocol as LocalExecutor (submit / poll / cancel / cost_hint(s) / idle_slots), so search and main do not
change. For every enabled host in the inventory (fleet.json, repo folder) main keeps one SSH session that runs the node
agent (evaluator/fleet/node.py --stdio); jobs go down and scored callsets come back as JSON lines on that session.

- joining: main pushes the self-test bundle, starts the agent, and accepts the node only if its self-test passed and
  its code digest equals main's
- rounds: publish_round() pushes each round's truth files (and, in replays, the BAM) to every node with rsync; on
  live rounds nodes fetch the BAM from the presigned URL themselves
- scheduling: one queue with dedup and a result cache; a family (configs that can share a GATK run) is always sent as
  one group; a node gets groups while it has fewer open groups than workers
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
from evaluator.local import MAX_BATCH_S, SOLO_S, extra_seconds, failed, job_key

log = logging.getLogger(__name__)
LOG_DIR = Path.home() / ".minos" / "fleet" / "logs"
KEEP_ROUNDS = 3                                              # rounds kept on main and nodes (older ones are pruned)


@dataclass
class _Published:
    round: Round
    answer_set: AnswerSet
    bam_url: Optional[str] = None
    bai_url: Optional[str] = None


@dataclass
class _Session:
    host: Host
    proc: Optional[subprocess.Popen] = None
    workers: int = 0
    ready: bool = False
    rounds_ready: set = field(default_factory=set)
    rounds_sent: set = field(default_factory=set)            # (round_id, answer-set version) shipped to the node
    rounds_failed: set = field(default_factory=set)          # round_ids the node could not prepare
    groups: List[set] = field(default_factory=list)          # job_ids per open group (one GATK run each)
    jobs: Dict[str, EvalJob] = field(default_factory=dict)   # outstanding job_id -> job
    done: int = 0                                            # results delivered (dashboard)
    gatk_runs: int = 0                                       # groups handed to it: one GATK run each (dashboard)
    gatk_runs_done: int = 0                                  # ... whose every result came back
    outq: "queue.Queue[Optional[str]]" = field(default_factory=queue.Queue)
    closed: bool = False
    finished: bool = False                                   # its session thread has ended
    rejected: str = ""
    logf: Optional[object] = None                            # the node's stderr, kept per host

    def open_groups(self) -> int:
        return sum(1 for g in self.groups if g)

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
        self._queue: list = []
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
        LOG_DIR.mkdir(parents=True, exist_ok=True)

    # ---- hosts --------------------------------------------------------------------------------------------------------
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
                    log.warning("%s: session error: %s", s.host.name, e)
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
        push(h, [bundle["bam"], bundle["bai"]], f"{h.repo}/cache/selftest/")
        args = h.extra.get("node_args", "") + (f" --workers {int(h.workers)}" if h.workers else "")
        env = " ".join(f"{k}={shlex.quote(str(v))}" for k, v in {"DATA_DIR": h.data_dir, **h.extra.get("env", {})}.items())
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
                s.rejected = f"self-test failed: {msg.get('digest')} != {msg.get('expected')}"
            elif msg.get("code") != self.code:
                s.rejected = f"code differs from main ({msg.get('code')} != {self.code}): run control sync"
            if s.rejected:
                log.error("%s: rejected - %s", s.host.name, s.rejected)
                s.hang_up()
                return
            s.workers, s.ready = int(msg["workers"]), True
            log.info("%s: ready, %d workers", s.host.name, s.workers)
            threading.Thread(target=self._send_rounds, args=(s,), daemon=True).start()
        elif t == "round_ready":
            s.rounds_ready.add(msg["round_id"])
            self._pump()
        elif t == "round_failed":
            log.error("%s: round %s unusable there: %s", s.host.name, msg["round_id"], msg.get("error"))
            s.rounds_failed.add(msg["round_id"])
            self._fail_if_unservable(msg["round_id"])
        elif t == "result":
            cs = codec.callset_from(msg["callset"])
            if cs.round_id in self._ended:
                return                                            # the round is over; its jobs were already freed
            job = s.jobs.pop(cs.job_id, None)
            for g in s.groups:
                g.discard(cs.job_id)
            open_before = len(s.groups)
            s.groups = [g for g in s.groups if g]
            s.gatk_runs_done += open_before - len(s.groups)
            if job is not None:
                s.done += 1
                self._finish(cs, job)
            self._pump()

    def _drop(self, s: _Session, why: str) -> None:
        """Session over: its unfinished jobs go back to the queue at once. Caller holds the lock."""
        if s.jobs:
            log.warning("%s: %s - %d jobs re-queued", s.host.name, why, len(s.jobs))
        for job in s.jobs.values():
            if job_key(job) not in self._done:
                heapq.heappush(self._queue, (-job.priority, next(self._seq), job))
        s.jobs, s.groups, s.ready = {}, [], False
        s.rounds_ready, s.rounds_sent, s.rounds_failed = set(), set(), set()
        if s.proc and s.proc.poll() is None:
            s.hang_up()
        self._cv.notify_all()
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
        if n:
            log.error("round %s: no node can prepare it - %d queued jobs failed", rid, n)
            self._cv.notify_all()

    # ---- rounds -------------------------------------------------------------------------------------------------------
    def publish_round(self, round: Round, answer_set: AnswerSet, bam_url: Optional[str] = None,
                      bai_url: Optional[str] = None) -> None:
        self.refresh_hosts()
        with self._lock:
            self._rounds.pop(round.round_id, None)
            self._rounds[round.round_id] = _Published(round, answer_set, bam_url, bai_url)
            while len(self._rounds) > KEEP_ROUNDS:              # oldest first (dicts keep insertion order)
                self._rounds.pop(next(iter(self._rounds)))
            sessions = [s for s in self._sessions.values() if s.ready]
        for s in sessions:
            threading.Thread(target=self._send_rounds, args=(s,), daemon=True).start()

    def _send_rounds(self, s: _Session) -> None:
        with self._lock:
            rounds = list(self._rounds.items())
        keep = [rid for rid, _ in rounds]
        for rid, p in rounds:
            key = (rid, p.answer_set.version)
            with self._lock:
                if key in s.rounds_sent or not s.ready:
                    continue
                s.rounds_sent.add(key)
                s.rounds_ready.discard(rid)                      # a new version: wait for the node to re-prepare
            try:
                s.send({**round_message(s.host, p.round, p.answer_set, p.bam_url, p.bai_url), "keep": keep})
            except Exception as e:  # noqa: BLE001 - e.g. rsync failed: this node cannot serve the round
                log.error("%s: could not ship round %s: %s", s.host.name, rid, e)
                with self._lock:
                    s.rounds_sent.discard(key)
                    s.rounds_failed.add(rid)
                    self._fail_if_unservable(rid)

    # ---- Executor protocol ----------------------------------------------------------------------------------------
    def set_finish_by(self, ts: Optional[float]) -> None:
        """Do not send a batch that would still be running at `ts` (the final submission). None clears it."""
        with self._lock:
            self._finish_by = ts

    def _room(self) -> float:
        room = MAX_BATCH_S
        if self._finish_by is not None:
            room = min(room, self._finish_by - self._now())
        return room

    def _branch(self, job: EvalJob) -> Optional[str]:
        f = getattr(self.runner, "branch_key", None)
        return f(job.config, job.fidelity) if f else None

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
                    heapq.heappush(self._queue, (-job.priority, next(self._seq), job))
            self._pump()
            self._cv.notify_all()

    def poll(self, timeout: float) -> List[CallSet]:
        deadline = time.time() + timeout
        with self._cv:
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
            for s in self._sessions.values():
                mine = [jid for jid, j in s.jobs.items() if j.round.round_id == round_id]
                if not mine:
                    continue
                for jid in mine:
                    self._inflight.pop(job_key(s.jobs.pop(jid)), None)
                for g in s.groups:
                    g.difference_update(mine)
                s.groups = [g for g in s.groups if g]
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
                    for p in [q[2] for q in self._queue] + jobs[:i]:
                        if p.round.round_id == job.round.round_id and self._branch(p) == bk:
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
            return max(0, free - len(self._queue))

    def stats(self) -> dict:
        """For the public dashboard: nodes are numbered, never named (no IPs)."""
        with self._lock:
            nodes = [dict(node=f"node {i}", ready=s.ready, workers=s.workers, busy=s.open_groups(), jobs=len(s.jobs),
                          done=s.done, gatk_runs=s.gatk_runs, gatk_runs_done=s.gatk_runs_done,
                          rejected=bool(s.rejected))
                     for i, s in enumerate(self._sessions.values(), 1)]
            return {"kind": "fleet", "queued": len(self._queue), "nodes": nodes, "done": sum(n["done"] for n in nodes),
                    "gatk_runs": sum(n["gatk_runs"] for n in nodes),
                    "gatk_runs_done": sum(n["gatk_runs_done"] for n in nodes),
                    "workers": sum(n["workers"] for n in nodes if n["ready"]),
                    "running": sum(n["busy"] for n in nodes if n["ready"])}

    def nodes(self) -> List[dict]:
        with self._lock:
            return [dict(host=n, ready=s.ready, workers=s.workers, open_groups=s.open_groups(), jobs=len(s.jobs),
                         rounds=len(s.rounds_ready), rejected=s.rejected) for n, s in self._sessions.items()]

    def close(self) -> None:
        """Stop every node (closing its stdin ends the agent, which kills its GATK containers)."""
        self._stopping = True
        with self._lock:
            for s in self._sessions.values():
                s.closed = True
                if s.proc and s.proc.poll() is None:
                    s.hang_up()

    # ---- scheduling -------------------------------------------------------------------------------------------------
    def _pump(self) -> None:
        """Hand queued groups to nodes with free workers that have the group's round ready. Caller holds the lock."""
        for s in self._sessions.values():
            while s.ready and s.open_groups() < s.workers and self._queue:
                groups = self._take_groups(s.rounds_ready)
                if not groups:
                    break
                g = groups[0]
                for j in g:
                    s.jobs[j.job_id] = j
                s.groups.append({j.job_id for j in g})
                s.gatk_runs += 1
                if not s.send({"type": "jobs", "groups": [[codec.job_to(j) for j in g]]}):
                    break                                      # the reader sees the session end and re-queues

    def _take_groups(self, rounds_ready: set) -> List[List[EvalJob]]:
        """The top queued job whose round the node has and that can finish in time, plus queued jobs that can
        share its GATK run without pushing the batch past the family cap or the time left."""
        skipped, too_long, group = [], [], None
        room = self._room()
        while self._queue and group is None:
            item = heapq.heappop(self._queue)
            job = item[2]
            if self._skip(job):
                continue
            if job.round.round_id not in rounds_ready:
                skipped.append(item)
                continue
            if SOLO_S > room:
                too_long.append(job)
                continue
            group, seconds, bk = [job], SOLO_S, self._branch(job)
            if bk is not None:
                keep = []
                while self._queue:
                    other_item = heapq.heappop(self._queue)
                    other = other_item[2]
                    if self._skip(other):
                        continue
                    extra = None
                    if (len(group) < self.max_batch and other.round.round_id == job.round.round_id
                            and self._branch(other) == bk):
                        extra = extra_seconds(self.runner, job, other)
                    if extra is not None and seconds + extra <= room:
                        group.append(other)
                        seconds += extra
                    else:
                        keep.append(other_item)
                for it in keep:
                    heapq.heappush(self._queue, it)
        for it in skipped:
            heapq.heappush(self._queue, it)
        for job in too_long:
            heapq.heappush(self._queue, (-job.priority, next(self._seq), job))
        return [group] if group else []

    def _skip(self, job: EvalJob) -> bool:
        k = job_key(job)
        if job.job_id in self._cancelled:
            waiters = self._waiters.pop(k, [])
            if waiters:
                self._inflight[k] = waiters[0]
                self._waiters[k] = waiters[1:]
                heapq.heappush(self._queue, (-waiters[0].priority, next(self._seq), waiters[0]))
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


def round_message(h: Host, r: Round, a: AnswerSet, bam_url: Optional[str], bai_url: Optional[str]) -> dict:
    """Push this round's files to the host and describe the round with the host's own paths."""
    rdir, adir = f"{h.data_dir}/rounds/{r.tag}", f"{h.repo}/cache/fleet/{r.tag}/v{a.version}"   # the node's own cache
    truth = [a.notes["truth_vcf"], a.notes["mutations_vcf"]]
    push(h, truth + [f + ".tbi" for f in truth if Path(f + ".tbi").exists()], adir + "/")
    if not (bam_url and bai_url):                              # replay: main's own copy of the BAM
        push(h, [r.bam, r.bam + ".bai"], rdir + "/")
    ref = f"{h.data_dir}/reference/{r.chrom}/{r.chrom}.fa"
    rnd = replace(r, bam=f"{rdir}/{Path(r.bam).name}", reference=ref, released_dir=None)
    notes = {"truth_vcf": f"{adir}/{Path(a.notes['truth_vcf']).name}",
             "mutations_vcf": f"{adir}/{Path(a.notes['mutations_vcf']).name}", "region": r.region, "reference": ref}
    return {"type": "round", "round": codec.round_to(rnd),
            "answer_set": {"version": a.version, "n_targets": a.n_targets, "notes": notes},
            "bam_url": bam_url, "bai_url": bai_url}
