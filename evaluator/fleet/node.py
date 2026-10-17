"""Node agent: runs on each sub-VPS inside an SSH session that main opens. Talks JSON lines on stdin/stdout.

    main opens:  ssh host 'cd /minos/minos-miner && .venv/bin/python -m evaluator.fleet.node --stdio'

main -> node   {"type": "round", "round": {...}, "bam_url": ..., "bai_url": ..., "keep": [...]}
                   the round itself: GATK can start once the BAM is here (a live BAM is fetched from the presigned
                   URL; a reference chromosome this node lacks is downloaded)
               {"type": "answer_set", "round_id", "answer_set": {"version", "n_targets", "notes"}}
                   its truth files (main pushed them first): scoring starts, with the results held so far
               {"type": "jobs", "groups": [[job, ...], ...], "limits"}  one group = one GATK run (a family); limits =
                                                                       main's batch limits in force, packed by them
               {"type": "abort", "round_id"}                          round over: drop its jobs, kill its containers
               {"type": "abort_jobs", "job_ids": [...]}               a duplicate run elsewhere delivered first:
                                                                      kill these runs, report nothing for them
node -> main   {"type": "ready", "ok": bool, "digest", "expected", "code", "workers"}   after the self-test
               {"type": "round_ready", "round_id"} | {"type": "round_failed", "round_id", "error"}
               {"type": "result", "callset": {...}}                   carries its score; the VCF stays here
A round is prepared on its own thread (a live BAM download takes minutes); rounds not in `keep` are forgotten.

The node lives exactly as long as the session: when stdin closes (main stopped or the connection dropped), it stops
and kills the GATK containers it started. Logs go to stderr (main keeps them per host).
"""
from __future__ import annotations

import argparse
import itertools
import json
import logging
import os
import queue
import shutil
import subprocess
import sys
import threading
import time
import uuid
from dataclasses import replace
from typing import Dict, Optional

from common.settings import SETTINGS, default_workers
from common.types import AnswerSet, CallSet, Regime, Round
from evaluator.fleet import codec, selftest
from evaluator.fleet.inventory import REPO, code_digest
from evaluator.pool import _proc_init, _proc_score_unit

log = logging.getLogger("node")
WATCHDOG_S = 120.0            # main pings every 20 s (coordinator.PING_S): this long without a line, main is gone
AGENTS = os.path.join(os.path.expanduser("~"), ".minos", "fleet", "agents")   # <prefix>.pid of every live agent


def register(prefix: str) -> None:
    """Record this agent's container prefix and pid, after killing the GATK jobs of agents that died without cleaning
    up (killed, or the node rebooted): their prefixes are on file with a pid that no longer runs. A live agent's jobs
    are left alone - another main may be using this node."""
    os.makedirs(AGENTS, exist_ok=True)
    for name in os.listdir(AGENTS):
        if not name.endswith(".pid"):
            continue
        path = os.path.join(AGENTS, name)
        try:
            with open(path) as f:
                pid = int(f.read().strip() or 0)
        except (OSError, ValueError):
            pid = 0
        alive = False
        if pid > 0:
            try:
                os.kill(pid, 0)
                alive = True
            except ProcessLookupError:
                alive = False
            except PermissionError:
                alive = True
        if not alive:
            log.warning("agent %s (pid %d) is gone: killing its leftover GATK jobs", name[:-4], pid)
            reap(name[:-4])
            try:
                os.remove(path)
            except OSError:
                pass
    with open(os.path.join(AGENTS, f"{prefix}.pid"), "w") as f:
        f.write(str(os.getpid()))


def unregister(prefix: str) -> None:
    try:
        os.remove(os.path.join(AGENTS, f"{prefix}.pid"))
    except OSError:
        pass



def _file_id(path: str) -> str:
    """sha1 of a (small) file, or 'missing': identifies a truth file by content."""
    import hashlib
    try:
        with open(path, "rb") as f:
            return hashlib.sha1(f.read()).hexdigest()[:16]
    except OSError:
        return "missing"


def fetch_reference(chrom: str, data_dir: str) -> None:
    """A chromosome this node has not got: its reference files (FASTA, index, dict, SDF) from the platform's
    endpoint, as setup fetches them (scripts/fetch_genome.sh skips files already present)."""
    log.info("fetching the reference for %s", chrom)
    env = {**os.environ, "CHROMS": chrom, "DATA_DIR": data_dir}
    subprocess.run(["bash", str(REPO / "scripts" / "fetch_genome.sh")], env=env, check=True, timeout=1800,
                   capture_output=True, text=True)


class Node:
    def __init__(self, workers: int, runner, scorer_factory, score_threads: Optional[int] = None,
                 score_dir: Optional[str] = None):
        """One scorer per GATK worker (`score_threads` overrides): main keeps a slot busy until its result is scored, so
        a worker's slot covers GATK and then scoring, and every finished GATK result is scored at once.
        `score_dir` set: the scorers are separate processes (each its own GIL) scoring with upstream's hap.py path.
        Unset (tests): threads sharing `scorer_factory`'s scorer."""
        from evaluator.local import LocalExecutor
        self.workers = workers
        self.runner, self.scorer_factory = runner, scorer_factory
        self.executor = LocalExecutor(runner, workers=workers)
        self._rounds: Dict[str, Round] = {}                  # round_id -> local Round (BAM and reference here)
        self._asets: Dict[str, AnswerSet] = {}               # round_id -> the answer set scoring uses (latest version)
        self._held: Dict[str, list] = {}                     # round_id -> scoring units made before its answer set
        self._held_lock = threading.Lock()
        self._in: "queue.Queue[dict]" = queue.Queue()
        self._to_score: "queue.Queue[list]" = queue.Queue()
        self._out_lock = threading.Lock()
        self._stop = threading.Event()
        from evaluator.pool import scorer_processes
        self.score_threads = scorer_processes(workers) if score_threads is None else score_threads
        self._scorers: list = []                             # thread mode (tests): one scorer shared by all threads
        self._pool = None
        if score_dir:
            import multiprocessing
            from concurrent.futures import ProcessPoolExecutor
            self._pool = ProcessPoolExecutor(self.score_threads, mp_context=multiprocessing.get_context("spawn"),
                                             initializer=_proc_init, initargs=(score_dir,))
        self._aborted: set = set()
        self._dropped: set = set()                           # job ids whose results nobody wants (abort_jobs)
        self._shared_data = os.environ.get("MINOS_NODE_SHARED_DATA") == "1"   # main's own data dir (host "local"):
                                                                            # its rounds are not ours to delete
        self._by_digest: Dict[tuple, tuple] = {}             # (round, answer-set version, records digest) -> result
        self._by_digest_lock = threading.Lock()
        self._group_of: Dict[str, int] = {}                  # job id -> the GATK batch (group) it was sent in
        self._group_ids = itertools.count()
        self._heard = time.time()                            # when main last wrote anything (watchdog)
        self._pinged = False                                 # main sends pings (so silence means it is gone)

    # ---- messages -----------------------------------------------------------------------------------------------
    def send(self, msg: dict) -> None:
        with self._out_lock:
            sys.stdout.write(json.dumps(msg) + "\n")
            sys.stdout.flush()

    def _reader(self) -> None:
        for line in sys.stdin:
            self._heard = time.time()
            line = line.strip()
            if line:
                try:
                    self._in.put(json.loads(line))
                except json.JSONDecodeError:
                    log.warning("bad line from main: %.200s", line)
        self._stop.set()                                      # main is gone

    def _silent(self) -> bool:
        """Main has written nothing (not even its ping, every PING_S) for WATCHDOG_S: it is gone even though the
        ssh channel is still open. Before main's first ping it is never called silent (an older main sends none)."""
        return self._pinged and time.time() - self._heard > WATCHDOG_S

    # ---- rounds -------------------------------------------------------------------------------------------------
    def prepare_round(self, msg: dict) -> None:
        rnd = codec.round_from(msg["round"])                  # main already set this node's paths
        try:
            for url, path in ((msg.get("bam_url"), rnd.bam), (msg.get("bai_url"), rnd.bam + ".bai")):
                if url and not os.path.exists(path):          # live: fetch the presigned BAM ourselves
                    os.makedirs(os.path.dirname(path), exist_ok=True)
                    subprocess.run(["curl", "-sSfL", "--retry", "4", "-o", path + ".part", url], check=True, timeout=1200)
                    os.replace(path + ".part", path)
            if not os.path.exists(rnd.reference):
                fetch_reference(rnd.chrom, SETTINGS.data_dir)
            for p in (rnd.bam, rnd.bam + ".bai", rnd.reference):
                if not os.path.exists(p):
                    raise FileNotFoundError(p)
            self._rounds[rnd.round_id] = rnd
            self.send({"type": "round_ready", "round_id": rnd.round_id})
        except Exception as e:  # noqa: BLE001 - reported to main, which then gives this node no work for the round
            log.exception("round %s not usable here", rnd.round_id)
            self.send({"type": "round_failed", "round_id": rnd.round_id, "error": str(e)[:300]})

    def set_answer_set(self, msg: dict) -> None:
        """The round's truth (a first version, or a newer one): scoring uses it from now on, starting with the units
        held since GATK ran ahead of it."""
        rid, a = msg["round_id"], msg["answer_set"]
        notes = dict(a["notes"])
        try:
            for p in (notes["truth_vcf"], notes.get("mutations_vcf")):
                if p and not os.path.exists(p):
                    raise FileNotFoundError(p)
        except Exception as e:  # noqa: BLE001 - main then stops sending this round here
            log.exception("answer set of %s not usable here", rid)
            self.send({"type": "round_failed", "round_id": rid, "error": str(e)[:300]})
            return
        notes["truth_id"] = _file_id(notes["truth_vcf"]) + "/" + _file_id(notes.get("mutations_vcf", ""))
        aset = AnswerSet(rid, a["version"], Regime.UNKNOWN, None, [], [], a["n_targets"], notes=notes)
        with self._held_lock:
            self._asets[rid] = aset
            held = self._held.pop(rid, [])
        for unit in held:
            self._to_score.put(unit)
        log.info("round %s: answer set v%d (%d held units to score)", rid, aset.version, len(held))

    def forget_rounds(self, keep: set) -> None:
        """Rounds main no longer keeps: drop their answer sets and scoring caches, and delete their BAM and callsets
        (a node serves for days on a small disk)."""
        for rid in [r for r in set(self._rounds) | set(self._asets) | set(self._held) if r not in keep]:
            rnd = self._rounds.pop(rid, None)
            with self._held_lock:
                self._asets.pop(rid, None)
                held = self._held.pop(rid, [])
            if rid not in self._aborted:                      # main may still wait for these: report, never drop
                for c in (c for unit in held for c in unit if c.job_id not in self._dropped):
                    self.send({"type": "result", "callset": codec.callset_to(
                        replace(c, ok=False, error="round forgotten on this node before its answer set came"))})
            for sc in self._scorers:
                drop = getattr(sc, "drop", None)
                if drop:
                    drop(rid)
            with self._by_digest_lock:
                for k in [k for k in self._by_digest if k[0] == rid]:
                    del self._by_digest[k]
            tag = Round.tag_of(rid)
            self._dropped = {j for j in self._dropped if not j.startswith(tag + "-")}
            if rnd is None:
                continue
            dirs = [os.path.join(SETTINGS.cache_dir, "callsets", rnd.tag), os.path.join(SETTINGS.cache_dir, "fleet", rnd.tag),
                    os.path.join(SETTINGS.cache_dir, "dcache", rnd.tag)]
            if not self._shared_data:
                dirs.append(os.path.dirname(rnd.bam))
            for d in dirs:
                safe = any(os.path.abspath(d).startswith(os.path.abspath(root) + os.sep)
                           for root in (SETTINGS.data_dir, SETTINGS.cache_dir))     # never outside our own folders
                if os.path.isdir(d) and safe:
                    shutil.rmtree(d, ignore_errors=True)

    def take_jobs(self, msg: dict) -> None:
        limits = msg.get("limits")                            # main's batch limits in force (settings on main)
        if isinstance(limits, dict):
            self.executor.limits = {k: float(v) for k, v in limits.items()
                                    if k in ("small", "big") and isinstance(v, (int, float)) and v >= 1.0}
        for group in msg["groups"]:
            jobs = [codec.job_from(d) for d in group]
            # A job main sends again (it reclaimed a run this node had been told to drop) keeps its id: it is wanted
            # now. Left in _dropped, its result was discarded here and main waited for it forever (2026-10-17).
            self._dropped.difference_update(j.job_id for j in jobs)
            local = self._rounds.get(jobs[0].round.round_id)
            if local is None:                                 # never sent before round_ready; a safe answer anyway
                for j in jobs:
                    self.send({"type": "result", "callset": codec.callset_to(
                        CallSet(j.round.round_id, j.config.key, "", "", 0.0, ok=False, error="round not prepared here",
                                job_id=j.job_id, fidelity=j.fidelity))})
                continue
            token = next(self._group_ids)
            for j in jobs:
                self._group_of[j.job_id] = token
            self.executor.submit([replace(j, round=local) for j in jobs])   # one call: a family stays one batch

    def abort_jobs(self, job_ids: list) -> None:
        """Their results were delivered by a duplicate run elsewhere: kill the runs here, report nothing for them."""
        ids = set(job_ids)
        self._dropped.update(ids)
        for jid in ids:
            self._group_of.pop(jid, None)
        self.executor.cancel_running(ids)

    # ---- scoring -------------------------------------------------------------------------------------------------
    def _scorer(self, i: int) -> None:
        """One feeder per scorer process (in thread mode, for tests: a scorer itself). It takes one GATK batch's
        callsets at a time and scores them one after another, so a worker's slot never runs more than one scoring:
        no burst of parallel hap.py runs when a family finishes (with --threads 4 each, they oversubscribed big nodes).
        A unit whose round has no answer set yet is held until it arrives."""
        scorer = None if self._pool else self._scorers[0]
        while not self._stop.is_set():
            try:
                unit = self._to_score.get(timeout=1.0)
            except queue.Empty:
                continue
            unit = [c for c in unit if c.job_id not in self._dropped and c.round_id not in self._aborted]
            gone = [c for c in unit if c.round_id not in self._rounds]
            for c in gone:                                     # main still waits for these: never drop them silently
                self.send({"type": "result", "callset": codec.callset_to(
                    replace(c, ok=False, error="round forgotten on this node before scoring"))})
            unit = [c for c in unit if c.round_id in self._rounds]
            if not unit:
                continue
            rid = unit[0].round_id
            with self._held_lock:
                aset = self._asets.get(rid)
                if aset is None:
                    self._held.setdefault(rid, []).append(unit)
                    continue
            done: Dict[int, tuple] = {}
            todo, first_of = [], {}
            for n, c in enumerate(unit):
                seen = self._seen(c, aset)
                if seen is not None:                           # identical records already scored here this round
                    done[n] = seen
                elif c.sha256 and c.sha256 in first_of:        # identical records twice in this unit: score once
                    continue
                else:
                    if c.sha256:
                        first_of[c.sha256] = n
                    todo.append(n)
            if todo:
                css = [unit[n] for n in todo]
                if self._pool is not None:
                    try:
                        pairs = self._pool.submit(_proc_score_unit, css, aset).result()
                    except Exception as e:  # noqa: BLE001 - a dead scorer process fails its callsets, not the node
                        pairs = [(replace(c, ok=False, error=f"scoring failed on node: {e}"), None) for c in css]
                else:
                    pairs = []
                    for c in css:
                        try:
                            pairs.append((c, scorer.score(c, aset)))
                        except Exception as e:  # noqa: BLE001
                            pairs.append((replace(c, ok=False, error=f"scoring failed on node: {e}"), None))
                for n, (c, r) in zip(todo, pairs):
                    done[n] = (c, r)
                    self._remember(c, aset, r)
            for n, c in enumerate(unit):
                if n not in done:                              # a duplicate within the unit: take its twin's result
                    done[n] = self._seen(c, aset) or (replace(c, ok=False, error="scoring failed on node"), None)
                c, r = done[n]
                if c.round_id not in self._aborted and c.job_id not in self._dropped:
                    self.send({"type": "result", "callset": codec.callset_to(replace(c, score=r) if r else c)})

    def _digest_key(self, c: CallSet, aset: AnswerSet) -> Optional[tuple]:
        """Identical records against the same truth CONTENT give the same score (never across two truths that
        happened to land on the same path)."""
        return (c.round_id, aset.version, aset.notes.get("truth_id", ""), c.sha256) if c.sha256 else None

    def _seen(self, c: CallSet, aset: AnswerSet) -> Optional[tuple]:
        """A callset whose records are byte-identical to one already scored here gets that score (or that failure):
        hap.py sees only the records, so the result is the same - no second hap.py run."""
        k = self._digest_key(c, aset)
        with self._by_digest_lock:
            hit = self._by_digest.get(k) if k else None
        if hit is None:
            return None
        r, error = hit
        if r is None:
            return replace(c, ok=False, error=error), None
        return c, replace(r, config_key=c.config_key, job_id=c.job_id, fidelity=c.fidelity)

    def _remember(self, c: CallSet, aset: AnswerSet, r) -> None:
        k = self._digest_key(c, aset)
        if k:
            with self._by_digest_lock:
                self._by_digest[k] = (r, c.error)

    # ---- main loop -------------------------------------------------------------------------------------------------
    def serve(self) -> None:
        threading.Thread(target=self._reader, daemon=True).start()
        if self._pool is None:
            self._scorers.append(self.scorer_factory(0))
        for i in range(self.score_threads):
            threading.Thread(target=self._scorer, args=(i,), daemon=True).start()
        while not self._stop.is_set():
            if self._silent():
                log.error("main silent for %.0f s: stopping (its GATK jobs are killed)", time.time() - self._heard)
                self._stop.set()
                break
            while True:
                try:
                    msg = self._in.get_nowait()
                except queue.Empty:
                    break
                if msg["type"] == "ping":
                    self._pinged = True
                elif msg["type"] == "round":
                    self._aborted.discard(msg["round"]["round_id"])   # main reopened a round it had ended
                    self.forget_rounds(set(msg.get("keep") or [msg["round"]["round_id"]]) | {msg["round"]["round_id"]})
                    threading.Thread(target=self.prepare_round, args=(msg,), daemon=True).start()
                elif msg["type"] == "answer_set":
                    self.set_answer_set(msg)
                elif msg["type"] == "jobs":
                    self.take_jobs(msg)
                elif msg["type"] == "abort":                  # main ended the round: stop its work, report nothing
                    self._aborted.add(msg["round_id"])
                    self.executor.end_round(msg["round_id"])
                elif msg["type"] == "abort_jobs":
                    self.abort_jobs(msg["job_ids"])
            units: Dict[object, list] = {}
            for cs in self.executor.poll(0.5):
                token = self._group_of.pop(cs.job_id, None)
                if cs.round_id in self._aborted or cs.job_id in self._dropped:
                    continue
                if cs.ok:                                      # one scoring unit per GATK batch delivered together
                    units.setdefault(token if token is not None else f"solo:{cs.job_id}", []).append(cs)
                else:
                    self.send({"type": "result", "callset": codec.callset_to(cs)})
            for unit in units.values():
                self._to_score.put(unit)


def reap(prefix: str) -> None:
    """Kill every GATK job this node started (names carry the node's prefix)."""
    from evaluator.sandbox import kill_prefix
    kill_prefix(f"{prefix}-")
    ids = subprocess.run(["docker", "ps", "-q", "--filter", f"name={prefix}-"], capture_output=True, text=True).stdout.split()
    if ids:
        subprocess.run(["docker", "kill", *ids], capture_output=True)


def main() -> None:
    ap = argparse.ArgumentParser(prog="python -m evaluator.fleet.node")
    ap.add_argument("--stdio", action="store_true", help="serve main over stdin/stdout (how main starts it over ssh)")
    ap.add_argument("--selftest-only", action="store_true", help="run the self-test, print the result, exit")
    ap.add_argument("--workers", type=int)
    ap.add_argument("--score-threads", type=int, help="scorer processes; default: one per GATK worker")
    ap.add_argument("--fake", action="store_true", help="tests only: FakeRunner/FakeScorer, no GATK, no self-test")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, stream=sys.stderr, format="%(asctime)s %(levelname)s node: %(message)s")
    from common.limits import raise_open_files
    raise_open_files()                        # hundreds of workers hold more files than the usual soft limit of 1024
    workers = a.workers or default_workers()
    prefix = f"minos-node-{uuid.uuid4().hex[:6]}"
    code = code_digest()
    if a.fake:                                # tests: no GATK; env knobs simulate slow runs and bad nodes
        from common.params import ParamRanges
        from evaluator.fake import FakeRunner, FakeScorer
        runner = FakeRunner(ParamRanges.from_file(str(REPO / "tests" / "fixtures" / "parameter_ranges.json")))
        delay = float(os.environ.get("MINOS_FAKE_DELAY", "0"))
        if delay:
            run = runner.run
            runner.run = lambda *x, **k: (time.sleep(delay), replace(run(*x, **k), runtime_s=delay))[1]
        scorer_factory = lambda i: FakeScorer(runner)  # noqa: E731
        ok = os.environ.get("MINOS_FAKE_SELFTEST") != "fail"
        result = {"ok": ok, "digest": "fake" if ok else "bad", "expected": "fake"}
        code = os.environ.get("MINOS_FAKE_CODE", code)
    else:
        from evaluator.gatk import DEFAULT_PATCH_JAR, GatkRunner, HappyScorer
        runner = GatkRunner(threads=1, memory_gb=4, patch_jar=DEFAULT_PATCH_JAR, container_prefix=prefix,
                            work_dir=os.path.join(SETTINGS.cache_dir, "callsets"))
        scorer_factory = lambda i: HappyScorer()  # noqa: E731
        result = None if a.selftest_only else selftest.cached_pass(SETTINGS.cache_dir)   # passed at sync: skip it
        if result is None:
            result = selftest.run(runner, SETTINGS.data_dir, SETTINGS.cache_dir)
            selftest.record(SETTINGS.cache_dir, result)
    if a.selftest_only:
        print(json.dumps(result)); return
    node = Node(workers, runner, scorer_factory, a.score_threads,
                score_dir=None if a.fake else os.path.join(SETTINGS.cache_dir, "dcache"))
    if not a.fake:
        register(prefix)
    import signal
    for sig in (signal.SIGTERM, signal.SIGHUP):               # a killed session or a hang-up: stop and clean up
        signal.signal(sig, lambda *_: node._stop.set())
    node.send({"type": "ready", "code": code, "workers": workers, **result})
    if not result["ok"]:
        log.error("self-test failed: %s != %s", result["digest"], result["expected"])
        if not a.fake:
            unregister(prefix)
        return
    try:
        node.serve()
    finally:
        if node._pool is not None:
            node._pool.shutdown(wait=False, cancel_futures=True)
        if not a.fake:
            reap(prefix)
            unregister(prefix)


if __name__ == "__main__":
    main()
