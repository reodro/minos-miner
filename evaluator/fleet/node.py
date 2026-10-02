"""Node agent: runs on each sub-VPS inside an SSH session that main opens. Talks JSON lines on stdin/stdout.

    main opens:  ssh host 'cd /minos/minos-miner && .venv/bin/python -m evaluator.fleet.node --stdio'

main -> node   {"type": "round", "round": {...}, "answer_set": {...}, "bam_url": ..., "bai_url": ..., "keep": [...]}
               {"type": "jobs", "groups": [[job, ...], ...]}          one group = one GATK run (a family)
               {"type": "abort", "round_id"}                          round over: drop its jobs, kill its containers
node -> main   {"type": "ready", "ok": bool, "digest", "expected", "code", "workers"}   after the self-test
               {"type": "round_ready", "round_id"} | {"type": "round_failed", "round_id", "error"}
               {"type": "result", "callset": {...}}                   carries its score; the VCF stays here
A round is prepared on its own thread (a live BAM download takes minutes); rounds not in `keep` are forgotten.

The node lives exactly as long as the session: when stdin closes (main stopped or the connection dropped), it stops
and kills the GATK containers it started. Logs go to stderr (main keeps them per host).
"""
from __future__ import annotations

import argparse
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
from common.types import AnswerSet, CallSet, Regime
from evaluator.fleet import codec, selftest
from evaluator.fleet.inventory import code_digest

log = logging.getLogger("node")

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
    logging.basicConfig(level=logging.INFO, stream=sys.stderr, format="%(asctime)s %(levelname)s node: %(message)s")
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
        self._rounds: Dict[str, tuple] = {}                  # round_id -> (local Round, AnswerSet)
        self._in: "queue.Queue[dict]" = queue.Queue()
        self._to_score: "queue.Queue[CallSet]" = queue.Queue()
        self._out_lock = threading.Lock()
        self._stop = threading.Event()
        self.score_threads = score_threads or workers
        self._scorers: list = []                             # thread mode (tests): one scorer shared by all threads
        self._pool = None
        if score_dir:
            import multiprocessing
            from concurrent.futures import ProcessPoolExecutor
            self._pool = ProcessPoolExecutor(self.score_threads, mp_context=multiprocessing.get_context("spawn"),
                                             initializer=_proc_init, initargs=(score_dir,))
        self._aborted: set = set()
        self._by_digest: Dict[tuple, tuple] = {}             # (round, answer-set version, records digest) -> result
        self._by_digest_lock = threading.Lock()

    # ---- messages -----------------------------------------------------------------------------------------------
    def send(self, msg: dict) -> None:
        with self._out_lock:
            sys.stdout.write(json.dumps(msg) + "\n")
            sys.stdout.flush()

    def _reader(self) -> None:
        for line in sys.stdin:
            line = line.strip()
            if line:
                try:
                    self._in.put(json.loads(line))
                except json.JSONDecodeError:
                    log.warning("bad line from main: %.200s", line)
        self._stop.set()                                      # main is gone

    # ---- rounds -------------------------------------------------------------------------------------------------
    def prepare_round(self, msg: dict) -> None:
        rnd = codec.round_from(msg["round"])                  # main already set this node's paths
        try:
            for url, path in ((msg.get("bam_url"), rnd.bam), (msg.get("bai_url"), rnd.bam + ".bai")):
                if url and not os.path.exists(path):          # live: fetch the presigned BAM ourselves
                    os.makedirs(os.path.dirname(path), exist_ok=True)
                    subprocess.run(["curl", "-sSfL", "--retry", "4", "-o", path + ".part", url], check=True, timeout=1200)
                    os.replace(path + ".part", path)
            for p in (rnd.bam, rnd.bam + ".bai", rnd.reference, msg["answer_set"]["notes"]["truth_vcf"]):
                if not os.path.exists(p):
                    raise FileNotFoundError(p)
            a = msg["answer_set"]
            aset = AnswerSet(rnd.round_id, a["version"], Regime.UNKNOWN, None, [], [], a["n_targets"], notes=a["notes"])
            self._rounds[rnd.round_id] = (rnd, aset)
            self.send({"type": "round_ready", "round_id": rnd.round_id})
        except Exception as e:  # noqa: BLE001 - reported to main, which then gives this node no work for the round
            log.exception("round %s not usable here", rnd.round_id)
            self.send({"type": "round_failed", "round_id": rnd.round_id, "error": str(e)[:300]})

    def forget_rounds(self, keep: set) -> None:
        """Rounds main no longer keeps: drop their answer sets and scoring caches, and delete their BAM and callsets
        (a node serves for days on a small disk)."""
        for rid in [r for r in self._rounds if r not in keep]:
            rnd, _ = self._rounds.pop(rid)
            for sc in self._scorers:
                drop = getattr(sc, "drop", None)
                if drop:
                    drop(rid)
            with self._by_digest_lock:
                for k in [k for k in self._by_digest if k[0] == rid]:
                    del self._by_digest[k]
            for d in (os.path.dirname(rnd.bam), os.path.join(SETTINGS.cache_dir, "callsets", rnd.tag),
                      os.path.join(SETTINGS.cache_dir, "fleet", rnd.tag), os.path.join(SETTINGS.cache_dir, "dcache", rnd.tag)):
                safe = any(os.path.abspath(d).startswith(os.path.abspath(root) + os.sep)
                           for root in (SETTINGS.data_dir, SETTINGS.cache_dir))     # never outside our own folders
                if os.path.isdir(d) and safe:
                    shutil.rmtree(d, ignore_errors=True)

    def take_jobs(self, msg: dict) -> None:
        for group in msg["groups"]:
            jobs = [codec.job_from(d) for d in group]
            local, _ = self._rounds[jobs[0].round.round_id]
            self.executor.submit([replace(j, round=local) for j in jobs])   # one call: a family stays one batch

    # ---- scoring -------------------------------------------------------------------------------------------------
    def _scorer(self, i: int) -> None:
        """One feeder per scorer process (in thread mode, for tests: a scorer itself); one callset at a time, like a
        validator: GATK's result goes straight to upstream's hap.py path."""
        scorer = None if self._pool else self._scorers[0]
        while not self._stop.is_set():
            try:
                c = self._to_score.get(timeout=1.0)
            except queue.Empty:
                continue
            if c.round_id not in self._rounds:                 # the round was pruned meanwhile: nothing to report
                continue
            aset = self._rounds[c.round_id][1]
            seen = self._seen(c, aset)
            if seen is not None:                               # identical records already scored here this round
                c, r = seen
            elif self._pool is not None:
                try:
                    c, r = self._pool.submit(_proc_score, c, aset).result()
                except Exception as e:  # noqa: BLE001 - a dead scorer process fails its callset, not the node
                    c, r = replace(c, ok=False, error=f"scoring failed on node: {e}"), None
                self._remember(c, aset, r)
            else:
                try:
                    r = scorer.score(c, aset)
                except Exception as e:  # noqa: BLE001
                    c, r = replace(c, ok=False, error=f"scoring failed on node: {e}"), None
                self._remember(c, aset, r)
            if c.round_id not in self._aborted:
                self.send({"type": "result", "callset": codec.callset_to(replace(c, score=r) if r else c)})

    def _digest_key(self, c: CallSet, aset: AnswerSet) -> Optional[tuple]:
        return (c.round_id, aset.version, c.sha256) if c.sha256 else None

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
            while True:
                try:
                    msg = self._in.get_nowait()
                except queue.Empty:
                    break
                if msg["type"] == "round":
                    self.forget_rounds(set(msg.get("keep") or [msg["round"]["round_id"]]) | {msg["round"]["round_id"]})
                    threading.Thread(target=self.prepare_round, args=(msg,), daemon=True).start()
                elif msg["type"] == "jobs":
                    self.take_jobs(msg)
                elif msg["type"] == "abort":                  # main ended the round: stop its work, report nothing
                    self._aborted.add(msg["round_id"])
                    self.executor.end_round(msg["round_id"])
            for cs in self.executor.poll(0.5):
                if cs.round_id in self._aborted:
                    continue
                if cs.ok:
                    self._to_score.put(cs)
                else:
                    self.send({"type": "result", "callset": codec.callset_to(cs)})


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
    workers = a.workers or default_workers()
    prefix = f"minos-node-{uuid.uuid4().hex[:6]}"
    code = code_digest()
    if a.fake:                                # tests: no GATK; env knobs simulate slow runs and bad nodes
        from common.params import ParamRanges
        from evaluator.fake import FakeRunner, FakeScorer
        from evaluator.fleet.inventory import REPO
        runner = FakeRunner(ParamRanges.from_file(str(REPO / "tests" / "fixtures" / "parameter_ranges.json")))
        delay = float(os.environ.get("MINOS_FAKE_DELAY", "0"))
        if delay:
            run = runner.run
            runner.run = lambda *x, **k: (time.sleep(delay), run(*x, **k))[1]
        scorer_factory = lambda i: FakeScorer(runner)  # noqa: E731
        ok = os.environ.get("MINOS_FAKE_SELFTEST") != "fail"
        result = {"ok": ok, "digest": "fake" if ok else "bad", "expected": "fake"}
        code = os.environ.get("MINOS_FAKE_CODE", code)
    else:
        from evaluator.gatk import DEFAULT_PATCH_JAR, GatkRunner, HappyScorer
        runner = GatkRunner(threads=1, memory_gb=4, patch_jar=DEFAULT_PATCH_JAR, container_prefix=prefix,
                            work_dir=os.path.join(SETTINGS.cache_dir, "callsets"))
        scorer_factory = lambda i: HappyScorer()  # noqa: E731
        result = selftest.run(runner, SETTINGS.data_dir, SETTINGS.cache_dir)
    if a.selftest_only:
        print(json.dumps(result)); return
    node = Node(workers, runner, scorer_factory, a.score_threads,
                score_dir=None if a.fake else os.path.join(SETTINGS.cache_dir, "dcache"))
    node.send({"type": "ready", "code": code, "workers": workers, **result})
    if not result["ok"]:
        log.error("self-test failed: %s != %s", result["digest"], result["expected"]); return
    try:
        node.serve()
    finally:
        if node._pool is not None:
            node._pool.shutdown(wait=False, cancel_futures=True)
        if not a.fake:
            reap(prefix)


if __name__ == "__main__":
    main()
