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

from common.settings import SETTINGS, default_score_threads, default_workers
from common.types import AnswerSet, CallSet, Regime
from evaluator.fleet import codec, selftest
from evaluator.fleet.inventory import code_digest

log = logging.getLogger("node")

_PROC_SCORER = None


def _proc_init(work_dir: str) -> None:
    """A scorer process: its own interpreter (own GIL) and its own decision cache. It exits as soon as its parent (the
    node agent, or any other caller) is gone - even after a hard kill - so no scorer is left running on its own."""
    global _PROC_SCORER
    parent = os.getppid()

    def watchdog():
        while os.getppid() == parent:
            time.sleep(2)
        os._exit(1)

    threading.Thread(target=watchdog, daemon=True).start()
    logging.basicConfig(level=logging.INFO, stream=sys.stderr, format="%(asctime)s %(levelname)s node: %(message)s")
    from evaluator.decision_cache import DecisionCacheScorer
    _PROC_SCORER = DecisionCacheScorer(work_dir=os.path.join(work_dir, f"p{os.getpid()}"))


def _proc_score(css: list, aset: AnswerSet, keep: set) -> list:
    """Score callsets of one round in a scorer process: [(callset, ScoreResult or None)], a failed one carrying its
    error. Caches of rounds the node no longer keeps are dropped first."""
    sc = _PROC_SCORER
    for rid in {k[0] for k in list(sc._caches)} - set(keep):
        sc.drop(rid)
    try:
        results = sc.score_many(css, aset) if len(css) > 1 else [sc.score(css[0], aset)]
        return list(zip(css, results))
    except Exception:  # noqa: BLE001 - score one by one so one bad callset cannot sink the batch
        pairs = []
        for c in css:
            try:
                pairs.append((c, sc.score(c, aset)))
            except Exception as e:  # noqa: BLE001
                pairs.append((replace(c, ok=False, error=f"scoring failed on node: {e}"), None))
        return pairs


class Node:
    def __init__(self, workers: int, runner, scorer_factory, score_threads: Optional[int] = None,
                 score_dir: Optional[str] = None):
        """`score_dir` set: score in `score_threads` separate processes (each its own GIL and decision cache) - the
        Python half of scoring needs ~1 core per few dozen GATK workers, more than one interpreter has on a big node.
        Unset (tests): score in threads with `scorer_factory`'s scorer."""
        from evaluator.local import LocalExecutor
        self.workers = workers
        self.runner, self.scorer_factory = runner, scorer_factory
        self.executor = LocalExecutor(runner, workers=workers)
        self._rounds: Dict[str, tuple] = {}                  # round_id -> (local Round, AnswerSet)
        self._in: "queue.Queue[dict]" = queue.Queue()
        self._to_score: "queue.Queue[CallSet]" = queue.Queue()
        self._out_lock = threading.Lock()
        self._stop = threading.Event()
        self.score_threads = score_threads or default_score_threads()
        self._scorers: list = []                             # thread mode (tests): one scorer shared by all threads
        self._pool = None
        if score_dir:
            import multiprocessing
            from concurrent.futures import ProcessPoolExecutor
            self._pool = ProcessPoolExecutor(self.score_threads, mp_context=multiprocessing.get_context("spawn"),
                                             initializer=_proc_init, initargs=(score_dir,))
        self._aborted: set = set()

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
            for d in (os.path.dirname(rnd.bam), os.path.join(SETTINGS.cache_dir, "callsets", rnd.tag),
                      os.path.join(SETTINGS.cache_dir, "fleet", rnd.tag)):
                if os.path.isdir(d) and os.path.abspath(d).startswith(os.path.abspath(SETTINGS.data_dir)):
                    shutil.rmtree(d, ignore_errors=True)

    def take_jobs(self, msg: dict) -> None:
        for group in msg["groups"]:
            jobs = [codec.job_from(d) for d in group]
            local, _ = self._rounds[jobs[0].round.round_id]
            self.executor.submit([replace(j, round=local) for j in jobs])   # one call: a family stays one batch

    # ---- scoring -------------------------------------------------------------------------------------------------
    def _scorer(self, i: int) -> None:
        """One feeder per scorer process (in thread mode, for tests: a scorer itself); batches of up to 8 callsets."""
        scorer = None if self._pool else self._scorers[0]
        while not self._stop.is_set():
            try:
                batch = [self._to_score.get(timeout=1.0)]
            except queue.Empty:
                continue
            while len(batch) < 8:
                try:
                    batch.append(self._to_score.get_nowait())
                except queue.Empty:
                    break
            for rid in {c.round_id for c in batch}:
                css = [c for c in batch if c.round_id == rid]
                if rid not in self._rounds:                     # the round was pruned meanwhile: nothing to report
                    continue
                aset = self._rounds[rid][1]
                if self._pool is not None:
                    try:
                        pairs = self._pool.submit(_proc_score, css, aset, set(self._rounds)).result()
                    except Exception as e:  # noqa: BLE001 - a dead scorer process fails its batch, not the node
                        pairs = [(replace(c, ok=False, error=f"scoring failed on node: {e}"), None) for c in css]
                    for c, r in pairs:
                        if c.round_id not in self._aborted:
                            self.send({"type": "result", "callset": codec.callset_to(replace(c, score=r) if r else c)})
                    continue
                many = getattr(scorer, "score_many", None)
                try:
                    results = many(css, aset) if many and len(css) > 1 else [scorer.score(c, aset) for c in css]
                    pairs = list(zip(css, results))
                except Exception:  # noqa: BLE001 - score one by one so one bad callset cannot sink the batch
                    pairs = []
                    for c in css:
                        try:
                            pairs.append((c, scorer.score(c, aset)))
                        except Exception as e:  # noqa: BLE001
                            pairs.append((replace(c, ok=False, error=f"scoring failed on node: {e}"), None))
                for c, r in pairs:
                    if c.round_id not in self._aborted:
                        self.send({"type": "result", "callset": codec.callset_to(replace(c, score=r) if r else c)})

    # ---- main loop -------------------------------------------------------------------------------------------------
    def serve(self) -> None:
        threading.Thread(target=self._reader, daemon=True).start()
        if self._pool is None:
            self._scorers.append(self.scorer_factory(0))        # shared: a cluster scored once is never re-run here
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
    """Kill every GATK container this node started (names carry the node's prefix)."""
    ids = subprocess.run(["docker", "ps", "-q", "--filter", f"name={prefix}-"], capture_output=True, text=True).stdout.split()
    if ids:
        subprocess.run(["docker", "kill", *ids], capture_output=True)


def main() -> None:
    ap = argparse.ArgumentParser(prog="python -m evaluator.fleet.node")
    ap.add_argument("--stdio", action="store_true", help="serve main over stdin/stdout (how main starts it over ssh)")
    ap.add_argument("--selftest-only", action="store_true", help="run the self-test, print the result, exit")
    ap.add_argument("--workers", type=int)
    ap.add_argument("--score-threads", type=int, help="scorer processes; default: ceil(cores * 0.1)")
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
        from evaluator.decision_cache import DecisionCacheScorer
        from evaluator.gatk import DEFAULT_PATCH_JAR, GatkRunner
        runner = GatkRunner(threads=1, memory_gb=4, patch_jar=DEFAULT_PATCH_JAR, container_prefix=prefix,
                            work_dir=os.path.join(SETTINGS.cache_dir, "callsets"))
        scorer_factory = lambda i: DecisionCacheScorer(work_dir=os.path.join(SETTINGS.cache_dir, "dcache_node"))  # noqa: E731
        result = selftest.run(runner, SETTINGS.data_dir)
    if a.selftest_only:
        print(json.dumps(result)); return
    node = Node(workers, runner, scorer_factory, a.score_threads,
                score_dir=None if a.fake else os.path.join(SETTINGS.cache_dir, "dcache_node"))
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
