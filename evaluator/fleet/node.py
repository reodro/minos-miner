"""Node agent: runs on each sub-VPS. Pulls job groups from main, runs GATK and scores the callsets here.

    MINOS_FLEET_TOKEN=... python -m evaluator.fleet.node --main https://MAIN:8707 --cacert fleet.crt

1. self-test: run the bundle's fixed config on its BAM slice; the output digest must equal the expected one
   (same GATK image, patch jar, reference and scoring stack as everywhere else) - otherwise refuse to join
2. register (name, workers = floor(cores x 0.95 - 1))
3. loop: when workers are free, ask main for that many job groups (a family is one group -> one GATK run);
   run them through LocalExecutor (batching, divergence repair, dedup as on a single machine); score every
   finished callset with DecisionCacheScorer against the round's truth files; send back callsets carrying their
   score. VCFs stay on this machine.
4. heartbeat every `heartbeat_s` so main keeps this node's leases.
"""
from __future__ import annotations

import argparse
import logging
import os
import queue
import socket
import threading
import time
from dataclasses import replace
from pathlib import Path
from typing import Dict, List, Optional

from common.settings import SETTINGS, default_workers
from common.types import FULL, AnswerSet, CallSet, Config, Regime, Round
from evaluator.fleet import codec

log = logging.getLogger(__name__)


class NodeAgent:
    def __init__(self, main_url: str, token: Optional[str] = None, name: Optional[str] = None,
                 workers: Optional[int] = None, verify=True, score_threads: int = 2, heartbeat_s: float = 30.0,
                 data_dir: Optional[str] = None, runner=None):
        import httpx
        self.token = token or os.environ.get("MINOS_FLEET_TOKEN") or ""
        if not self.token:
            raise ValueError("node needs the fleet token (MINOS_FLEET_TOKEN)")
        self.name = name or socket.gethostname()
        self.workers = workers or default_workers()
        self.data_dir = Path(data_dir or SETTINGS.data_dir)
        self.http = httpx.Client(base_url=main_url.rstrip("/"), headers={"Authorization": f"Bearer {self.token}"},
                                 verify=verify, timeout=120)
        self.verify, self.heartbeat_s, self.score_threads = verify, heartbeat_s, score_threads
        if runner is None:
            from evaluator.gatk import DEFAULT_PATCH_JAR, GatkRunner
            runner = GatkRunner(threads=1, memory_gb=4, patch_jar=DEFAULT_PATCH_JAR,
                                work_dir=str(self.data_dir / "cache" / "callsets"))
        self.runner = runner
        self.node_id: Optional[str] = None
        self._rounds: Dict[str, tuple] = {}                    # round_id -> (local Round, local AnswerSet)
        self._lease_of: Dict[str, str] = {}                    # job_id -> lease_id
        self._to_score: "queue.Queue[CallSet]" = queue.Queue()
        self._scored: "queue.Queue[CallSet]" = queue.Queue()
        self._lock = threading.Lock()
        self._stop = threading.Event()

    # ---- talking to main ----------------------------------------------------------------------------------------
    def _post(self, path: str, body: dict) -> dict:
        for attempt in range(6):
            try:
                r = self.http.post(path, json=body)
                r.raise_for_status()
                return r.json()
            except Exception as e:  # noqa: BLE001 - main restarting or the network blinking: retry
                if attempt == 5:
                    raise
                log.warning("POST %s failed (%s); retrying", path, e)
                time.sleep(2 * (attempt + 1))
        return {}

    def _download(self, url_or_path: str, dest: Path, direct: bool = False) -> None:
        if dest.exists():
            return
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_name(dest.name + ".part")
        if direct:                                              # presigned URL: no fleet token
            import httpx
            with httpx.stream("GET", url_or_path, timeout=600, follow_redirects=True) as r, open(tmp, "wb") as f:
                r.raise_for_status()
                for chunk in r.iter_bytes(1 << 20):
                    f.write(chunk)
        else:
            with self.http.stream("GET", url_or_path) as r, open(tmp, "wb") as f:
                r.raise_for_status()
                for chunk in r.iter_bytes(1 << 20):
                    f.write(chunk)
        os.replace(tmp, dest)

    def reference_for(self, chrom: str) -> str:
        return str(self.data_dir / "reference" / chrom / f"{chrom}.fa")

    # ---- self-test ----------------------------------------------------------------------------------------------------
    def selftest(self) -> bool:
        from evaluator.gatk import records_digest
        m = self.http.get("/v1/selftest").raise_for_status().json()
        d = self.data_dir / "selftest"
        for name in ("input.bam", "input.bam.bai"):
            self._download(f"/v1/selftest/{name}", d / f"{m['tag']}_{name}")
        chrom, span = m["slice"].split(":")
        a, b = (int(x) for x in span.split("-"))
        rnd = Round(f"selftest-{m['tag']}", chrom, a, b, str(d / f"{m['tag']}_input.bam"), self.reference_for(chrom))
        cs = self.runner.run(rnd, Config("gatk", m["config"]), FULL)
        got = records_digest(Path(cs.vcf)) if cs.ok else f"failed: {cs.error}"
        ok = got == m["expected_sha256"]
        (log.info if ok else log.error)("self-test %s: got %s, expected %s", "passed" if ok else "FAILED", got,
                                        m["expected_sha256"])
        return ok

    # ---- rounds -----------------------------------------------------------------------------------------------------
    def ensure_round(self, round_id: str) -> tuple:
        if round_id in self._rounds:
            return self._rounds[round_id]
        from urllib.parse import quote
        m = self.http.get(f"/v1/rounds/{quote(round_id, safe='')}").raise_for_status().json()
        main_round = codec.round_from(m["round"])
        tag = main_round.tag
        rd = self.data_dir / "rounds" / tag
        for name in ("input.bam", "input.bam.bai"):
            if name in m["urls"]:
                self._download(m["urls"][name], rd / name, direct=True)
            else:
                self._download(f"/v1/files/{quote(round_id, safe='')}/{name}", rd / name)
        ad = self.data_dir / "cache" / "fleet" / tag / f"v{m['version']}"
        for name in m["files"]:
            if name.startswith(("truth", "mutations")):
                self._download(f"/v1/files/{quote(round_id, safe='')}/{name}", ad / name)
        ref = self.reference_for(main_round.chrom)
        local = replace(main_round, bam=str(rd / "input.bam"), reference=ref, released_dir=None)
        aset = AnswerSet(round_id, m["version"], Regime.UNKNOWN, None, [], [], m["n_targets"],
                         notes={"truth_vcf": str(ad / "truth.vcf.gz"), "mutations_vcf": str(ad / "mutations.vcf.gz"),
                                "region": main_round.region, "reference": ref})
        self._rounds[round_id] = (local, aset)
        return local, aset

    # ---- scoring threads --------------------------------------------------------------------------------------------
    def _scorer_loop(self, i: int) -> None:
        from evaluator.decision_cache import DecisionCacheScorer
        scorer = DecisionCacheScorer(work_dir=str(self.data_dir / "cache" / f"dcache_node{i}"))
        while not self._stop.is_set():
            try:
                first = self._to_score.get(timeout=1.0)
            except queue.Empty:
                continue
            batch = [first]
            while len(batch) < 8:
                try:
                    batch.append(self._to_score.get_nowait())
                except queue.Empty:
                    break
            by_round: Dict[str, List[CallSet]] = {}
            for cs in batch:
                by_round.setdefault(cs.round_id, []).append(cs)
            for rid, css in by_round.items():
                _, aset = self._rounds[rid]
                try:
                    results = scorer.score_many(css, aset) if len(css) > 1 else [scorer.score(css[0], aset)]
                    pairs = list(zip(css, results))
                except Exception:  # noqa: BLE001 - one bad callset must not sink the others: one by one
                    pairs = []
                    for cs in css:
                        try:
                            pairs.append((cs, scorer.score(cs, aset)))
                        except Exception as e:  # noqa: BLE001
                            pairs.append((replace(cs, ok=False, error=f"scoring failed on node: {e}"), None))
                for cs, r in pairs:
                    self._scored.put(replace(cs, score=r) if r is not None else cs)

    # ---- heartbeat -------------------------------------------------------------------------------------------------
    def _heartbeat_loop(self) -> None:
        while not self._stop.wait(self.heartbeat_s):
            with self._lock:
                leases = sorted(set(self._lease_of.values()))
            try:
                self._post("/v1/heartbeat", {"node_id": self.node_id, "lease_ids": leases})
            except Exception as e:  # noqa: BLE001
                log.warning("heartbeat failed: %s", e)

    # ---- main loop -------------------------------------------------------------------------------------------------
    def run(self, max_idle_s: Optional[float] = None) -> None:
        """Serve until stopped (or, for tests, until nothing has arrived for `max_idle_s`)."""
        from evaluator.local import LocalExecutor
        if not self.selftest():
            raise SystemExit("self-test failed: this node must not evaluate configs")
        self.node_id = self._post("/v1/register", {"name": self.name, "workers": self.workers})["node_id"]
        executor = LocalExecutor(self.runner, workers=self.workers)
        threads = [threading.Thread(target=self._scorer_loop, args=(i,), daemon=True) for i in range(self.score_threads)]
        threads.append(threading.Thread(target=self._heartbeat_loop, daemon=True))
        for t in threads:
            t.start()
        last_activity, backoff = time.time(), 0.0
        try:
            while not self._stop.is_set():
                free = executor.idle_slots()
                if free > 0 and time.time() >= backoff:
                    w = self._post("/v1/work", {"node_id": self.node_id, "slots": free})
                    if w["groups"]:
                        last_activity = time.time()
                        for group in w["groups"]:
                            jobs = [codec.job_from(d) for d in group]
                            local, _ = self.ensure_round(jobs[0].round.round_id)
                            jobs = [replace(j, round=local) for j in jobs]
                            with self._lock:
                                for j in jobs:
                                    self._lease_of[j.job_id] = w["lease_id"]
                            executor.submit(jobs)
                    else:
                        backoff = time.time() + 1.0             # nothing queued on main: ask again in a second
                for cs in executor.poll(0.5):
                    if cs.ok:
                        self._to_score.put(cs)
                    else:
                        self._scored.put(cs)
                self._send_results()
                if max_idle_s is not None and time.time() - last_activity > max_idle_s and not self._lease_of:
                    break
        finally:
            self._stop.set()

    def _send_results(self) -> None:
        out: Dict[str, List[dict]] = {}
        while True:
            try:
                cs = self._scored.get_nowait()
            except queue.Empty:
                break
            with self._lock:
                lease = self._lease_of.pop(cs.job_id, None)
            if lease:
                out.setdefault(lease, []).append(codec.callset_to(cs))
        for lease, css in out.items():
            self._post("/v1/results", {"node_id": self.node_id, "lease_id": lease, "callsets": css})

    def stop(self) -> None:
        self._stop.set()


def main() -> None:
    ap = argparse.ArgumentParser(prog="python -m evaluator.fleet.node")
    ap.add_argument("--main", required=True, help="https://MAIN_HOST:8707")
    ap.add_argument("--name")
    ap.add_argument("--workers", type=int)
    ap.add_argument("--cacert", help="main's certificate (self-signed); omit to use the system CAs")
    ap.add_argument("--insecure-http", action="store_true", help="allow http:// (localhost tests only)")
    ap.add_argument("--score-threads", type=int, default=2)
    ap.add_argument("--max-idle-s", type=float, help="exit after this long without work (tests)")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s node: %(message)s")
    if a.main.startswith("http://") and not a.insecure_http:
        raise SystemExit("refusing to send the fleet token over plain http (use https, or --insecure-http on localhost)")
    NodeAgent(a.main, name=a.name, workers=a.workers, verify=a.cacert or True,
              score_threads=a.score_threads).run(max_idle_s=a.max_idle_s)


if __name__ == "__main__":
    main()
