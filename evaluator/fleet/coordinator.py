"""FleetExecutor: the Executor on main when GATK runs on sub-VPSs. Main never runs GATK or hap.py in this mode.

Same protocol as LocalExecutor (submit / poll / cancel / cost_hint(s) / idle_slots), so search and main do not
change. Behind it, an HTTP API (FastAPI) that node agents pull from:

    POST /v1/register     node announces itself (name, workers) after passing its self-test -> node_id
    POST /v1/work         "I have N free slots" -> a lease with up to N job groups (a family is one group)
    POST /v1/results      callsets (with scores computed on the node) for a lease
    POST /v1/heartbeat    keeps the node's leases alive
    GET  /v1/rounds/{id}  round manifest: region, answer-set version/N, where to get the BAM and truth files
    GET  /v1/files/{id}/{name}   files main serves (truth VCFs; the BAM only when no presigned URL is given)
    GET  /v1/selftest     the self-test bundle (a small BAM slice + config + expected digest)

A lease not renewed within `lease_s` is dropped and its unfinished jobs go back to the queue. Every request must
carry `Authorization: Bearer <token>`; serve over HTTPS (certfile/keyfile) whenever it leaves localhost.
"""
from __future__ import annotations

import heapq
import itertools
import json
import logging
import os
import secrets
import threading
import time
import uuid
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from common.settings import SETTINGS
from common.types import AnswerSet, CallSet, EvalJob, Round
from evaluator.fleet import codec
from evaluator.local import failed, job_key

log = logging.getLogger(__name__)
SELFTEST_FILE = Path(__file__).with_name("selftest.json")


@dataclass
class _Node:
    node_id: str
    name: str
    workers: int
    last_seen: float
    busy_groups: int = 0


@dataclass
class _Lease:
    lease_id: str
    node_id: str
    jobs: Dict[str, EvalJob]                  # job_id -> job, removed as results arrive
    groups: List[set]                         # job_ids per group (one GATK run each); a group frees a slot when empty
    expires: float

    def open_groups(self) -> int:
        return sum(1 for g in self.groups if g)


@dataclass
class _RoundAssets:
    round: Round
    version: int
    n_targets: Optional[int]
    files: Dict[str, str] = field(default_factory=dict)       # name -> local path served by main
    urls: Dict[str, str] = field(default_factory=dict)        # name -> URL nodes fetch directly (presigned BAM)


class FleetExecutor:
    def __init__(self, planner, token: Optional[str] = None, host: str = "0.0.0.0", port: int = 8707,
                 lease_s: float = 900.0, max_batch: int = 32, certfile: Optional[str] = None,
                 keyfile: Optional[str] = None):
        """planner: an object with branch_key/variant_cost (a GatkRunner with patch_jar set) used only to group
        families and price them - nothing is run on main."""
        self.runner = planner
        self.token = token or os.environ.get("MINOS_FLEET_TOKEN") or ""
        if not self.token:
            raise ValueError("FleetExecutor needs a token (MINOS_FLEET_TOKEN)")
        self.host, self.port, self.lease_s, self.max_batch = host, port, lease_s, max_batch
        self.certfile, self.keyfile = certfile, keyfile
        self._lock = threading.Lock()
        self._cv = threading.Condition(self._lock)
        self._seq = itertools.count()
        self._queue: list = []                                  # (-priority, seq, job)
        self._inflight: Dict[Tuple, EvalJob] = {}
        self._waiters = defaultdict(list)
        self._done: Dict[Tuple, CallSet] = {}
        self._ready: List[CallSet] = []
        self._cancelled = set()
        self._leases: Dict[str, _Lease] = {}
        self._nodes: Dict[str, _Node] = {}
        self._rounds: Dict[str, _RoundAssets] = {}
        self._server = None
        self._thread: Optional[threading.Thread] = None

    # ---- round assets (main calls this once the answer set exists) ----------------------------------------------
    def publish_round(self, round: Round, answer_set: AnswerSet, bam_url: Optional[str] = None,
                      bai_url: Optional[str] = None) -> None:
        files = {"truth.vcf.gz": answer_set.notes["truth_vcf"], "mutations.vcf.gz": answer_set.notes["mutations_vcf"]}
        for name in ("truth.vcf.gz", "mutations.vcf.gz"):
            if os.path.exists(files[name] + ".tbi"):
                files[name + ".tbi"] = files[name] + ".tbi"
        urls = {}
        if bam_url and bai_url:
            urls = {"input.bam": bam_url, "input.bam.bai": bai_url}
        else:                                                  # replay / tests: main serves its own copy
            files["input.bam"] = round.bam
            files["input.bam.bai"] = round.bam + ".bai"
        with self._lock:
            self._rounds[round.round_id] = _RoundAssets(round, answer_set.version, answer_set.n_targets, files, urls)

    # ---- Executor protocol ----------------------------------------------------------------------------------------
    def _branch(self, job: EvalJob) -> Optional[str]:
        f = getattr(self.runner, "branch_key", None)
        return f(job.config, job.fidelity) if f else None

    def submit(self, jobs: Iterable[EvalJob]) -> None:
        with self._lock:
            for job in jobs:
                k = job_key(job)
                if k in self._done:
                    self._ready.append(_rekey(self._done[k], job.job_id))
                elif k in self._inflight:
                    self._waiters[k].append(job)
                else:
                    self._inflight[k] = job
                    heapq.heappush(self._queue, (-job.priority, next(self._seq), job))
            self._cv.notify_all()

    def poll(self, timeout: float) -> List[CallSet]:
        deadline = time.time() + timeout
        with self._cv:
            self._expire_leases()
            while not self._ready and time.time() < deadline:
                self._cv.wait(max(0.0, deadline - time.time()))
                self._expire_leases()
            out, self._ready = self._ready, []
            return out

    def cancel(self, job_ids: Iterable[str]) -> None:
        with self._lock:
            self._cancelled.update(job_ids)

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
                    peers = [q[2] for q in self._queue] + jobs[:i]
                    for p in peers:
                        if p.round.round_id == job.round.round_id and self._branch(p) == bk:
                            v = self.runner.variant_cost(p.config, job.config, job.fidelity)
                            if v is not None:
                                c = min(c, v)
                out.append(c)
            return out

    def idle_slots(self) -> int:
        with self._lock:
            live = [n for n in self._nodes.values() if time.time() - n.last_seen < self.lease_s]
            return max(0, sum(n.workers for n in live) - sum(n.busy_groups for n in live) - len(self._queue))

    # ---- scheduling -------------------------------------------------------------------------------------------------
    def _take_groups(self, slots: int) -> List[List[EvalJob]]:
        """Up to `slots` groups; a group is the top job plus every queued job that can share its GATK run."""
        groups = []
        while len(groups) < slots and self._queue:
            _, _, job = heapq.heappop(self._queue)
            if self._skip(job):
                continue
            group, bk = [job], self._branch(job)
            if bk is not None:
                keep = []
                while self._queue:
                    item = heapq.heappop(self._queue)
                    other = item[2]
                    if self._skip(other):
                        continue
                    if (len(group) < self.max_batch and other.round.round_id == job.round.round_id
                            and self._branch(other) == bk):
                        group.append(other)
                    else:
                        keep.append(item)
                for item in keep:
                    heapq.heappush(self._queue, item)
            groups.append(group)
        return groups

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
        return k in self._done                                  # finished by a late result from an expired lease

    def _finish(self, cs: CallSet, job: EvalJob) -> None:
        k = job_key(job)
        if k in self._done:
            return
        cs = _rekey(cs, job.job_id)
        self._inflight.pop(k, None)
        if cs.ok:
            self._done[k] = cs
        self._ready.append(cs)
        self._ready += [_rekey(cs, w.job_id) for w in self._waiters.pop(k, [])]

    def _expire_leases(self) -> None:
        now = time.time()
        for lid in [l for l, lease in self._leases.items() if lease.expires < now]:
            lease = self._leases.pop(lid)
            node = self._nodes.get(lease.node_id)
            if node:
                node.busy_groups = max(0, node.busy_groups - lease.open_groups())
            for job in lease.jobs.values():
                if job_key(job) not in self._done:
                    heapq.heappush(self._queue, (-job.priority, next(self._seq), job))
            log.warning("lease %s of node %s expired: %d jobs re-queued", lid, lease.node_id, len(lease.jobs))

    # ---- HTTP API ---------------------------------------------------------------------------------------------------
    def app(self):
        from fastapi import Depends, FastAPI, Header, HTTPException
        from fastapi.responses import FileResponse

        api = FastAPI(title="minos fleet")

        def auth(authorization: str = Header(default="")):
            if not secrets.compare_digest(authorization, f"Bearer {self.token}"):
                raise HTTPException(401, "bad token")

        @api.post("/v1/register", dependencies=[Depends(auth)])
        def register(body: dict):
            node_id = body.get("node_id") or uuid.uuid4().hex[:12]
            with self._lock:
                self._nodes[node_id] = _Node(node_id, body["name"], int(body["workers"]), time.time())
                self._cv.notify_all()
            log.info("node %s (%s) registered with %s workers", node_id, body["name"], body["workers"])
            return {"node_id": node_id, "lease_s": self.lease_s}

        @api.post("/v1/work", dependencies=[Depends(auth)])
        def work(body: dict):
            with self._lock:
                node = self._nodes.get(body["node_id"])
                if node is None:
                    raise HTTPException(409, "register first")
                node.last_seen = time.time()
                self._expire_leases()
                groups = self._take_groups(int(body.get("slots", 1)))
                if not groups:
                    return {"lease_id": None, "groups": []}
                lease = _Lease(uuid.uuid4().hex[:12], node.node_id, {j.job_id: j for g in groups for j in g},
                               [{j.job_id for j in g} for g in groups], time.time() + self.lease_s)
                self._leases[lease.lease_id] = lease
                node.busy_groups += len(groups)
                return {"lease_id": lease.lease_id, "groups": [[codec.job_to(j) for j in g] for g in groups]}

        @api.post("/v1/results", dependencies=[Depends(auth)])
        def results(body: dict):
            with self._cv:
                lease = self._leases.get(body["lease_id"])
                node = self._nodes.get(body["node_id"])
                if node:
                    node.last_seen = time.time()
                for d in body["callsets"]:
                    cs = codec.callset_from(d)
                    job = lease.jobs.pop(cs.job_id, None) if lease else None
                    if lease and job is not None:
                        for g in lease.groups:
                            if cs.job_id in g:
                                g.discard(cs.job_id)
                                if not g and node:           # that GATK run is over: its worker is free again
                                    node.busy_groups = max(0, node.busy_groups - 1)
                    if job is None:                           # late result for an expired lease: still useful
                        job = next((j for j in self._inflight.values() if j.job_id == cs.job_id), None)
                    if job is not None:
                        self._finish(cs, job)
                if lease and not lease.jobs:
                    self._leases.pop(lease.lease_id, None)
                self._cv.notify_all()
            return {"ok": True}

        @api.post("/v1/heartbeat", dependencies=[Depends(auth)])
        def heartbeat(body: dict):
            with self._lock:
                node = self._nodes.get(body["node_id"])
                if node is None:
                    raise HTTPException(409, "register first")
                node.last_seen = time.time()
                for lid in body.get("lease_ids", []):
                    if lid in self._leases:
                        self._leases[lid].expires = time.time() + self.lease_s
            return {"ok": True}

        @api.get("/v1/rounds/{round_id}", dependencies=[Depends(auth)])
        def round_manifest(round_id: str):
            with self._lock:
                a = self._rounds.get(round_id)
            if a is None:
                raise HTTPException(404, "round not published")
            return {"round": codec.round_to(a.round), "version": a.version, "n_targets": a.n_targets,
                    "files": sorted(a.files), "urls": a.urls}

        @api.get("/v1/files/{round_id}/{name}", dependencies=[Depends(auth)])
        def round_file(round_id: str, name: str):
            with self._lock:
                a = self._rounds.get(round_id)
            if a is None or name not in a.files:
                raise HTTPException(404, "no such file")
            return FileResponse(a.files[name])

        @api.get("/v1/selftest", dependencies=[Depends(auth)])
        def selftest():
            return selftest_manifest()

        @api.get("/v1/selftest/{name}", dependencies=[Depends(auth)])
        def selftest_file(name: str):
            path = selftest_files().get(name)
            if path is None:
                raise HTTPException(404, "no such file")
            return FileResponse(path)

        return api

    def start(self) -> "FleetExecutor":
        """Serve the API in a background thread."""
        import uvicorn
        cfg = uvicorn.Config(self.app(), host=self.host, port=self.port, log_level="warning",
                             ssl_certfile=self.certfile, ssl_keyfile=self.keyfile)
        self._server = uvicorn.Server(cfg)
        self._thread = threading.Thread(target=self._server.run, daemon=True)
        self._thread.start()
        for _ in range(100):
            if self._server.started:
                break
            time.sleep(0.05)
        return self

    def stop(self) -> None:
        if self._server:
            self._server.should_exit = True
        if self._thread:
            self._thread.join(timeout=10)

    def nodes(self) -> List[dict]:
        with self._lock:
            return [dict(node_id=n.node_id, name=n.name, workers=n.workers, busy_groups=n.busy_groups,
                         seen_s_ago=round(time.time() - n.last_seen, 1)) for n in self._nodes.values()]


def _rekey(cs: CallSet, job_id: str) -> CallSet:
    from dataclasses import replace
    return replace(cs, job_id=job_id)


# ---- self-test bundle: a small BAM slice of a stored round, a fixed config and the exact expected output ------------
def selftest_spec() -> dict:
    return json.loads(SELFTEST_FILE.read_text())


def selftest_files() -> Dict[str, str]:
    """The bundle's BAM slice (+ index), cut from main's copy of the round the first time it is asked for."""
    spec = selftest_spec()
    d = Path(SETTINGS.cache_dir) / "selftest"
    d.mkdir(parents=True, exist_ok=True)
    bam = d / f"{spec['tag']}_{spec['slice'].replace(':', '_').replace('-', '_')}.bam"
    if not (bam.exists() and Path(str(bam) + ".bai").exists()):
        import pysam
        chrom, span = spec["slice"].split(":")
        a, b = (int(x) for x in span.split("-"))
        src = os.path.join(SETTINGS.rounds_dir, spec["tag"], "input.bam")
        tmp = bam.with_suffix(".part.bam")
        with pysam.AlignmentFile(src) as inp, pysam.AlignmentFile(str(tmp), "wb", template=inp) as out:
            for read in inp.fetch(chrom, max(0, a - 1 - spec["margin"]), b + spec["margin"]):
                out.write(read)
        pysam.index(str(tmp))
        os.replace(tmp, bam)
        os.replace(str(tmp) + ".bai", str(bam) + ".bai")
    return {"input.bam": str(bam), "input.bam.bai": str(bam) + ".bai"}


def selftest_manifest() -> dict:
    spec = selftest_spec()
    return {k: spec[k] for k in ("tag", "slice", "config", "expected_sha256")}
