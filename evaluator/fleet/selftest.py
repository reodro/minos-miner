"""Node self-test: a fixed config on a small BAM slice must match selftest.json.

A node is accepted only when both hold: the VCF digest equals expected_sha256, and HappyScorer
(the validator scoring path) returns expected_score. A node that can call variants but cannot
score does not become ready. Main cuts the slice and copies truth/mutations once, then pushes
that bundle to every node.

The test runs when the node is provisioned or synced (control sync), and a pass is recorded with what it
vouches for: the code digest, the patch jar and the expected digest. A session start then reuses the pass
instead of spending 40-130 s on the test again; any change to those three runs it afresh.
"""
from __future__ import annotations

import gzip
import json
import os
import shutil
from pathlib import Path
from typing import Dict, Optional

SPEC_FILE = Path(__file__).with_name("selftest.json")
_TRUTH = (("truth", "truth.vcf.gz"), ("truth_tbi", "truth.vcf.gz.tbi"),
          ("mutations", "mutations.vcf.gz"), ("mutations_tbi", "mutations.vcf.gz.tbi"))
_SCORE_TOL = 1e-4


def spec() -> dict:
    return json.loads(SPEC_FILE.read_text())


def bundle_name() -> str:
    s = spec()
    return f"{s['tag']}_{s['slice'].replace(':', '_').replace('-', '_')}.bam"


def bundle_files(bundle: Dict[str, str]) -> list:
    """Local paths main pushes: slice BAM, its index, and the truth the scorer needs."""
    order = ("bam", "bai", "truth", "truth_tbi", "mutations", "mutations_tbi")
    return [bundle[k] for k in order if bundle.get(k)]


def reject_reason(result: dict) -> str:
    """Why a ready message is not accepted. A score miss is not reported as a digest miss."""
    if result.get("digest") != result.get("expected"):
        return f"self-test failed: {result.get('digest')} != {result.get('expected')}"
    if result.get("error"):
        return f"self-test failed: {result['error']}"
    return f"self-test failed: score {result.get('score')} != {result.get('expected_score')}"


def make_bundle(cache_dir: str, rounds_dir: str) -> Dict[str, str]:
    """On main: slice BAM (+ index) and truth/mutations under cache_dir/selftest.

    The BAM is cut from the stored round when missing. Truth files are copied from that round
    when main has them; scoring fails closed on the node if they are absent.
    """
    s = spec()
    d = Path(cache_dir) / "selftest"
    d.mkdir(parents=True, exist_ok=True)
    bam = d / bundle_name()
    if not (bam.exists() and Path(str(bam) + ".bai").exists()):
        import pysam
        chrom, span = s["slice"].split(":")
        a, b = (int(x) for x in span.split("-"))
        tmp = bam.with_suffix(".part.bam")
        with pysam.AlignmentFile(os.path.join(rounds_dir, s["tag"], "input.bam")) as inp, \
                pysam.AlignmentFile(str(tmp), "wb", template=inp) as out:
            for read in inp.fetch(chrom, max(0, a - 1 - s["margin"]), b + s["margin"]):
                out.write(read)
        pysam.index(str(tmp))
        os.replace(tmp, bam)
        os.replace(str(tmp) + ".bai", str(bam) + ".bai")
    out = {"bam": str(bam), "bai": str(bam) + ".bai"}
    src_dir = Path(rounds_dir) / s["tag"]
    for key, name in _TRUTH:
        src, dest = src_dir / name, d / name
        if src.exists() and (not dest.exists() or dest.stat().st_size != src.stat().st_size):
            shutil.copy2(src, dest)
        if dest.exists():
            out[key] = str(dest)
    return out


def _count_records(path: Path) -> int:
    n = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line and line[0] != "#":
                n += 1
    return n


def _score(runner_cs, s: dict, cache_dir: str, reference: str) -> dict:
    """HappyScorer on the self-test VCF. Any failure is a failed node, not a skipped check."""
    from common.types import AnswerSet, Regime
    from evaluator.gatk import HappyScorer
    d = Path(cache_dir) / "selftest"
    truth, mutations = d / "truth.vcf.gz", d / "mutations.vcf.gz"
    if not truth.exists() or not mutations.exists():
        raise FileNotFoundError("self-test bundle has no truth/mutations (run: python -m evaluator.fleet.control sync)")
    n_targets = _count_records(mutations)
    expected_n = s.get("n_targets")
    if expected_n is not None and int(expected_n) != n_targets:
        raise ValueError(f"mutations count {n_targets} != n_targets {expected_n}")
    if n_targets < 1:
        raise ValueError("self-test mutations VCF has no records")
    answer = AnswerSet(f"selftest-{s['tag']}", 1, Regime.UNKNOWN, None, [], [], n_targets,
                       notes={"truth_vcf": str(truth), "mutations_vcf": str(mutations),
                              "region": s["slice"], "reference": reference})
    scored = HappyScorer().score(runner_cs, answer)
    if scored.score is None:
        raise RuntimeError("HappyScorer returned no score")
    return {"score": float(scored.score), "n_targets": n_targets}


PASS_FILE = "pass.json"


def fingerprint() -> dict:
    """What a pass vouches for: the code, the GATK patch jar and the self-test's expected digest."""
    import hashlib
    from evaluator.fleet.inventory import code_digest
    from evaluator.gatk import DEFAULT_PATCH_JAR
    jar = Path(DEFAULT_PATCH_JAR)
    return {"code": code_digest(), "jar": hashlib.sha256(jar.read_bytes()).hexdigest()[:16] if jar.exists() else None,
            "expected": spec()["expected_sha256"]}


def cached_pass(cache_dir: str) -> Optional[dict]:
    """The recorded pass's result when it still matches this node's fingerprint, else None."""
    try:
        rec = json.loads((Path(cache_dir) / "selftest" / PASS_FILE).read_text())
    except (OSError, ValueError):
        return None
    if rec.get("fingerprint") != fingerprint() or not (rec.get("result") or {}).get("ok"):
        return None
    return {**rec["result"], "cached": True}


def record(cache_dir: str, result: dict) -> None:
    """Keep a passing result (a failure removes any old pass)."""
    path = Path(cache_dir) / "selftest" / PASS_FILE
    if result.get("ok"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"fingerprint": fingerprint(), "result": result}))
    elif path.exists():
        path.unlink()


def run(runner, data_dir: str, cache_dir: str) -> dict:
    """On a node: reproduce the VCF digest, then the validator score. {"ok", "digest", "expected", ...}."""
    from common.types import FULL, Config, Round
    from evaluator.gatk import records_digest
    s = spec()
    base = {"expected": s["expected_sha256"], "expected_score": s.get("expected_score")}
    bam = Path(cache_dir) / "selftest" / bundle_name()
    if not bam.exists():
        return {**base, "ok": False, "digest": "no bundle (run: python -m evaluator.fleet.control sync)"}
    chrom, span = s["slice"].split(":")
    a, b = (int(x) for x in span.split("-"))
    reference = os.path.join(data_dir, "reference", chrom, f"{chrom}.fa")
    rnd = Round(f"selftest-{s['tag']}", chrom, a, b, str(bam), reference)
    cs = runner.run(rnd, Config("gatk", s["config"]), FULL)
    got = records_digest(Path(cs.vcf)) if cs.ok else f"failed: {cs.error}"
    if got != s["expected_sha256"]:
        return {**base, "ok": False, "digest": got}
    try:
        measured = _score(cs, s, cache_dir, reference)
    except Exception as e:  # noqa: BLE001 - a node that cannot score must not become ready
        return {**base, "ok": False, "digest": got, "error": f"scoring failed: {e}"}
    expected = s.get("expected_score")
    if expected is None:
        return {**base, "ok": False, "digest": got, **measured,
                "error": "self-test has no expected_score"}
    ok = abs(float(measured["score"]) - float(expected)) <= _SCORE_TOL
    out = {**base, "ok": ok, "digest": got, "score": measured["score"]}
    if not ok:
        out["error"] = f"score {measured['score']} != {expected}"
    return out
