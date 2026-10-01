"""Node self-test: a fixed config on a small BAM slice must reproduce selftest.json's digest exactly.

It checks the whole stack a node evaluates with - GATK image, patch jar, reference, runner flags - in ~11 s.
Main cuts the slice from its copy of the round once and pushes it to every node.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict

SPEC_FILE = Path(__file__).with_name("selftest.json")


def spec() -> dict:
    return json.loads(SPEC_FILE.read_text())


def bundle_name() -> str:
    s = spec()
    return f"{s['tag']}_{s['slice'].replace(':', '_').replace('-', '_')}.bam"


def make_bundle(cache_dir: str, rounds_dir: str) -> Dict[str, str]:
    """On main: the slice BAM (+ index) under cache_dir/selftest - rebuildable, cut from the stored round when missing."""
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
    return {"bam": str(bam), "bai": str(bam) + ".bai"}


def run(runner, data_dir: str, cache_dir: str) -> dict:
    """On a node: run the bundle (from cache_dir/selftest, pushed by main) on the reference in data_dir, compare
    digests. {"ok", "digest", "expected"}."""
    from common.types import FULL, Config, Round
    from evaluator.gatk import records_digest
    s = spec()
    bam = Path(cache_dir) / "selftest" / bundle_name()
    if not bam.exists():
        return {"ok": False, "digest": "no bundle (run: python -m evaluator.fleet.control sync)", "expected": s["expected_sha256"]}
    chrom, span = s["slice"].split(":")
    a, b = (int(x) for x in span.split("-"))
    rnd = Round(f"selftest-{s['tag']}", chrom, a, b, str(bam), os.path.join(data_dir, "reference", chrom, f"{chrom}.fa"))
    cs = runner.run(rnd, Config("gatk", s["config"]), FULL)
    got = records_digest(Path(cs.vcf)) if cs.ok else f"failed: {cs.error}"
    return {"ok": got == s["expected_sha256"], "digest": got, "expected": s["expected_sha256"]}
