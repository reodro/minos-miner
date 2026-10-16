"""True score of the replay's current best, for the dashboard.

The search scores against its own answer set. That number is not a place on the published field.
This process scores the current best callset on the released truth and writes dashboard/static/true_live.json.
It does not change what the search proposes.
"""
from __future__ import annotations

import json
import os
import shlex
import subprocess
import time
from pathlib import Path

from common.settings import SETTINGS
from common.types import FULL, CallSet, Round

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "static" / "true_live.json"
LIVE = Path(SETTINGS.status_dir) / "live.json"
WORK = Path(SETTINGS.cache_dir) / "true_live"


def _write(row: dict) -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.with_suffix(".json.part")
    tmp.write_text(json.dumps(row))
    os.replace(tmp, OUT)
    print(row.get("state"), row.get("tag"), row.get("rank"), flush=True)


def _read_out() -> dict:
    try:
        return json.loads(OUT.read_text())
    except (OSError, ValueError):
        return {}


def _indexed(path: str) -> str:
    p = Path(path)
    if not Path(str(p) + ".tbi").exists() and not Path(str(p) + ".csi").exists():
        subprocess.run(["bcftools", "index", "-f", str(p)], check=True, timeout=120, capture_output=True)
    return str(p)


def _pull(tag: str, key: str) -> str | None:
    """The best config's VCF, from this machine or the fleet node that ran it."""
    local = Path(SETTINGS.cache_dir) / "callsets" / tag / f"{key}_{FULL}.vcf.gz"
    if local.exists():
        return _indexed(str(local))
    from evaluator.fleet.inventory import ensure_agent, load
    hosts = load()
    ensure_agent(hosts)
    dest = WORK / tag / f"{key}_{FULL}.vcf.gz"
    for h in hosts:
        remote = f"{h.repo}/cache/callsets/{tag}/{key}_{FULL}.vcf.gz"
        probe = subprocess.run(h.ssh(f"test -f {shlex.quote(remote)}"), capture_output=True, timeout=30)
        if probe.returncode != 0:
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        if h.local:
            subprocess.run(["cp", "-a", remote, str(dest)], check=True, timeout=60)
        else:
            ssh = " ".join(["ssh", *h.ssh_opts()])
            subprocess.run(
                ["rsync", "-a", "--partial", "-e", ssh, f"{h.user}@{h.host}:{remote}", str(dest)],
                check=True, timeout=180, capture_output=True)
        if dest.exists():
            return _indexed(str(dest))
    return None


def _score(live: dict, vcf: str, key: str) -> float:
    from answer_set.oracle import OracleProvider
    from backtest.data import round_dir
    from evaluator.gatk import HappyScorer
    rnd_live = live["round"]
    chrom, span = rnd_live["region"].split(":")
    start_s, end_s = span.split("-")
    tag = rnd_live["tag"]
    rnd = Round(rnd_live["id"], chrom, int(start_s), int(end_s),
                os.path.join(round_dir(tag), "input.bam"),
                SETTINGS.reference_fasta(chrom),
                released_dir=round_dir(tag))
    oracle = OracleProvider().build(rnd)
    res = HappyScorer().score(CallSet(rnd.round_id, key, vcf, "0", 0.0), oracle)
    if res.score is None:
        raise RuntimeError("true score missing")
    return float(res.score)


def once() -> None:
    from dashboard.field import place, published_scores
    try:
        live = json.loads(LIVE.read_text())
    except (OSError, ValueError):
        return
    rnd = live.get("round") or {}
    search = live.get("search") or {}
    tag, key = rnd.get("tag"), search.get("best_config")
    if rnd.get("kind") != "replay" or not tag or not key:
        return
    prev = _read_out()
    if prev.get("tag") == tag and prev.get("config") == key and prev.get("state") == "ready":
        return
    _write({"tag": tag, "config": key, "state": "scoring", "local": search.get("best")})
    try:
        vcf = _pull(tag, key)
        if not vcf:
            _write({"tag": tag, "config": key, "state": "missing", "local": search.get("best")})
            return
        score = _score(live, vcf, key)
        field = published_scores(str(rnd["id"]))
        row = place(score, field, replay=True)
        _write({"tag": tag, "config": key, "state": "ready", "true_score": round(score, 4),
                "local": search.get("best"), "rank": row["rank"], "field": row["field"],
                "field_best": row["field_best"], "gap_to_best": row["gap_to_best"], "at": time.time()})
    except Exception:  # noqa: BLE001 - the board keeps the last ready row; the message stays off the public file
        if prev.get("state") == "ready" and prev.get("tag") == tag and prev.get("config") == key:
            return
        _write({"tag": tag, "config": key, "state": "error"})


def main() -> None:
    while True:
        once()
        time.sleep(20)


if __name__ == "__main__":
    main()
