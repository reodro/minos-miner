"""Released rounds on disk: data/rounds/<tag>/{input.bam,truth.vcf.gz,mutations.vcf.gz,...,round.json,scores.json}.

    python3 -m backtest.data --last 20        # fetch the 20 newest finalized rounds
"""
from __future__ import annotations

import argparse
import json
import os
import urllib.request
from typing import List, Optional

from common.platform import Platform
from common.settings import SETTINGS
from common.types import Round

FILES = ("input.bam", "input.bam.bai", "truth.vcf.gz", "truth.vcf.gz.tbi", "mutations.vcf.gz", "mutations.vcf.gz.tbi")


def round_dir(tag: str) -> str:
    return os.path.join(SETTINGS.rounds_dir, tag)


def fetch(round_id: str, platform: Optional[Platform] = None, with_bam: bool = True) -> str:
    """Download one released round (files, draw record, every score row). Idempotent."""
    p = platform or Platform()
    ver = p.verification(round_id)
    if "reveal" not in ver:
        raise RuntimeError(f"{round_id} not revealed yet")
    tag = round_id[:19].replace("-", "").replace(":", "") + "Z"
    d = round_dir(tag)
    os.makedirs(d, exist_ok=True)
    for f in ver["reveal"]["files"]:
        if f["name"] not in FILES or (not with_bam and f["name"].startswith("input.bam")):
            continue
        dest = os.path.join(d, f["name"])
        if not os.path.exists(dest):
            urllib.request.urlretrieve(f["url"], dest + ".part")
            os.replace(dest + ".part", dest)
    region = p.leaderboard(round_id)["round"]["region"]
    with open(os.path.join(d, "round.json"), "w") as fh:
        json.dump({"round_id": round_id, "region": region, "verification": ver}, fh)
    with open(os.path.join(d, "scores.json"), "w") as fh:
        json.dump(p.scores(round_id), fh)
    return d


def load(tag: str, reference: str = "") -> Round:
    d = round_dir(tag)
    with open(os.path.join(d, "round.json")) as fh:
        meta = json.load(fh)
    chrom, start, end = Round.parse_region(meta["region"])
    return Round(meta["round_id"], chrom, start, end, bam=os.path.join(d, "input.bam"),
                 reference=reference, released_dir=d)


def field_scores(tag: str) -> List[float]:
    """The real field's final scores for the round (0-100), one per submission, best first."""
    with open(os.path.join(round_dir(tag), "scores.json")) as fh:
        rows = json.load(fh)["scores"]
    by_sub = {}
    for r in rows:
        if r.get("combined_final") is not None:
            by_sub[r["submission_id"]] = 100.0 * r["combined_final"]
    return sorted(by_sub.values(), reverse=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--last", type=int, default=20)
    ap.add_argument("--no-bam", action="store_true")
    a = ap.parse_args()
    p = Platform()
    for r in p.rounds(limit=100):
        if a.last <= 0:
            break
        if r.get("is_finalized"):
            print(fetch(r["round_id"], p, with_bam=not a.no_bam))
            a.last -= 1


if __name__ == "__main__":
    main()
