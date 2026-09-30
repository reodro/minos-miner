"""Released rounds on disk: data/rounds/<tag>/{input.bam,truth.vcf.gz,mutations.vcf.gz,...,round.json,scores.json}.

    python3 -m backtest.data --last 20        # fetch the 20 newest finalized rounds
"""
from __future__ import annotations

import argparse
import gzip
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
                 reference=reference, released_dir=d)   # n_targets stays None: unknown live, like here


def field_scores(tag: str) -> List[float]:
    """The real field's final scores for the round (0-100), one per submission, best first."""
    with open(os.path.join(round_dir(tag), "scores.json")) as fh:
        rows = json.load(fh)["scores"]
    by_sub = {}
    for r in rows:
        if r.get("combined_final") is not None:
            by_sub[r["submission_id"]] = 100.0 * r["combined_final"]
    return sorted(by_sub.values(), reverse=True)


def discover(since: str, platform: Optional[Platform] = None, max_missing: int = 20) -> List[str]:
    """Finalized round ids back to `since` (ISO time), newest first.

    The API lists only the latest 100 rounds, so older ids are found by stepping back
    72 minutes and probing +-5 minutes around each slot (the grid drifts by a minute
    across restarts). Stops after `max_missing` consecutive empty slots.
    """
    import datetime as dt
    import urllib.error
    p = platform or Platform()
    recent = [r for r in p.rounds(limit=100) if r.get("is_finalized")]
    ids = [r["round_id"] for r in recent]
    t = dt.datetime.fromisoformat(recent[-1]["round_id"])
    stop = dt.datetime.fromisoformat(since)
    missing = 0
    while t > stop and missing < max_missing:
        expected = t - dt.timedelta(minutes=72)
        found = None
        for off in (0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5):
            rid = (expected + dt.timedelta(minutes=off)).isoformat()
            try:
                if p.leaderboard(rid)["round"].get("is_finalized"):
                    found = rid
                    break
            except urllib.error.HTTPError:
                continue
        if found:
            ids.append(found)
            t = dt.datetime.fromisoformat(found)
            missing = 0
        else:
            t = expected
            missing += 1
    return [i for i in ids if dt.datetime.fromisoformat(i) >= stop]


def main():
    from concurrent.futures import ThreadPoolExecutor
    ap = argparse.ArgumentParser()
    ap.add_argument("--last", type=int, default=20, help="newest N finalized rounds")
    ap.add_argument("--since", help="every finalized round since this ISO time (overrides --last)")
    ap.add_argument("--no-bam", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    p = Platform()
    if a.since:
        ids = discover(a.since, p)
    else:
        ids = [r["round_id"] for r in p.rounds(limit=100) if r.get("is_finalized")][:a.last]
    print(f"{len(ids)} rounds", flush=True)

    def one(rid):
        try:
            return fetch(rid, p, with_bam=not a.no_bam)
        except Exception as e:  # one bad round must not stop the rest
            return f"FAILED {rid}: {e}"

    with ThreadPoolExecutor(a.workers) as pool:
        for line in pool.map(one, ids):
            print(line, flush=True)


if __name__ == "__main__":
    main()
