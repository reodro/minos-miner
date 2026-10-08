"""Seeds for the opening: the legacy portfolio plus what our own reviewed rounds add. Owner: search.

  legacy   search/portfolio/legacy.json: the legacy pipeline's portfolio (749 configs from 79 rounds, 2026-09-06 to
           09-18; per round up to 10 diverse configs within 0.3 of its best, spaced apart in normalised setting
           space - his portfolio_select.py). legacy_winners.json (624) is kept as a backup and not read.
  ours     $DATA_DIR/portfolio/ours.jsonl: the same rule applied to each of our reviewed runs (add_run), newest
           first. Not in git: it grows with every reviewed round.

A member has a profile: "A" for rounds with 250 bp reads (HG005), "B" for every other donor (his rule). pick()
takes the round's profile first, best-ranked and newest first, skips a member too close to one already taken,
and fills up with the other profile only when the round's own runs short. Members carry only configs: a new
round has a new donor and region, so their old scores mean nothing here. A reviewed run also keeps a high-core
group and a high false-positive group, tagged, so the next opening can draw both.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from common.settings import SETTINGS
from common.types import Config

HERE = Path(__file__).resolve().parent / "portfolio"
LEGACY = HERE / "legacy.json"
OURS = Path(SETTINGS.data_dir) / "portfolio" / "ours.jsonl"
LONG_READS = 240               # most reads at least this long: profile A
MIN_DISTANCE = 0.05            # members closer than this (mean normalised difference) are one member
SCORE_WINDOW, PER_RUN_MAX = 0.3, 10
GROUP_MAX = 5                 # high-core and high-fp members kept beside the near-best score set
CORE_BAND, FP_BAND = 0.01, 0.25   # how far below the best of that measure a member may sit
_TAG = re.compile(r"(2026\d{4}T\d{6,7})")


def bam_profile(bam: str, chrom: str, start: int, end: int, reads: int = 300) -> Optional[str]:
    """"A" when most reads of the region are at least LONG_READS long, else "B"; None when the BAM cannot be read."""
    try:
        import pysam
        lengths: Dict[int, int] = {}
        with pysam.AlignmentFile(bam) as f:
            for i, r in enumerate(f.fetch(chrom, max(0, start - 1), end)):
                lengths[r.query_length] = lengths.get(r.query_length, 0) + 1
                if i >= reads:
                    break
    except Exception:  # noqa: BLE001 - no BAM, no index, no pysam: no profile
        return None
    if not lengths:
        return None
    return "A" if max(lengths, key=lengths.get) >= LONG_READS else "B"


def _round_of(src: str) -> str:
    m = _TAG.search(src or "")
    return m.group(1)[:16] if m else ""


def members(legacy: Path = LEGACY, ours: Path = OURS) -> List[dict]:
    """Every member, ours first (newest round first), then the legacy ones (newest round first, best rank first)."""
    out: List[dict] = []
    try:
        lines = ours.read_text().splitlines()
    except OSError:
        lines = []
    mine = [json.loads(l) for l in lines if l.strip()]
    mine.sort(key=lambda m: (_round_of(m.get("src", "")), -int(m.get("rank", 1))), reverse=True)
    out += mine
    try:
        doc = json.loads(legacy.read_text())
    except (OSError, ValueError):
        doc = {}
    old = doc.get("members", doc) if isinstance(doc, dict) else doc
    old = sorted(old, key=lambda m: (int(m.get("rank") or 1), _neg(_round_of(m.get("src", "")))))
    return out + old


def _neg(tag: str) -> Tuple[int, ...]:
    return tuple(-ord(ch) for ch in tag)              # newest round first within one rank


def distance(space, a: Config, b: Config) -> float:
    """Mean normalised difference over the space's settings: |unit a - unit b| for numbers, 0/1 for choices."""
    total, n = 0.0, 0
    for name in space.names:
        s = space.specs[name]
        x, y = a.params.get(name, s.default), b.params.get(name, s.default)
        if s.allowed or s.type == "bool":
            total += 0.0 if x == y else 1.0
        else:
            total += abs(space._to_unit(s, x) - space._to_unit(s, y))
        n += 1
    return total / max(1, n)


def _member_group(m: dict) -> str:
    return m.get("group") if m.get("group") in ("score", "core", "fp") else "score"


def pick(space, profile: Optional[str], n: int, pool: Optional[Sequence[dict]] = None) -> List[Config]:
    """Up to n seeds: the round's profile first, then the other. Within a profile, score, core and fp members
    are taken in turn, each at least MIN_DISTANCE from those taken."""
    if n <= 0:
        return []
    pool = list(members() if pool is None else pool)
    own = [m for m in pool if profile is None or m.get("profile") in (profile, None)]
    rest = [m for m in pool if m not in own]
    out: List[Config] = []

    def take(group: List[dict]) -> None:
        buckets: Dict[str, List[dict]] = {"score": [], "core": [], "fp": []}
        for m in group:
            buckets[_member_group(m)].append(m)
        at = {k: 0 for k in buckets}
        moved = True
        while len(out) < n and moved:
            moved = False
            for name in ("score", "core", "fp"):
                while at[name] < len(buckets[name]):
                    m = buckets[name][at[name]]
                    at[name] += 1
                    c = space._assemble({k: v for k, v in (m.get("config") or {}).items() if k in space.names})
                    if c is None or any(c.key == x.key or distance(space, c, x) < MIN_DISTANCE for x in out):
                        continue
                    out.append(c)
                    moved = True
                    break
                if len(out) >= n:
                    return

    take(own)
    if len(out) < n:
        take(rest)
    return out


def select_run(space, configs: Dict[str, dict], tag: str, profile: Optional[str],
               known: Iterable[Config] = ()) -> List[dict]:
    """A run's contribution: configs within SCORE_WINDOW of its best local score, spaced apart (greedy
    farthest-point from the best), at most PER_RUN_MAX, none closer than MIN_DISTANCE to `known`."""
    rows = [(k, v) for k, v in configs.items() if v.get("local") is not None and v.get("config")]
    if not rows:
        return []
    top = max(v["local"] for _, v in rows)
    near = sorted(((k, v) for k, v in rows if v["local"] >= top - SCORE_WINDOW), key=lambda kv: -kv[1]["local"])
    cand = []
    for k, v in near:                                   # only configs this space accepts (pinned settings ours)
        c = space._assemble({n: x for n, x in Config.from_dict(v["config"]).params.items() if n in space.names})
        if c is not None:
            cand.append((c, v["local"]))
    taken: List[Tuple[Config, float]] = []
    others = list(known)
    while cand and len(taken) < PER_RUN_MAX:
        if not taken:
            c, s = cand.pop(0)
        else:                                           # the candidate farthest from everything taken
            i = max(range(len(cand)), key=lambda j: min(distance(space, cand[j][0], t) for t, _ in taken))
            c, s = cand.pop(i)
        if any(distance(space, c, t) < MIN_DISTANCE for t, _ in taken) or \
                any(distance(space, c, o) < MIN_DISTANCE for o in others):
            continue
        taken.append((c, s))
    rows = [{"src": f"minos:{tag}", "config": c.params, "profile": profile, "rank": i + 1, "local": round(s, 4),
             "group": "score", "at": round(time.time(), 1)} for i, (c, s) in enumerate(taken)]
    have = {c.key for c, _ in taken}

    def more(metric: str, group: str) -> None:
        cand = []
        for k, v in configs.items():
            raw = v.get(metric)
            if raw is None or v.get("local") is None or not v.get("config"):
                continue
            c = space._assemble({n: x for n, x in Config.from_dict(v["config"]).params.items() if n in space.names})
            if c is None or c.key in have:
                continue
            cand.append((c, float(raw), v["local"]))
        if not cand:
            return
        best = max(m for _, m, _ in cand)
        band = CORE_BAND if metric == "core" else FP_BAND
        cand = [t for t in cand if t[1] >= best - band]
        cand.sort(key=lambda t: -t[1])
        chosen: List[Tuple[Config, float]] = []
        while cand and len(chosen) < GROUP_MAX:
            if not chosen:
                c, _metric, local = cand.pop(0)
            else:
                i = max(range(len(cand)), key=lambda j: min(distance(space, cand[j][0], t) for t, _ in chosen))
                c, _metric, local = cand.pop(i)
            if any(distance(space, c, t) < MIN_DISTANCE for t, _ in taken) or \
                    any(distance(space, c, t) < MIN_DISTANCE for t, _ in chosen) or \
                    any(distance(space, c, o) < MIN_DISTANCE for o in others):
                continue
            chosen.append((c, local))
            have.add(c.key)
        base = len(rows)
        for i, (c, local) in enumerate(chosen):
            rows.append({"src": f"minos:{tag}", "config": c.params, "profile": profile, "rank": base + i + 1,
                         "local": round(local, 4), "group": group, "at": round(time.time(), 1)})

    more("core", "core")
    more("fp_per_target", "fp")
    return rows


def add_run(space, configs: Dict[str, dict], tag: str, profile: Optional[str], path: Path = OURS) -> int:
    """Append a reviewed run's members to ours.jsonl; returns how many."""
    known = []
    for m in members(ours=path):
        c = space._assemble({k: v for k, v in (m.get("config") or {}).items() if k in space.names})
        if c is not None:
            known.append(c)
    new = select_run(space, configs, tag, profile, known)
    if new:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "a") as f:
            for m in new:
                f.write(json.dumps(m) + "\n")
    return len(new)
