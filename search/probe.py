"""A one-off fleet probe of the centre-jump policies (test only: MINOS_PROBE=<json>). No search runs.

For each centre the JSON names (a config key, or its first 8 characters, from this round's results store) and the
start of the run it came from, both policies read that run's scores from before the centre was scored, as a climber
standing there would have:

    v1   nearest first, interleaved with expected value x nearness (leaps.order)
    v2   the bandit's priors (leaps.pick), Thompson sampling of win rate x gain per win

Each sends its first n jumps (tags P1, P2; a config both chose is sent once as P12) and the probe waits for them.
The summary (sent, scored, beat the centre by climb.accept_min, best gain) is logged and written next to the replay
logs as probe_<tag>_<HHMM>.json.

    {"n": 40, "centres": [{"key": "13b3d031", "label": "17:06 peak", "from": "2026-10-15 17:07:01"}, ...]}
"""
from __future__ import annotations

import asyncio
import json
import logging
import os
import random
import sqlite3
import time
from typing import Dict, List

from common.tuning import tune
from common.types import FULL, Config
from search import leaps
from search.space import ACTIVITY_PARAMS, ASSEMBLY_PARAMS, FULL_RUN_PARAMS
from search.steps import Ladders

log = logging.getLogger(__name__)


class Past:
    """What a climber at a past centre could read: that run's scores before the centre, on this round's space."""

    def __init__(self, space, rows, seed: int):
        self.space = space
        self.ladders = Ladders(space)
        self.rng = random.Random(seed)
        self.done: Dict[str, tuple] = {c.key: (c, s) for c, s in rows}

    def clock(self) -> float:
        return 0.0

    def score(self, c: Config):
        hit = self.done.get(c.key)
        return None if hit is None else hit[1]

    def scored(self):
        return list(self.done.values())

    def fresh(self, c: Config) -> bool:
        return c is not None and c.key not in self.done


def _ts(text: str) -> float:
    return time.mktime(time.strptime(text, "%Y-%m-%d %H:%M:%S"))


async def run(client, space, path: str) -> None:
    spec = json.load(open(path))
    n = int(spec.get("n", 40))
    db = getattr(client.store, "path", None)
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True, timeout=60)
    names_all = [x for x in FULL_RUN_PARAMS + ASSEMBLY_PARAMS + ACTIVITY_PARAMS if x in space.names]
    plan: List[dict] = []
    send: Dict[str, Config] = {}
    tags: Dict[str, set] = {}
    for i, cen in enumerate(spec["centres"]):
        row = con.execute("select config_key, config, score, created from results where config_key like ? and score is not null "
                          "and created >= ? order by created limit 1", (cen["key"] + "%", _ts(cen["from"]))).fetchone()
        if row is None:
            log.info("engine: probe: centre %s not found after %s", cen["key"], cen["from"])
            continue
        key, raw, cscore, created = row
        centre = Config.from_dict(json.loads(raw))
        hist = []
        for k2, raw2, s2 in con.execute("select config_key, config, score from results where created >= ? and created < ? "
                                        "and score is not null and (gate_passed is null or gate_passed != 0)",
                                        (_ts(cen["from"]), created)):
            hist.append((Config.from_dict(json.loads(raw2)), float(s2)))
        past = Past(space, hist, seed=1000 + i)
        names = [x for x in names_all if past.ladders.movable(x)]
        v1 = leaps.order(past, centre, names, lambda c: False)[:n]
        v2 = [g[3] for g in leaps.pick(past, centre, names, lambda c: False, n, past.rng)]
        entry = {"label": cen.get("label", key[:8]), "key": key, "score": cscore, "history": len(hist),
                 "v1": [c.key for c in v1], "v2": [c.key for c in v2]}
        plan.append(entry)
        for policy, configs in (("P1", v1), ("P2", v2)):
            for c in configs:
                send[c.key] = c
                tags.setdefault(c.key, set()).add(policy)
        log.info("engine: probe: %s %s %.4f, %d earlier scores; v1 %d jumps, v2 %d, %d in both", entry["label"], key[:8],
                 cscore, len(hist), len(v1), len(v2), len(set(entry["v1"]) & set(entry["v2"])))
    con.close()
    futs = []
    for k, c in send.items():
        tag = "P12" if len(tags[k]) == 2 else next(iter(tags[k]))
        futs += client.propose_many([c], tag=tag, batch="small", priority=1)
    log.info("engine: probe: %d jumps sent", len(futs))
    pending = set(futs)
    while pending:
        t = client.time_left()
        if t is not None and t <= 5:
            break
        _, pending = await asyncio.wait(pending, timeout=2.0, return_when=asyncio.FIRST_COMPLETED)
    accept = float(tune("climb.accept_min"))
    out = []
    for entry in plan:
        row = {"label": entry["label"], "key": entry["key"], "score": entry["score"], "history": entry["history"]}
        for policy in ("v1", "v2"):
            scores = []
            for k in entry[policy]:
                r = client.latest.get((k, FULL))
                if r is not None and r.score is not None and r.gate_passed is not False:
                    scores.append(float(r.score))
            wins = [s - entry["score"] for s in scores if s - entry["score"] >= accept]
            row[policy] = {"sent": len(entry[policy]), "scored": len(scores), "wins": len(wins),
                           "best_gain": round(max((s - entry["score"] for s in scores), default=float("nan")), 4),
                           "best": round(max(scores, default=float("nan")), 4)}
        out.append(row)
        log.info("engine: probe: %s %.4f | v1 %d/%d beat it, best %+.4f | v2 %d/%d beat it, best %+.4f", row["label"],
                 row["score"], row["v1"]["wins"], row["v1"]["scored"], row["v1"]["best_gain"], row["v2"]["wins"],
                 row["v2"]["scored"], row["v2"]["best_gain"])
    runs = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test", "runs")
    os.makedirs(runs, exist_ok=True)
    dest = os.path.join(runs, f"probe_{client.round.tag}_{time.strftime('%H%M')}.json")
    with open(dest, "w") as f:
        json.dump(out, f, indent=1)
    log.info("engine: probe: summary in %s", dest)
