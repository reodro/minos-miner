"""A one-off focused search from given centres (test only: MINOS_POLISH=<json>). No islands, no Step A, no GA.

2026-10-17: round 20261016T130600Z ended 7th at 91.1563, 0.005 behind #1 (91.1613). Its winning step was a jump in
two genotyping settings, and around every recent pick the calling threshold had been tried at 1 of 74 values. This mode
climbs from the best configs of that run for the replay's minutes, to show how far the neighbourhood goes:

    centre    the best config so far (scored again in this run, so it can be the pick)
    cheap     with the centre's run: the genotyping settings over their whole ladder (every rung near the centre,
              every third rung beyond) and the PairHMM settings within 4 rungs; riders share the centre's GATK run
    heavy     every read, assembly and activity setting alone at 1, 2, 3, 5 and 8 rungs either way
    pairs     pairs and triples of the best single moves (a near-tie counts: small gains add up)
    jumps     centre jumps from each line's centre (search/leaps.py, the engine's arms and priors, over the scored history
              in spec "history_dbs"): single-setting jumps of 9+ rungs and two-setting jumps; they run before the
              random fill (2026-10-17: the round's best is a strict local optimum; the 7th-place step was a jump)
    fill      random two- and three-setting moves within 3 rungs, on workers nothing else needs
    move      any higher score becomes the centre at once; the newest centre's work is queued first
    keep      a line may freeze settings ("freeze": {setting: value}) and guard a class ("keep_class": "indel_8bp"):
              its moves never change a frozen setting, and a config moves the line only when that class scores at least
              as well as at the line's start (2026-10-17: unguarded, the "indel 8bp +2" line climbed back to the round's
              best by undoing the very change that won the two indels)
    lines     each spec centre climbs on its own line: a config planned from a line that beats that line's best moves
              the line (2026-10-17: the change that gains two 8 bp+ indels costs 0.10 elsewhere at the round's best, so
              it is climbed as its own line, where calling threshold and the indel settings can repair the cost)

A config the earlier run already scored (spec "known") is not sent again. Progress is logged as "engine: polish: ...";
the summary is written next to the replay logs as polish_<tag>_<HHMM>_result.json.

    {"target": 91.1613, "centres": [{"key": ..., "config": {...}, "score": ...}], "known": {"<key>": <score>}}
"""
from __future__ import annotations

import asyncio
import itertools
import json
import logging
import os
import random
import time
from typing import Dict, List, Optional, Tuple

from common.types import FULL, Config
from search.space import ACTIVITY_PARAMS, ASSEMBLY_PARAMS, FULL_RUN_PARAMS, GENOTYPING_PARAMS, PAIRHMM_PARAMS
from search.steps import Ladders

log = logging.getLogger(__name__)
HEAVY_STEPS = (1, -1, 2, -2, 3, -3, 5, -5, 8, -8)
STOP_S = 200.0            # nothing new is sent once a run could not come back before the end
TICK_S = 2.0


class Polish:
    def __init__(self, space, spec: dict, rng: Optional[random.Random] = None):
        self.space = space
        self.lad = Ladders(space)
        self.rng = rng or random.Random(7)
        self.known: Dict[str, float] = {k: float(v) for k, v in spec.get("known", {}).items()}
        self.target = spec.get("target")
        self.cheap = [n for n in GENOTYPING_PARAMS + PAIRHMM_PARAMS if n in space.names and self.lad.movable(n)]
        self.heavy = [n for n in FULL_RUN_PARAMS + ASSEMBLY_PARAMS + ACTIVITY_PARAMS
                      if n in space.names and self.lad.movable(n)]
        self.centres: List[Tuple[Config, float]] = [(Config.from_dict(c["config"]), float(c["score"]))
                                                    for c in spec["centres"]]
        self.lines: List[dict] = [{"best_c": c, "best_s": sc, "label": spec["centres"][i].get("label", c.key[:8]),
                                   "start": sc, "freeze": dict(spec["centres"][i].get("freeze", {})),
                                   "keep_class": spec["centres"][i].get("keep_class"),
                                   "keep_min": spec["centres"][i].get("keep_min_f1")}
                                  for i, (c, sc) in enumerate(self.centres)]
        self.line_of: Dict[str, int] = {}
        self.sent: set = set()
        self.queue: List[Tuple[int, List[Config], List[str]]] = []      # (priority, batch, labels)
        self.moves: Dict[str, List[Tuple[str, object, Config]]] = {}    # centre key -> its single moves
        self.combined: set = set()
        self.wave = 0
        self.history: List[dict] = []
        self.past = None                                                 # scored history the jump arms read
        self.jumped: Dict[str, int] = {}                                 # centre key -> jump passes planned

    # ---- configs -----------------------------------------------------------------------------------------------------
    def with_(self, base: Config, changes: Dict[str, object]) -> Optional[Config]:
        return self.space._assemble({**base.params, **changes})

    def new(self, c: Optional[Config]) -> bool:
        return c is not None and c.key not in self.sent and c.key not in self.known

    def cheap_values(self, centre: Config, n: str) -> List[object]:
        rungs = self.lad.rungs[n]
        i = int(round(self.lad.index(n, centre.params.get(n))))
        if n in GENOTYPING_PARAMS:
            idx = set(range(i - 6, i + 7)) | set(range(0, len(rungs), 3))
        else:
            idx = set(range(i - 4, i + 5))
        return [rungs[j] for j in sorted(idx, key=lambda j: (abs(j - i), j)) if 0 <= j < len(rungs) and j != i]

    # ---- planning ----------------------------------------------------------------------------------------------------
    def plan_centre(self, centre: Config, line: Optional[int] = None) -> int:
        """The centre with its cheap sweep (riders), then every heavy single. Returns how many configs were queued."""
        self.wave += 1
        prio = 100 + self.wave                    # the newest centre is taken first
        moves: List[Tuple[str, object, Config]] = []
        geno, pair = [], []
        for n in self.cheap:
            for v in self.cheap_values(centre, n):
                c = self.with_(centre, {n: v})
                if c is None or c.key == centre.key:
                    continue
                moves.append((n, v, c))
                if self.new(c):
                    (geno if n in GENOTYPING_PARAMS else pair).append(c)
        head = [centre] if centre.key not in self.sent else []
        n_q = 0
        for chunk in [geno[i:i + 60] for i in range(0, len(geno), 60)] or [[]]:
            batch = head + chunk
            head = []
            if batch:
                self.push(prio, batch, line)
                n_q += len(batch)
        for i in range(0, len(pair), 6):          # a PairHMM rider costs about 0.12 of a run
            self.push(prio, [centre] + pair[i:i + 6] if centre.key not in self.sent else pair[i:i + 6], line)
            n_q += len(pair[i:i + 6])
        for n in self.heavy:
            for k in HEAVY_STEPS:
                v = self.lad.shift(n, centre.params.get(n), k)
                c = None if v is None else self.with_(centre, {n: v})
                if c is None or c.key == centre.key:
                    continue
                moves.append((n, v, c))
                if self.new(c):
                    self.push(prio - 1, [c], line)
                    n_q += 1
        self.moves[centre.key] = moves
        return n_q

    def plan_extra(self, configs: List[Config], line: Optional[int] = 0) -> int:
        """Configs a teammate planned (spec "extra"), first in the queue. Those sharing the heavy and PairHMM settings
        go out as one run (a head and its genotyping riders, ~1% of a run each)."""
        cheap_geno = set(GENOTYPING_PARAMS)
        groups: Dict[tuple, List[Config]] = {}
        for c in configs:
            if self.new(c):
                key = tuple(sorted((k, repr(v)) for k, v in c.params.items() if k not in cheap_geno))
                groups.setdefault(key, []).append(c)
        n = 0
        for batch in groups.values():
            for i in range(0, len(batch), 40):
                self.push(200, batch[i:i + 40], line)
                n += len(batch[i:i + 40])
        return n

    def plan_pairs(self, centre: Config, score_of, line: Optional[int] = None) -> int:
        """Pairs and triples of the best single moves from `centre` that are back (near-ties included)."""
        scored = [(score_of(c), n, v) for n, v, c in self.moves.get(centre.key, ()) if score_of(c) is not None]
        base = score_of(centre)
        if base is None or len(scored) < 2:
            return 0
        best_per: Dict[str, Tuple[float, object]] = {}
        for s, n, v in scored:
            if s >= base - 0.003 and (n not in best_per or s > best_per[n][0]):
                best_per[n] = (s, v)
        top = sorted(best_per.items(), key=lambda kv: -kv[1][0])[:8]
        out = 0
        for r in (2, 3):
            for combo in itertools.combinations(top, r):
                if r == 3 and any(x[0] not in dict(top[:4]) for x in combo):
                    continue
                c = self.with_(centre, {n: sv[1] for n, sv in combo})
                if self.new(c) and c.key not in self.combined:
                    self.combined.add(c.key)
                    self.push(99 + self.wave, [c], line)
                    out += 1
        return out

    def plan_jumps(self, centre: Config, line: Optional[int], n: int) -> int:
        """Up to n centre jumps from `centre`: about half single-setting jumps of 9+ rungs, the rest two-setting jumps,
        both drawn from the engine's jump arms. Nothing frozen on the line moves; nothing known is sent again."""
        if self.past is None:
            return 0
        from search import leaps
        names = [x for x in self.heavy + self.cheap if not (line is not None and x in self.lines[line]["freeze"])]

        def skip(c: Config) -> bool:
            return not self.new(c) or not self.allowed(c, line)
        out = 0
        try:
            singles = [g for g in leaps.pick(self.past, centre, names, skip, 4 * n, self.rng) if g[1] >= 2][:max(1, n // 2)]
            multis = leaps.mutations(self.past, centre, names, skip, n - len(singles), self.rng, sizes=(2,))
        except Exception:  # noqa: BLE001 - jumps are an extra; a failure leaves the rest of the polish running
            log.exception("engine: polish: centre jumps failed")
            return 0
        for g in singles + multis:
            c = g[3]
            if not skip(c):
                self.push(50, [c], line)
                out += 1
        self.jumped[centre.key] = self.jumped.get(centre.key, 0) + 1
        return out

    def plan_fill(self, centre: Config, n: int, line: Optional[int] = None) -> int:
        names = self.heavy + self.cheap
        out = 0
        for _ in range(n * 4):
            if out >= n:
                break
            k = self.rng.choice((2, 2, 3))
            changes = {}
            for name in self.rng.sample(names, k):
                v = self.lad.shift(name, centre.params.get(name), self.rng.choice((1, -1, 2, -2, 3, -3)))
                if v is not None:
                    changes[name] = v
            c = self.with_(centre, changes) if changes else None
            if self.new(c):
                self.push(10, [c], line)
                out += 1
        return out

    def allowed(self, c: Config, line: Optional[int]) -> bool:
        """A config a frozen line may send: every frozen setting at its frozen value."""
        if line is None:
            return True
        fz = self.lines[line]["freeze"] if line < len(self.lines) else {}
        return all(c.params.get(k) == v for k, v in fz.items())

    def keeps(self, line: int, result) -> bool:
        """The line's guarded class scores at least its "keep_min_f1" in this result."""
        ln = self.lines[line]
        cls, need = ln.get("keep_class"), ln.get("keep_min")
        if not cls or need is None:
            return True
        cc = (getattr(result, "counts", None) or {}).get(cls)
        if cc is None:
            return False
        den = 2 * cc.tp + cc.fp + cc.fn
        return den > 0 and 2 * cc.tp / den >= need - 1e-12

    def push(self, prio: int, batch: List[Config], line: Optional[int] = None) -> None:
        batch = [c for c in batch if c is not None and self.allowed(c, line)]
        if not batch:
            return
        for c in batch:
            self.sent.add(c.key)
            if line is not None:
                self.line_of.setdefault(c.key, line)
        self.queue.append((prio, batch, ["small"] + ["big"] * (len(batch) - 1)))

    def take(self, slots: int) -> List[Tuple[int, List[Config], List[str]]]:
        self.queue.sort(key=lambda q: -q[0])
        out, self.queue = self.queue[:max(0, slots)], self.queue[max(0, slots):]
        return out


async def run(client, space, path: str) -> None:
    spec = json.load(open(path))
    p = Polish(space, spec)
    t0 = time.time()

    def score_of(c: Config) -> Optional[float]:
        r = client.latest.get((c.key, FULL))
        if r is not None and r.score is not None and r.gate_passed is not False:
            return float(r.score)
        return p.known.get(c.key)

    best_c, best_s = max(p.centres, key=lambda cs: cs[1])
    log.info("engine: polish: %d centres from the earlier run, best %s %.4f; target %s; %d configs already known",
             len(p.centres), best_c.key[:8], best_s, p.target, len(p.known))
    for ln_id, items in (spec.get("extra_lines") or {}).items():     # planned configs for one line (its guard applies)
        n = p.plan_extra([Config.from_dict(d) for d in items], int(ln_id))
        log.info("engine: polish: %d planned configs for line %s queued first", n, ln_id)
    extra = [Config.from_dict(d) for d in spec.get("extra", [])]
    if extra:
        n = p.plan_extra(extra, 0)
        log.info("engine: polish: %d planned configs from the spec queued first (%d runs)", n,
                 sum(1 for q in p.queue if q[0] == 200))
    dbs = [d for d in spec.get("history_dbs", []) if os.path.isfile(d)]
    if dbs:
        import sqlite3
        from search.probe import Past
        rows: Dict[str, Tuple[Config, float]] = {}
        for d in dbs:
            con = sqlite3.connect(f"file:{d}?mode=ro", uri=True, timeout=60)
            for k, raw, sc in con.execute("select config_key, config, score from results where score is not null "
                                          "and (gate_passed is null or gate_passed != 0)"):
                if k not in rows or sc > rows[k][1]:
                    rows[k] = (Config.from_dict(json.loads(raw)), float(sc))
            con.close()
        p.past = Past(space, list(rows.values()), seed=11)
        log.info("engine: polish: jump arms read %d scored configs", len(rows))
    for i, (c, s) in enumerate(p.centres):
        n = p.plan_centre(c, i)
        nj = p.plan_jumps(c, i, 60)
        log.info("engine: polish: line %d (%s) from %s %.4f: %d configs and %d centre jumps queued", i,
                 p.lines[i]["label"], c.key[:8], s, n, nj)
    paired: set = set()
    seen: set = set()
    rr = 0
    last_log = 0.0
    futs: List[asyncio.Future] = []
    while True:
        left = client.time_left()
        if left is not None and left <= STOP_S:
            break
        # move at once: a line moves on its own best; the round's best is the pick
        for (key, fid), r in list(client.latest.items()):
            if fid != FULL or key in seen or r.score is None or r.gate_passed is False:
                continue
            seen.add(key)
            c = client.configs.get(key)
            if c is None:
                continue
            sc = float(r.score)
            if p.past is not None:
                p.past.done[key] = (c, sc)
            line = p.line_of.get(key)
            if line is not None and sc > p.lines[line]["best_s"] + 1e-9 and p.keeps(line, r):
                ln = p.lines[line]
                gain_l = sc - ln["best_s"]
                ln["best_c"], ln["best_s"] = c, sc
                n = p.plan_centre(c, line) + p.plan_jumps(c, line, 60)
                log.info("engine: polish: line %d (%s) new best %s %.4f (%+.4f; %+.4f from its start) at t+%.1fm; "
                         "%d configs queued", line, ln["label"], key[:8], sc, gain_l, sc - ln["start"],
                         (time.time() - t0) / 60, n)
            if sc > best_s + 1e-9:
                gain = sc - best_s
                best_c, best_s = c, sc
                p.history.append({"t_min": round((time.time() - t0) / 60, 1), "key": key, "score": round(best_s, 4),
                                  "gain": round(gain, 4), "line": line})
                log.info("engine: polish: new best %s %.4f (%+.4f) at t+%.1fm%s", key[:8], best_s, gain,
                         (time.time() - t0) / 60, " - ABOVE THE TARGET" if p.target and best_s > p.target else "")
        if best_c.key in {k for (k, f) in client.latest if f == FULL}:
            client.pick_key = best_c.key
        for i, ln in enumerate(p.lines):
            c = ln["best_c"]
            if c.key not in paired and len([m for m in p.moves.get(c.key, ()) if score_of(m[2]) is not None]) >= 20:
                paired.add(c.key)
                n = p.plan_pairs(c, score_of, i)
                if n:
                    log.info("engine: polish: %d pairs and triples of the best moves around %s", n, c.key[:8])
        slots = client.idle_slots()
        slots = 64 if slots is None else int(slots)
        if slots > 0 and not p.queue:
            rr = (rr + 1) % len(p.lines)                  # jumps, then random moves, take turns among the lines
            ln = p.lines[rr]
            if p.jumped.get(ln["best_c"].key, 0) < 6 and p.plan_jumps(ln["best_c"], rr, max(20, slots)):
                pass
            else:
                p.plan_fill(ln["best_c"], slots, rr)
        for prio, batch, labels in p.take(slots):
            futs += client.propose_many(batch, tag="P", batches=labels, priority=prio)
        if time.time() - last_log >= 60:
            last_log = time.time()
            done = sum(1 for f in futs if f.done())
            log.info("engine: polish: t+%.0fm best %.4f (target %s, gap %+.4f); %d sent, %d back, %d queued",
                     (time.time() - t0) / 60, best_s, p.target, (best_s - p.target) if p.target else float("nan"),
                     len(futs), done, sum(len(q[1]) for q in p.queue))
        pending = [f for f in futs if not f.done()]
        if pending:
            await asyncio.wait(pending, timeout=TICK_S, return_when=asyncio.FIRST_COMPLETED)
        else:
            await asyncio.sleep(TICK_S)
    out = {"start": round(max(s for _, s in p.centres), 4), "best": round(best_s, 4), "best_key": best_c.key,
           "target": p.target, "beat_target": bool(p.target and best_s > p.target), "sent": len(futs),
           "history": p.history,
           "lines": [{"label": ln["label"], "start": round(ln["start"], 4), "best": round(ln["best_s"], 4),
                      "best_key": ln["best_c"].key} for ln in p.lines]}
    log.info("engine: polish: done: %.4f -> %.4f (target %s)%s; %d configs sent", out["start"], best_s, p.target,
             " - the target is beaten" if out["beat_target"] else "", len(futs))
    runs = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test", "runs")
    os.makedirs(runs, exist_ok=True)
    with open(os.path.join(runs, f"polish_{client.round.tag}_{time.strftime('%H%M')}_result.json"), "w") as f:
        json.dump(out, f, indent=1)
