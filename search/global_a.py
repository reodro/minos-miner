"""Island mode (2026-10-20, islands.global_a): one global Step A for the whole fleet, then one island per segment.

Step A searches the influential settings (a broad range and a big score effect over 30 rounds, plus the narrow
high-impact pcr_indel_model and min_pruning) over their whole published range, then hands Step B the best distinct
segments. Each island climbs one segment with its influential settings held in a band around the segment's values,
so two islands cannot end on one hill (the two islands of v9.18 shared pcr_indel_model in 5 of 5 live rounds).

    A1  from 8 starts (the portfolio's best for each pcr_indel_model x min_pruning 2/3), every value of each
        influential setting, one setting at a time                                                       one wave
    A2  each setting's local peaks per start; candidates = the best PEAKS over all starts (pcr: all values)
    A3  every combination of the candidates, the other settings from the matching start                 one wave
    A4  the best TOP combinations that differ in an influential setting: +-3 then +-1 rungs on RETUNE  two waves
    band  island i climbs segment i (and i + k, i + 2k ... when a hill is found) with each influential setting within
        BAND rungs of the segment's value and pcr_indel_model fixed

Step A stops at islands.global_a_share of the round, whatever stage it is in: the waves still running are dropped and
the segments are built from what has scored. The same algorithm runs offline in test/bench/island_stepa.py.
"""
from __future__ import annotations

import asyncio
import copy
import dataclasses
import itertools
import json
import logging
import math
from typing import Dict, List, Optional, Sequence, Tuple

from common.tuning import tune
from common.types import Config
from search import portfolio
from search.steps import Ladders

log = logging.getLogger(__name__)

BROAD = ("base_quality_score_threshold", "active_probability_threshold", "min_base_quality_score",
         "phred_scaled_global_read_mismapping_rate")
PCR, MP = "pcr_indel_model", "min_pruning"
INFLUENTIAL = BROAD + (PCR, MP)
EXTRA = ("standard_min_confidence_threshold_for_calling", "min_mapping_quality_score")   # islands.global_a_extra
RETUNE = ("pair_hmm_gap_continuation_penalty", "recover_all_dangling_branches", "assembly_region_padding",
          "heterozygosity", "min_mapping_quality_score", "standard_min_confidence_threshold_for_calling")
PEAKS, SEP, TOP, BAND = 3, 2, 20, 2
# Two peaks of a setting are one hill unless the curve between them drops by more than VALLEY; a hill's band is every
# value within VALLEY of its top. 13:06 (2026-10-20): active_probability_threshold 0.001-0.023, base_quality_score_
# threshold 6-39 and the mismapping rate 26-60 were plateaus with bumps of at most 0.34, then cliffs of 2-86 points.
VALLEY = 0.5
HILL_DROP = 1.0          # a hill whose top is more than this below the best hill cannot be a winning island
PASSES = (3, 1, 1)       # A4 re-tune passes (rungs): few segments after A3, so the budget goes into the re-tune
PASSES_SHORT = ((3, 1), (1,))   # islands.global_a_short: +-3 and +-1 in one wave, then +-1 (two waves, not three)


def enabled() -> bool:
    return int(tune("islands.global_a")) > 0


RUNS_PER_START = 40      # GATK runs one start's A1 sweep needs after batching (live 07:22, 2026-10-21: ~330 for 8 starts)
MAX_STARTS = 32
TAIL_BACK = 0.6          # islands.global_a_tail: the share of a wave that must be back before an emptying fleet ends it


def n_starts(cap: int) -> int:
    """A1 is one wave: as many starts as the fleet runs at once (8 starts kept ~330 of 971 workers busy)."""
    return max(8, min(MAX_STARTS, int(cap or 0) // RUNS_PER_START))


def generalists() -> List[Dict[str, object]]:
    """Cross-round configs for A1 (test/bench/generalist.py, 2026-10-21: best mean TRUE over 15 released rounds; on
    held-out rounds they seldom crashed while half of the portfolio starts did on 2 of 5 rounds)."""
    from common.settings import SETTINGS
    try:
        return list(json.load(open(f"{SETTINGS.data_dir}/portfolio/generalists_a1.json")).get("configs") or [])
    except (OSError, ValueError):
        return []


def starts(space, profile: Optional[str], tag: str, n: int = 8) -> Dict[str, Config]:
    """The portfolio's best member for each pcr_indel_model x min_pruning 2/3 (else its top member with those set),
    then the cross-round generalists, then the next best distinct portfolio members, up to n starts."""
    ranked = portfolio.pick(space, profile, 300, exclude_round=tag)
    out: Dict[str, Config] = {}
    if not ranked:
        return out
    for pcr in space.specs[PCR].allowed:
        for mp in (2, 3):
            c = next((c for c in ranked if c.params.get(PCR) == pcr and c.params.get(MP) == mp), None)
            c = c or space._assemble(dict(ranked[0].params, **{PCR: pcr, MP: mp}))
            if c is not None:
                out[f"{pcr}/{mp}"] = c
    keys = {c.key for c in out.values()}
    for i, params in enumerate(generalists()):            # cross-round configs: they seldom crash (held-out test)
        if len(out) >= n:
            break
        c = space._assemble(params)
        if c is not None and c.key not in keys:
            keys.add(c.key)
            out[f"g{i}"] = c
    for i, c in enumerate(ranked):
        if len(out) >= n:
            break
        if c.key not in keys:
            keys.add(c.key)
            out[f"p{i}"] = c
    return out


def local_peaks(ys: Sequence[Optional[float]]) -> List[int]:
    idx = [i for i, y in enumerate(ys) if y is not None]
    cand = []
    for k, i in enumerate(idx):
        left = ys[idx[k - 1]] if k > 0 else float("-inf")
        right = ys[idx[k + 1]] if k + 1 < len(idx) else float("-inf")
        if ys[i] >= left and ys[i] >= right:
            cand.append(i)
    out: List[int] = []
    for i in sorted(cand, key=lambda i: -ys[i]):
        if all(abs(i - j) >= SEP for j in out):
            out.append(i)
        if len(out) == PEAKS:
            break
    return out


def hills(ys: Sequence[Optional[float]]) -> List[Tuple[int, int, int]]:
    """The hills of a curve, highest first, as (top, lo, hi) indices: local maxima merged unless the curve between them
    drops by more than VALLEY; a hill spans the contiguous values within VALLEY of its top."""
    idx = [i for i, y in enumerate(ys) if y is not None]
    if not idx:
        return []
    peaks = []
    for k, i in enumerate(idx):
        left = ys[idx[k - 1]] if k > 0 else float("-inf")
        right = ys[idx[k + 1]] if k + 1 < len(idx) else float("-inf")
        if ys[i] >= left and ys[i] >= right:
            peaks.append(i)
    tops: List[int] = []
    for p in sorted(peaks, key=lambda i: -ys[i]):
        joined = False
        for t in tops:
            a, b = min(p, t), max(p, t)
            low = min(ys[i] for i in idx if a <= i <= b)
            if low >= ys[p] - VALLEY:
                joined = True
                break
        if not joined:
            tops.append(p)
    out = []
    tops = [t for t in tops if ys[t] >= ys[tops[0]] - HILL_DROP]
    for t in tops:
        k = idx.index(t)
        lo = hi = k
        while lo > 0 and ys[idx[lo - 1]] >= ys[t] - VALLEY:
            lo -= 1
        while hi + 1 < len(idx) and ys[idx[hi + 1]] >= ys[t] - VALLEY:
            hi += 1
        out.append((t, idx[lo], idx[hi]))
    return out


def banded(space, centre: Dict[str, object], bands: Optional[Dict[str, Tuple[object, object]]] = None):
    """A copy of `space` with pcr_indel_model fixed and each other influential setting inside its hill around `centre`
    (`bands`: setting -> (lo, hi) values), else within BAND rungs of `centre`."""
    lad = Ladders(space)
    new = copy.copy(space)
    new.specs = dict(space.specs)
    for n in INFLUENTIAL:
        s = new.specs.get(n)
        if s is None or n not in lad.rungs or centre.get(n) is None:
            continue
        rungs = lad.rungs[n]
        if n == PCR:
            new.specs[n] = dataclasses.replace(s, allowed=(centre[n],))
            continue
        if bands and n in bands:
            # the hill, and at least BAND rungs either side of the centre: hills drawn from every scored config can be
            # one value wide around a best config that stands out (live 01:11: min_pruning 2..2, apt 0.023..0.023)
            lo, hi = (lad.index(n, v) for v in bands[n])
            i = lad.index(n, centre[n])
            lo, hi = min(lo, i - BAND), max(hi, i + BAND)
            band = [v for v in rungs if lo - 1e-9 <= lad.index(n, v) <= hi + 1e-9] or [centre[n]]
        else:
            i = int(round(lad.index(n, centre[n])))
            band = rungs[max(0, i - BAND):i + BAND + 1]
        if s.allowed:
            new.specs[n] = dataclasses.replace(s, allowed=tuple(v for v in s.allowed if v in band) or (centre[n],))
        else:
            new.specs[n] = dataclasses.replace(s, min=min(band), max=max(band))
    new.api = {n: {"min": s.min, "max": s.max, "allowed": s.allowed, "type": s.type} for n, s in new.specs.items()}
    new.names = [n for n, s in sorted(new.specs.items())
                 if not new._fixed(s) and n not in new.pinned and n in space.names]
    return new


class GlobalA:
    def __init__(self, client, space, profile: Optional[str], length_s: float, cap: int = 0):
        self.client, self.space, self.profile = client, space, profile
        self.cap = cap if int(tune("islands.global_a_fill")) > 0 else 0   # fill the waves after A1 up to this
        # the fleet in workers: how many A1 starts fit (with evaluator.split_k a node offers workers // k slots, each
        # finishing a run ~k times sooner, so the same starts fit in the same time)
        self.cap_all = cap * max(1, int(tune("evaluator.split_k")))
        self.lad = Ladders(space)
        self.t0 = client.clock()
        self.stop_at = self.t0 + length_s * float(tune("islands.global_a_share"))
        self.scores: Dict[str, Tuple[Config, float]] = {}
        self.spans: Dict[str, List[Tuple[object, object, object]]] = {}     # setting -> hills (top, lo, hi values)
        self.report: Dict[str, object] = {"stages": {}}

    def now(self) -> float:
        return self.client.clock()

    def assemble(self, params) -> Optional[Config]:
        c = self.space._assemble(params)
        return None if c is None or self.client.ranges.errors(c) else c

    async def wave(self, name: str, configs: List[Config]) -> Dict[str, float]:
        """Run the distinct configs. The wave ends once islands.global_a_quorum of them are back (scored or failed;
        live 07:22, 2026-10-21: A1 waited on its last 53 configs with 6 GATK runs going) or at the stop time; the rest
        keep running and still reach the store."""
        uniq = list({c.key: c for c in configs}.values())
        futs = self.client.propose_many(uniq, tag="A")
        left = self.stop_at - self.now()
        pending = set(futs)
        need = len(futs) - math.ceil(float(tune("islands.global_a_quorum")) * len(futs))   # may still be pending
        while len(pending) > need and left > 0 and not self.emptying(len(futs) - len(pending), len(futs)):
            _, pending = await asyncio.wait(pending, timeout=min(left, 5.0), return_when=asyncio.FIRST_COMPLETED)
            left = self.stop_at - self.now()
        out: Dict[str, float] = {}
        for c, f in zip(uniq, futs):
            if f.done() and not f.cancelled() and f.exception() is None:
                o = f.result()
                if o is not None and o.ok:
                    out[c.key] = o.score
                    self.scores[c.key] = (c, o.score)
        self.report["stages"][name] = {"sent": len(uniq), "scored": len(out), "cut": bool(pending),
                                       "minute": round((self.now() - self.t0) / 60, 1),
                                       "best": round(max(out.values()), 4) if out else None}
        self.publish()
        log.info("global A %s: %d sent, %d scored%s", name, len(uniq), len(out), ", cut at the stop time" if pending
                 else "")
        return out

    def emptying(self, back: int, sent: int) -> bool:
        """islands.global_a_tail (idle fix 2, 2026-10-22): the wave ends once TAIL_BACK of it is back and at least this
        share of the workers is idle (round 04:58: busy 993 -> 735 -> 317 -> 114 from t+8 to t+14 min while A1 waited for
        its quorum); its last runs keep going and still reach the store. 0 = off."""
        share = float(tune("islands.global_a_tail"))
        if share <= 0 or sent <= 0 or back < TAIL_BACK * sent:
            return False
        idle, cap = self.client.idle_slots(), self.client.capacity()
        return bool(idle is not None and cap and idle >= share * cap)

    async def wave_inline(self, name: str, st: Dict[str, Config], pool: Dict[str, Config], sweep_of) -> Dict[str, float]:
        """islands.global_a_inline (idle fix 1, 2026-10-22): the starts go out first (priority) TOGETHER with their sweeps,
        instead of a probe that kept 977 of 993 workers idle for 2.8 min (round 04:58; 3 crashed starts in 17 probes).
        A start whose own run fails is replaced by the next spare, which goes out with its sweep at once; the failed
        start's sweep is not used. The wave ends like wave(): quorum, the emptying fleet, or the stop time."""
        spares = [k for k in pool if k not in st]
        live: Dict[str, Config] = dict(st)
        futs: Dict[str, asyncio.Future] = {}
        cfgs: Dict[str, Config] = {}

        def send(sname: str, s: Config) -> None:
            first = [s] if s.key not in futs else []
            for c, f in zip(first, self.client.propose_many(first, tag="A", priority=10)):
                futs[c.key], cfgs[c.key] = f, c
            sweep = list({c.key: c for c in sweep_of(sname, s) if c.key not in futs}.values())
            for c, f in zip(sweep, self.client.propose_many(sweep, tag="A")):
                futs[c.key], cfgs[c.key] = f, c
        for k, s in st.items():
            send(k, s)
        crashed: List[str] = []
        checked: set = set()
        quorum = float(tune("islands.global_a_quorum"))
        left = self.stop_at - self.now()
        while left > 0:
            for k, s in list(live.items()):
                f = futs.get(s.key)
                if f is None or not f.done() or k in checked:
                    continue
                checked.add(k)
                o = f.result() if not f.cancelled() and f.exception() is None else None
                if o is None or not o.ok:
                    crashed.append(k)
                    del live[k]
                    if spares:
                        r = spares.pop(0)
                        live[r] = pool[r]
                        send(r, pool[r])
            pending = [f for f in futs.values() if not f.done()]
            back = len(futs) - len(pending)
            if not pending or back >= quorum * len(futs) or self.emptying(back, len(futs)):
                break
            await asyncio.wait(pending, timeout=min(left, 5.0), return_when=asyncio.FIRST_COMPLETED)
            left = self.stop_at - self.now()
        out: Dict[str, float] = {}
        for key, f in futs.items():
            if f.done() and not f.cancelled() and f.exception() is None:
                o = f.result()
                if o is not None and o.ok:
                    out[key] = o.score
                    self.scores[key] = (cfgs[key], o.score)
        cut = any(not f.done() for f in futs.values())
        self.report["probe"] = {"starts": len(st), "crashed": crashed, "replaced": len(crashed), "inline": True}
        self.report["stages"][name] = {"sent": len(futs), "scored": len(out), "cut": cut,
                                       "minute": round((self.now() - self.t0) / 60, 1),
                                       "best": round(max(out.values()), 4) if out else None}
        self.publish()
        log.info("global A %s (starts in flight): %d sent, %d scored, %d start(s) crashed and replaced (%s)%s", name,
                 len(futs), len(out), len(crashed), ", ".join(crashed) or "-", ", cut" if cut else "")
        return out

    def out_of_time(self) -> bool:
        return self.now() >= self.stop_at

    async def run(self, tag: str) -> List[Tuple[Config, float]]:
        n = n_starts(self.cap_all)
        pool = starts(self.space, self.profile, tag, n + 16)          # 16 spares replace crashing starts
        st = dict(list(pool.items())[:n])
        if not st:
            return []
        if int(tune("islands.global_a_scan")):
            # Build B (test, 2026-10-22): the generalized scan (search/global_g.py: 36 sub-domains over every searched
            # setting, riders, the starts as anchors) fills the scores instead of the one-setting sweep; hills,
            # segments, heads, bands and seats below stay as they are
            from search.global_g import GeneralizedA
            g = GeneralizedA(self.client, self.space, max(60.0, self.stop_at - self.now()),
                             runs=self.cap_all * max(1, int(tune("islands.global_g_waves"))), anchors=list(st.values()))
            await g.run()
            self.scores.update(g.scores)
            self.report["stages"]["A1"] = dict(g.report, scan="generalized")
            self.publish()
            return await self.after_a1(st)
        probe_s = float(tune("islands.global_a_probe"))
        inline = int(tune("islands.global_a_inline")) > 0
        if probe_s > 0 and not inline:
            st = await self.probe(st, pool, probe_s)
        # A1
        extra = int(tune("islands.global_a_extra")) > 0
        main8 = set(list(st)[:8])

        def sweep_of(sname: str, s: Config) -> List[Config]:
            names = INFLUENTIAL + (EXTRA if extra and sname in main8 else ())   # the friend's two settings, 8 main starts
            out = []
            for n in names:
                for v in self.lad.rungs.get(n, ()):
                    c = self.assemble(dict(s.params, **{n: v}))
                    if c is not None and c.params.get(n) == v:
                        out.append(c)
            return out
        if inline:
            got = await self.wave_inline("A1", st, pool, sweep_of)
        else:
            await self.wave("A1", [c for sname, s in st.items() for c in sweep_of(sname, s)])
        return await self.after_a1(st)

    async def after_a1(self, st: Dict[str, Config]) -> List[Tuple[Config, float]]:
        """Hills, then segments (one-wave Step A) or the A2/A3 waves, from every config scored so far."""
        if self.out_of_time():
            return self.segments(list(self.scores.values()))
        # A2
        cand: Dict[str, list] = {}
        for n in INFLUENTIAL:
            rungs = self.lad.rungs.get(n, [])
            if n == PCR:
                cand[n] = list(rungs)
                continue
            # the best score at each value over EVERY config scored so far, not only this setting's sweep: live
            # 01:11 (2026-10-21) the min_pruning sweep peaked at 4 while the round's best config (86.03) had 2, so the
            # best config was off the hills and headed no island
            env: Dict[str, float] = {}
            for c, sc in self.scores.values():
                key = json.dumps(c.params.get(n))
                env[key] = max(env.get(key, float("-inf")), sc)
            ys = [env.get(json.dumps(v)) for v in rungs]
            hs = hills(ys)[:PEAKS]
            self.spans[n] = [(rungs[t], rungs[lo], rungs[hi]) for t, lo, hi in hs]
            cand[n] = [rungs[t] for t, _, _ in hs] or [st[next(iter(st))].params.get(n)]
        self.report["candidates"] = cand
        self.report["hills"] = {n: [list(x) for x in v] for n, v in self.spans.items()}
        self.publish()
        log.info("global A candidates: %s", json.dumps(cand, default=str))
        if int(tune("islands.global_a_waves")) == 1:
            # one-wave Step A (user, 2026-10-21): islands straight from A1 (each headed by its best A1 config on its
            # hills); no A3 / re-tune waves: Step B's climbers do that on the whole fleet (live: those small waves left
            # 19-99 of 486 workers busy for ~10 minutes)
            if int(tune("islands.global_a_heads")):
                # Build B2 (test): the best complete configs head the islands, kept islands.seat_apart settings apart,
                # so a winning combination of values is not broken up into per-value hills
                return self.segments(self.top_configs())
            on_hills = [(c, s) for c, s in self.scores.values() if self.on_hills(c.params)]
            return self.segments(self.distinct(sorted(on_hills, key=lambda t: -t[1])))
        # A3
        grid = []
        for combo in itertools.product(*(cand[n] for n in INFLUENTIAL)):
            p = dict(zip(INFLUENTIAL, combo))
            base = st.get(f"{p[PCR]}/{2 if p[MP] <= 2 else 3}") or next(iter(st.values()))
            c = self.assemble(dict(base.params, **p))
            if c is not None:
                grid.append(c)
        short = int(tune("islands.global_a_short")) > 0
        if short:                           # only the islands no A1 point stands on get their A3 combination
            have = {self.island_key(c.params) for c, _ in self.scores.values() if self.on_hills(c.params)}
            grid = [c for c in grid if self.island_key(c.params) not in have]
        if grid or self.cap:
            # the free workers take neighbours of the best island heads so far (live 2026-10-21: after A1, 19 of 486
            # workers were busy for ~10 minutes)
            # centres: the 20 best scored configs (hills drawn from every config are tight around the best one, so
            # few configs stand on all of them yet: one centre gave 63 neighbours for 300 workers in the simulation)
            top = sorted(self.scores.values(), key=lambda t: -t[1])[:20]
            fill = self.fill([c for c, _ in top], [c for c, _ in top], len(grid), {c.key for c in grid})
            if grid or fill:
                await self.wave("A3", grid + [m for _, m in fill])
        # each island (pcr_indel_model + one hill per setting) is headed by its best scored config from A1 and A3,
        # so an island survives a crashed A3 combination; a config off every hill of a setting heads nothing
        on_hills = [(c, s) for c, s in self.scores.values()
                    if all(self.hill_of(n, c.params.get(n)) is not None for n in self.spans)]
        picked = self.distinct(sorted(on_hills, key=lambda t: -t[1]))
        # A4
        seg = [[c, s] for c, s in picked]
        heads = [c for c, _ in seg]
        for p, sizes in enumerate(PASSES_SHORT if short else tuple((x,) for x in PASSES)):
            if self.out_of_time():
                break
            moves: List[Tuple[int, Config]] = []
            for i, (c, _) in enumerate(seg):
                for n in RETUNE + BROAD + (MP,):            # influential settings move only inside their hill
                    rungs = self.lad.rungs.get(n, [])
                    if len(rungs) < 2:
                        continue
                    for k in [d for size in sizes for d in (size, -size)]:
                        j = int(round(self.lad.index(n, c.params.get(n)))) + k
                        if 0 <= j < len(rungs):
                            m = self.assemble(dict(c.params, **{n: rungs[j]}))
                            if m is not None and m.params.get(n) == rungs[j] and self.in_band(heads[i], m):
                                moves.append((i, m))
            fill = self.fill([c for c, _ in seg], heads, len(moves), {m.key for _, m in moves})
            moves += fill
            got = await self.wave(f"A4_{p + 1}_pm{'+'.join(map(str, sizes))}", [m for _, m in moves])
            for i, m in moves:
                if m.key in got and got[m.key] > seg[i][1]:
                    seg[i] = [m, got[m.key]]
        return self.segments([(c, s) for c, s in seg])

    async def probe(self, st: Dict[str, Config], pool: Dict[str, Config], probe_s: float) -> Dict[str, Config]:
        """islands.global_a_probe (2026-10-21): the starts run alone first, while the answer set builds; a start whose
        GATK run fails (the assembly crash "Graph must have ref source and sink vertices" kills every config of its
        sweep: 35% of A1 in round 04:58) is replaced by the next spare before its sweep is sent. A crash comes back as
        soon as GATK dies; a start still running after probe_s counts as fine."""
        names = list(st)
        futs = self.client.propose_many([st[k] for k in names], tag="A")
        t_end = self.now() + probe_s
        pending = set(futs)
        while pending and self.now() < min(t_end, self.stop_at):
            _, pending = await asyncio.wait(pending, timeout=5.0, return_when=asyncio.FIRST_COMPLETED)
        spares = [k for k in pool if k not in st]
        out: Dict[str, Config] = {}
        crashed = []
        for k, f in zip(names, futs):
            o = f.result() if f.done() and not f.cancelled() and f.exception() is None else None
            if f.done() and (o is None or not o.ok):
                crashed.append(k)
                if spares:
                    r = spares.pop(0)
                    out[r] = pool[r]
            else:
                out[k] = st[k]
                if o is not None and o.ok:
                    self.scores[st[k].key] = (st[k], o.score)
        self.report["probe"] = {"starts": len(names), "crashed": crashed, "replaced": len(crashed)}
        log.info("global A probe: %d starts, %d crashed and replaced (%s)", len(names), len(crashed),
                 ", ".join(crashed) or "-")
        return out

    def publish(self) -> None:
        """The report for the dashboard (main/status.py: search.global_a), refreshed after every wave."""
        try:
            self.report["stop_minute"] = round((self.stop_at - self.t0) / 60, 1)
            self.report["fill"] = bool(self.cap)
            self.client.notes["global_a"] = json.loads(json.dumps(self.report, default=str))
        except Exception:  # noqa: BLE001 - the dashboard must never break the search
            log.exception("global A: publishing the report failed")

    def fill(self, centres: List[Config], heads: List[Config], used: int, taken) -> List[Tuple[int, Config]]:
        """Extra configs for the free workers of a wave, counted in GATK runs: each centre's neighbours on every
        searched setting at +-1, +-2, +-3 rungs, inside the centre's island; never a config already scored or in the
        wave. Neighbours on settings that need their own GATK run (search/space.py FULL_RUN, activity, assembly) are
        taken first, one run each, centres in turn, until cap - used runs; the batchable ones (genotyping, PairHMM:
        they ride inside a run) of the same centres follow on top (live 02:26, 2026-10-21: 449 configs, mostly
        batchable, became ~99 runs on 486 workers)."""
        from search.space import GENOTYPING_PARAMS, PAIRHMM_PARAMS
        room = self.cap - used
        if room <= 0 or not centres:
            return []
        cheap_names = set(GENOTYPING_PARAMS) | set(PAIRHMM_PARAMS)
        seen = set(taken) | set(self.scores)
        costly: List[List[Tuple[int, Config]]] = []
        cheap: List[List[Tuple[int, Config]]] = []
        for i, c in enumerate(centres):
            mine_costly: List[Tuple[int, Config]] = []
            mine_cheap: List[Tuple[int, Config]] = []
            for size in (1, 2, 3):
                for n in self.space.names:
                    rungs = self.lad.rungs.get(n, [])
                    if len(rungs) < 2:
                        continue
                    for k in (size, -size):
                        j = int(round(self.lad.index(n, c.params.get(n)))) + k
                        if not 0 <= j < len(rungs):
                            continue
                        m = self.assemble(dict(c.params, **{n: rungs[j]}))
                        if m is None or m.key in seen or m.params.get(n) != rungs[j] or not self.in_band(heads[i], m):
                            continue
                        seen.add(m.key)
                        (mine_cheap if n in cheap_names else mine_costly).append((i, m))
            costly.append(mine_costly)
            cheap.append(mine_cheap)
        out: List[Tuple[int, Config]] = []
        runs = 0
        while runs < room and any(costly):
            for mine in costly:
                if mine and runs < room:
                    out.append(mine.pop(0))
                    runs += 1
        centres_used = {i for i, _ in out} or set(range(len(centres)))
        for i in sorted(centres_used):
            out += cheap[i][:3 * room // max(1, len(centres_used))]
        return out

    def on_hills(self, p: Dict[str, object]) -> bool:
        return all(self.hill_of(n, p.get(n)) is not None for n in self.spans)

    def island_key(self, p: Dict[str, object]) -> tuple:
        """An island: pcr_indel_model and the hill of each setting."""
        return (p.get(PCR),) + tuple(self.hill_of(n, p.get(n)) for n in sorted(self.spans))

    def hill_of(self, n: str, v) -> Optional[int]:
        """Which hill of setting n holds value v (None: outside every hill, or no hills known)."""
        x = self.lad.index(n, v)
        for k, (_, lo, hi) in enumerate(self.spans.get(n, ())):
            if self.lad.index(n, lo) - 1e-9 <= x <= self.lad.index(n, hi) + 1e-9:
                return k
        return None

    def apart(self, p: Dict[str, object], q: Dict[str, object]) -> bool:
        """Two configs are different islands when pcr_indel_model differs or they sit on different hills of a
        setting. Without hills (Step A cut before A2) a setting counts SEP rungs as apart."""
        for n in INFLUENTIAL:
            if n not in self.lad.rungs:
                continue
            if n == PCR:
                if p.get(n) != q.get(n):
                    return True
            elif n in self.spans:
                if self.hill_of(n, p.get(n)) != self.hill_of(n, q.get(n)):
                    return True
            elif abs(self.lad.index(n, p.get(n)) - self.lad.index(n, q.get(n))) >= SEP:
                return True
        return False

    def bands_for(self, head: Config) -> Dict[str, Tuple[object, object]]:
        """The hill spans an island headed by `head` climbs in."""
        out = {}
        for n in self.spans:
            k = self.hill_of(n, head.params.get(n))
            if k is not None:
                out[n] = self.spans[n][k][1:]
        return out

    def gap(self, n: str, p: Dict[str, object], q: Dict[str, object]) -> float:
        if n == PCR:
            return 0.0 if p.get(n) == q.get(n) else float("inf")
        return abs(self.lad.index(n, p.get(n)) - self.lad.index(n, q.get(n)))

    def top_configs(self) -> List[Tuple[Config, float]]:
        """The best scored configs, each at least islands.seat_apart settings from those kept before it (TOP)."""
        apart = max(1, int(tune("islands.seat_apart")))
        out: List[Tuple[Config, float]] = []
        for c, s in sorted(self.scores.values(), key=lambda t: -t[1]):
            if all(sum(c.params.get(n) != d.params.get(n) for n in self.space.names) >= apart for d, _ in out):
                out.append((c, s))
            if len(out) == TOP:
                break
        return out

    def heads(self, segs: List[Tuple[Config, float]], k: int) -> List[Tuple[Config, float]]:
        """Up to k island heads, best first, whose bands do not overlap: another pcr_indel_model or another hill of
        a setting (without hills: more than 2 * BAND rungs apart)."""
        if int(tune("islands.global_a_heads")):          # Build B2: the top configs as they are
            return list(segs[:k])
        out: List[Tuple[Config, float]] = []
        for c, s in segs:
            if all(self.disjoint(c.params, h.params) for h, _ in out):
                out.append((c, s))
            if len(out) == k:
                break
        return out

    def disjoint(self, p: Dict[str, object], q: Dict[str, object]) -> bool:
        for n in INFLUENTIAL:
            if n not in self.lad.rungs:
                continue
            if n in self.spans:
                a, b = self.hill_of(n, p.get(n)), self.hill_of(n, q.get(n))
                if a != b and a is not None and b is not None:
                    return True
            elif self.gap(n, p, q) > 2 * BAND:
                return True
        return False

    def in_band(self, head: Config, c: Config) -> bool:
        for n in INFLUENTIAL:
            if n not in self.lad.rungs:
                continue
            if n in self.spans:
                if self.hill_of(n, c.params.get(n)) != self.hill_of(n, head.params.get(n)):
                    return False
            elif self.gap(n, head.params, c.params) > BAND:
                return False
        return True

    def distinct(self, ranked: List[Tuple[Config, float]]) -> List[Tuple[Config, float]]:
        out: List[Tuple[Config, float]] = []
        for c, s in ranked:
            if all(self.apart(c.params, d.params) for d, _ in out):
                out.append((c, s))
            if len(out) == TOP:
                break
        return out

    def segments(self, pairs: List[Tuple[Config, float]]) -> List[Tuple[Config, float]]:
        segs = self.distinct(sorted(pairs, key=lambda t: -t[1]))
        self.report["segments"] = [{"key": c.key, "score": round(s, 4),
                                    "influential": {n: c.params.get(n) for n in INFLUENTIAL}} for c, s in segs]
        self.report["done_minute"] = round((self.now() - self.t0) / 60, 1)
        self.publish()
        log.info("global A: %d segments: %s", len(segs),
                 ", ".join(f"{c.key[:8]} {s:.4f} {c.params.get(PCR)}/{c.params.get(MP)}" for c, s in segs[:6]))
        return segs
