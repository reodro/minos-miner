"""Step A: candidates (search/v9/DESIGN.md, Step A).

    wave 0      the defaults, the official seed, the portfolio (up to candidates.portfolio_share of the workers: the
                round's read-length profile first, deduplicated, never the round being replayed) and a Latin hypercube
                for the rest; each point carries candidates.riders calling-threshold riders (one GATK run)
    3 GA gens   parents from three groups of the scored configs, chosen far apart by callset distance:
                  both  the Pareto front on (core up, false positives down)
                  core  highest core, within candidates.core_band of the best core
                  fp    fewest false positives, core within candidates.core_band of the best core
                generations 1-2 breed within each group (a third of the children each); generation 3 crosses
                core x fp (half) and both x (core or fp) (half). Uniform crossover; candidates.mutate_frac of the
                children get candidates.mutate_axes settings redrawn anywhere on their ladder
    pick        n = max(1, workers // candidates.workers_per_candidate): the best score first, then the best of the
                fp, core and both groups in turn, each at least candidates.min_distance (callset distance) from the
                ones taken

A generation launches once candidates.launch_back of the previous wave's points are back; Step A ends when the last
generation is that far back. candidates.max_min (30 minutes, scaled to the round) is only a safety limit: the
candidates are never picked from a generation whose results are not back.
"""
from __future__ import annotations

import asyncio
import logging
from typing import Dict, List, Optional, Sequence, Tuple

from common.tuning import tune
from common.types import FULL, Config
from search import portfolio
from search.axes import latin_hypercube
from search.space import CALL_CONF_GRID

log = logging.getLogger(__name__)

OFFICIAL_SEED = {
    "min_base_quality_score": 16, "min_mapping_quality_score": 46, "base_quality_score_threshold": 22,
    "standard_min_confidence_threshold_for_calling": 80.0, "pcr_indel_model": "NONE", "min_pruning": 2,
    "max_alternate_alleles": 6, "min_dangling_branch_length": 20, "recover_all_dangling_branches": False,
    "max_num_haplotypes_in_population": 128, "active_probability_threshold": 0.002, "min_assembly_region_size": 50,
    "max_assembly_region_size": 270, "assembly_region_padding": 100, "pair_hmm_gap_continuation_penalty": 7,
    "phred_scaled_global_read_mismapping_rate": 45, "heterozygosity": 0.002154, "indel_heterozygosity": 0.000125,
    "contamination_fraction_to_filter": 0.04, "max_reads_per_alignment_start": 116,
}
CONF = "standard_min_confidence_threshold_for_calling"


class Point:
    __slots__ = ("config", "score", "core", "fp")

    def __init__(self, config: Config, score: float, core: float, fp: float):
        self.config, self.score, self.core, self.fp = config, score, core, fp


class StepA:
    def __init__(self, ctx):
        self.ctx = ctx
        self.waves: List[dict] = []            # what each wave sent and got back, for the notes
        self.picked: List[dict] = []           # the candidates Step B starts from

    # ---- the scored pool --------------------------------------------------------------------------------------
    def pool(self) -> List[Point]:
        out = []
        for (key, fid), r in self.ctx.client.latest.items():
            if fid != FULL or r.score is None or r.score <= 0 or r.gate_passed is False or r.core is None:
                continue
            c = self.ctx.client.configs.get(key)
            if c is not None:
                out.append(Point(c, r.score, r.core, float(r.region_fp)))
        return out

    def diverse(self, ranked: Sequence[Point], n: int) -> List[Point]:
        """Up to n of `ranked` (best first), greedily farthest from those taken (callset distance), among the best 4n."""
        cand = list(ranked[:4 * n])
        if not cand:
            return []
        taken = [cand.pop(0)]
        while cand and len(taken) < n:
            i = max(range(len(cand)), key=lambda j: (min(self.ctx.dist.between(cand[j].config, t.config)[0]
                                                         for t in taken), -j))
            taken.append(cand.pop(i))
        return taken

    def groups(self) -> Dict[str, List[Point]]:
        pts = self.pool()
        if not pts:
            return {"both": [], "core": [], "fp": []}
        size, band = int(tune("candidates.group_size")), float(tune("candidates.core_band"))
        best_core = max(p.core for p in pts)
        near = [p for p in pts if p.core >= best_core - band]
        return {"both": self.diverse(pareto_layers(pts, size), size),
                "core": self.diverse(sorted(near, key=lambda p: (-p.core, p.fp)), size),
                "fp": self.diverse(sorted(near, key=lambda p: (p.fp, -p.core)), size)}

    # ---- breeding -----------------------------------------------------------------------------------------------
    def child(self, a: Config, b: Config) -> Optional[Config]:
        ctx = self.ctx
        p = dict(a.params)
        for n in ctx.space.names:
            if ctx.rng.random() < 0.5:
                p[n] = b.params.get(n, p.get(n))
        if ctx.rng.random() < float(tune("candidates.mutate_frac")):
            movable = [n for n in ctx.space.names if ctx.ladders.movable(n)]
            for n in ctx.rng.sample(movable, min(int(tune("candidates.mutate_axes")), len(movable))):
                p[n] = ctx.rng.choice(ctx.ladders.rungs[n])
        return ctx.space._assemble(p)

    def breed(self, gen: int, n: int) -> List[Config]:
        ctx = self.ctx
        g = self.groups()
        pairs: List[Tuple[List[Point], List[Point]]] = []
        if gen < 3:
            live = [x for x in (g["both"], g["core"], g["fp"]) if len(x) >= 2]
            plan = [(x, x) for x in live]
        else:
            other = g["core"] + g["fp"]
            plan = [(g["core"], g["fp"]), (g["both"], other)]
            plan = [(x, y) for x, y in plan if x and y]
        if not plan:                                           # too few scored configs: breed from the whole pool
            pool = sorted(self.pool(), key=lambda p: -p.score)[:40]
            plan = [(pool, pool)] if len(pool) >= 2 else []
        if not plan:
            return []
        out, seen, tries = [], set(), 0
        while len(out) < n and tries < 20 * n:
            tries += 1
            x, y = plan[len(out) % len(plan)]
            a = ctx.rng.choice(x)
            b = ctx.rng.choice([q for q in y if q.config.key != a.config.key] or y)
            c = self.child(a.config, b.config)
            if c is None or c.key in seen or not ctx.fresh(c):
                continue
            seen.add(c.key)
            out.append(c)
        ctx.say("generation %d: %d children (groups both %d, core %d, fp %d)", gen, len(out), len(g["both"]),
                len(g["core"]), len(g["fp"]))
        return out

    # ---- waves --------------------------------------------------------------------------------------------------
    def seeds(self) -> List[Config]:
        ctx = self.ctx
        out: List[Config] = []
        for p in ({}, OFFICIAL_SEED):
            c = ctx.space._assemble(dict(p))
            if c is not None and all(c.key != x.key for x in out):
                out.append(c)
        want = int(ctx.cap * float(tune("candidates.portfolio_share")))
        n = 0
        for c in portfolio.pick(ctx.space, ctx.profile, want, exclude_round=ctx.client.round.tag):
            if all(c.key != x.key for x in out):
                out.append(c)
                n += 1
        ctx.say("seeds: defaults, the official seed and %d from the portfolio (profile %s)", n, ctx.profile or "unknown")
        return out

    def riders(self, c: Config, k: int) -> List[Config]:
        if k <= 0 or CONF not in self.ctx.space.names:
            return []
        out = []
        for v in CALL_CONF_GRID:
            if len(out) >= k:
                break
            if v != c.params.get(CONF):
                r = self.ctx.space._assemble({**c.params, CONF: v})
                if r is not None and self.ctx.fresh(r):
                    out.append(r)
        return out

    async def wave(self, heads: List[Config], riders: int, label: str, until: float) -> None:
        """Send a wave and return once candidates.launch_back of its points are back (or at `until`)."""
        ctx = self.ctx
        futs = []
        for c in heads:
            if not ctx.fresh(c):
                continue
            batch = [c] + self.riders(c, riders)
            fs = ctx.propose(batch, "A", ["big"] * len(batch))
            futs.append(fs[0])
        top = ctx.best()
        w = {"wave": len(self.waves), "label": label, "sent": len(futs), "back": 0, "started_s": round(ctx.elapsed(), 1),
             "best_at_launch": None if top is None else round(top[1], 4), "futs": futs}
        self.waves.append(w)
        ctx.say("wave %d (%s) at t+%.1fm: %d points", w["wave"], label, ctx.elapsed() / 60, len(futs))
        await ctx.wait_share(futs, float(tune("candidates.launch_back")), until)

    async def run(self) -> List[Config]:
        ctx = self.ctx
        until = ctx.t0 + ctx.T(60 * float(tune("candidates.max_min")))
        wave0 = self.seeds()
        wave0 += [c for c in latin_hypercube(ctx.space, ctx.rng, max(0, ctx.cap - len(wave0)))
                  if all(c.key != x.key for x in wave0)]
        await self.wave(wave0, int(tune("candidates.riders")), "seeds + portfolio + hypercube", until)
        for gen in range(1, int(tune("candidates.generations")) + 1):
            if ctx.clock() >= until:
                ctx.say("step A: safety limit (candidates.max_min) before generation %d", gen)
                break
            await self.wave(self.breed(gen, ctx.cap), 0, f"generation {gen}", until)
        return self.pick()

    def wave_rows(self) -> List[dict]:
        return [{**{k: v for k, v in w.items() if k != "futs"}, "back": sum(1 for f in w["futs"] if f.done())}
                for w in self.waves]

    def pick(self) -> List[Config]:
        ctx = self.ctx
        n = max(1, ctx.cap // int(tune("candidates.workers_per_candidate")))
        pts = sorted(self.pool(), key=lambda p: -p.score)
        if not pts:
            ctx.say("step A: nothing scored; starting from the defaults")
            c = ctx.space._assemble({})
            return [c] if c is not None else []
        g = self.groups()
        order = [pts[0]] + [p for i in range(max(len(v) for v in g.values()) if g else 0)
                            for name in ("fp", "core", "both") for p in g[name][i:i + 1]] + pts[1:]
        mc, ms = float(tune("candidates.min_distance")), float(tune("candidates.min_setting_distance"))
        taken: List[Point] = []
        for p in order:
            if len(taken) >= n:
                break
            if all(p.config.key != t.config.key and ctx.dist.far(p.config, t.config, mc, ms) for t in taken):
                taken.append(p)
        self.picked = [{"key": t.config.key, "score": round(t.score, 4), "core": round(t.core, 4), "fp": t.fp,
                        "group": next((name for name in ("both", "core", "fp") if any(q.config.key == t.config.key
                                                                                       for q in g[name])), "score")}
                       for t in taken]
        ctx.say("step A: %d candidate(s) from %d scored configs: %s", len(taken), len(pts),
                ", ".join(f"{t.config.key[:8]} {t.score:.4f} (core {t.core:.4f}, fp {t.fp:.0f})" for t in taken))
        return [t.config for t in taken]


def pareto_layers(pts: Sequence[Point], at_least: int) -> List[Point]:
    """The Pareto front on (core up, fp down), then the next fronts, until at least `at_least` points; best score
    first within the result."""
    left = list(pts)
    out: List[Point] = []
    while left and len(out) < at_least:
        left.sort(key=lambda p: (-p.core, p.fp))
        front, rest, best_fp = [], [], float("inf")
        for p in left:
            if p.fp < best_fp:
                front.append(p)
                best_fp = p.fp
            else:
                rest.append(p)
        out += front
        left = rest
    return sorted(out, key=lambda p: -p.score)
