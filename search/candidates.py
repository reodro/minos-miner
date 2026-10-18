"""Step A: candidates (search/DESIGN.md, Step A).

    wave 0      the defaults, the official seed, validator-edge seeds (min_pruning 1, active 0.0001,
                haplotypes 256/512, calling 0/20, contamination 0.15, ploidy 1/3, soft-clips on),
                the portfolio (up to candidates.portfolio_share of the
                workers: the round's read-length profile first, deduplicated, never the round being replayed)
                and a Latin hypercube for the rest; each point carries candidates.riders calling-threshold
                riders (one GATK run).
                Wave 0 ends once half of it is back and the best score has stood still for 90 seconds, at the
                latest at candidates.launch_back (candidates.wave0_settle 1, the default since 2026-10-18: on round
                05:54 of that day the best was set at minute 8 and the 80% wait ran to minute 14.6). 0: the 80% wait.
    1 GA gen    60 parents, the best scores that are pairwise at least candidates.min_distance apart.
                Three groups are built inside those 60:
                  both  the Pareto front on (core up, false positives down)
                  core  highest core, within candidates.core_band of the best core
                  fp    fewest false positives, core within candidates.core_band of the best core
                generation 1 breeds within each group. MINOS_RF_UCB adds a forest share of that generation.
                Generations 2 and 3 are not run.
    pick        the best score first, then the best of the fp, core and both groups in turn, each at least
                candidates.min_distance (callset distance) from the ones taken. The last opening seat is the
                edge family (validator-legal values the published API does not propose), even when it is
                up to 20 behind. Family A (PCR NONE/CONSERVATIVE × large assembly) is a local opt ~1.6
                below the king on 20261010T044200Z; a PCR flip on this shelf is still this table.
                Step B climbs the first n = max(1, workers // candidates.workers_per_candidate) at once.
"""
from __future__ import annotations

import asyncio
import json
import logging
import math
import os
import time
from typing import Dict, List, Optional, Sequence, Tuple

from common.tuning import tune
from common.types import FULL, Config
from search import portfolio
from search.axes import latin_hypercube
from search.space import CALL_CONF_GRID

log = logging.getLogger(__name__)

PARENTS = 60          # parents are the best scores that stay at least the seat distance apart
SEAT_GAP = 0.10        # a near seat must be within this of the best; the far seat (B27) may sit up to 1.0 below
FAR_SEAT_FLOOR = 1.0   # same cap as climb.REPLACEMENT_FLOOR: a far seat further than this is not opened
FAMILY_SEAT_FLOOR = 20.0  # edge seeds can open well below this shelf; 3.0 dropped family A, 1.0 drops everything
ELITE_MIN_DIFFS = 4    # same as climb.SIMILAR_SETTINGS: an elite-anchored child must leave the elite's hill
# Published API box. A reserved-seat config sits outside it on at least one of these axes.
API_MIN_PRUNING = 2
API_MIN_DANGLING = 2
API_MIN_ACTIVE = 0.001
API_MAX_HAPLOTYPES = 128
API_MAX_CONTAM = 0.05
API_MIN_CALLING = 30.0
API_MAX_AR = 700
API_MIN_READS = 25
API_PLOIDY = 2
# Scouts snap onto validator-edge rungs, not family A (PCR NONE/CONSERVATIVE × large assembly).
EDGE_SCOUT_CHOICES = {
    "min_pruning": (1,),
    "min_dangling_branch_length": (1,),
    "active_probability_threshold": (0.0001, 0.0002, 0.0005),
    "max_num_haplotypes_in_population": (160, 192, 256, 384, 512),
    "contamination_fraction_to_filter": (0.08, 0.15, 0.25, 0.5),
    "standard_min_confidence_threshold_for_calling": (0.0, 10.0, 20.0),
    "max_assembly_region_size": (800, 900, 1000),
    "max_reads_per_alignment_start": (0,),
    "sample_ploidy": (1, 3),
    "dont_use_soft_clipped_bases": (True,),
}
EDGE_SCOUT_AXES = 4


def n_seats(ctx) -> int:
    """How many hills climb at once. Rank 1 needs the edge family beside the best B, so two seats
    whenever there are at least two workers. `workers_per_candidate` only adds seats above that."""
    cap = max(1, int(getattr(ctx, "cap", 1) or 1))
    n = max(1, cap // int(tune("candidates.workers_per_candidate")))
    if cap >= 2:
        return max(2, n)
    return n
STEP_B_MUTATE_START = 0.3   # share of Step B children mutated at the start of the climb
STEP_B_MUTATE_END = 0.7     # that share at the end of Step B; linear in between
STEP_C_MUTATE = (0.3, 0.6)  # Step C generation 1, then generation 2


def step_b_mutation(ctx) -> Optional[float]:
    """Step B only. 0.3 when the climb starts, 0.7 when Step B's clock ends, linear in between.

    None outside Step B. Step C sets its own share: generation 1 is 0.3, generation 2 is 0.6.
    """
    if not str(getattr(ctx, "phase", "")).startswith("B"):
        return None
    span = getattr(ctx, "b_span", None)
    if not span:
        return STEP_B_MUTATE_START
    start, end = span
    if end <= start:
        return STEP_B_MUTATE_END
    p = (ctx.clock() - start) / (end - start)
    if p < 0:
        p = 0.0
    elif p > 1:
        p = 1.0
    return STEP_B_MUTATE_START + (STEP_B_MUTATE_END - STEP_B_MUTATE_START) * p


def mutation_share(ctx) -> float:
    """The share of children that move one rung. Step B ramps 0.3 to 0.7. Step C generation 1 is 0.3 and generation 2 is 0.6."""
    frac = step_b_mutation(ctx)
    if frac is not None:
        return frac
    gen = getattr(ctx, "step_c_gen", None)
    if gen == 1:
        return STEP_C_MUTATE[0]
    if gen == 2:
        return STEP_C_MUTATE[1]
    return float(tune("candidates.mutate_frac"))


def _mutate_note(ctx) -> str:
    if step_b_mutation(ctx) is None and not getattr(ctx, "step_c_gen", None):
        return ""
    return f", mutation {mutation_share(ctx):.2f}"

# The configuration that held a high score on 20261010T044200Z. It is the Step A
# metric: after the interval is chosen from scored configs, record whether this
# configuration lies inside it. It is not a seed and it does not set the interval.
# The edge-family seeds below are built on it, so a replay of REFERENCE_ROUND skips them.
REFERENCE_ROUND = "20261010T044200Z"
REFERENCE_PEAK = {
    "min_base_quality_score": 24, "min_mapping_quality_score": 46, "base_quality_score_threshold": 33,
    "standard_min_confidence_threshold_for_calling": 80.0, "pcr_indel_model": "AGGRESSIVE", "min_pruning": 2,
    "max_alternate_alleles": 3, "min_dangling_branch_length": 20, "recover_all_dangling_branches": False,
    "max_num_haplotypes_in_population": 91, "active_probability_threshold": 0.022558,
    "min_assembly_region_size": 70, "max_assembly_region_size": 102, "assembly_region_padding": 15,
    "pair_hmm_gap_continuation_penalty": 17, "phred_scaled_global_read_mismapping_rate": 31,
    "heterozygosity": 0.005012, "indel_heterozygosity": 0.000835, "contamination_fraction_to_filter": 0.01875,
    "max_reads_per_alignment_start": 25,
}

OFFICIAL_SEED = {
    "min_base_quality_score": 16, "min_mapping_quality_score": 46, "base_quality_score_threshold": 22,
    "standard_min_confidence_threshold_for_calling": 80.0, "pcr_indel_model": "NONE", "min_pruning": 2,
    "max_alternate_alleles": 6, "min_dangling_branch_length": 20, "recover_all_dangling_branches": False,
    "max_num_haplotypes_in_population": 128, "active_probability_threshold": 0.002, "min_assembly_region_size": 50,
    "max_assembly_region_size": 270, "assembly_region_padding": 100, "pair_hmm_gap_continuation_penalty": 7,
    "phred_scaled_global_read_mismapping_rate": 45, "heterozygosity": 0.002154, "indel_heterozygosity": 0.000125,
    "contamination_fraction_to_filter": 0.04, "max_reads_per_alignment_start": 116,
}
def _edge_seed(**over) -> dict:
    return {**REFERENCE_PEAK, **over}


# Validator-legal values the published API does not propose. Family A (PCR NONE/CONSERVATIVE ×
# large assembly) is a local opt on this shelf; those seeds stay in wave 0 via OFFICIAL_SEED only.
INDEL_FAMILY_SEEDS = (
    _edge_seed(min_pruning=1),
    _edge_seed(active_probability_threshold=0.0001),
    _edge_seed(max_num_haplotypes_in_population=256),
    _edge_seed(max_num_haplotypes_in_population=512),
    _edge_seed(contamination_fraction_to_filter=0.15),
    _edge_seed(standard_min_confidence_threshold_for_calling=20.0),
    _edge_seed(standard_min_confidence_threshold_for_calling=0.0),
    _edge_seed(min_dangling_branch_length=1),
    _edge_seed(max_assembly_region_size=1000),
    _edge_seed(sample_ploidy=1),
    _edge_seed(sample_ploidy=3),
    _edge_seed(dont_use_soft_clipped_bases=True),
    _edge_seed(max_reads_per_alignment_start=0),
    _edge_seed(min_pruning=1, active_probability_threshold=0.0001,
               max_num_haplotypes_in_population=256, dont_use_soft_clipped_bases=True,
               standard_min_confidence_threshold_for_calling=20.0),
)


def in_indel_family(c: Config) -> bool:
    """The reserved far-seat box: at least one setting outside the published API box."""
    p = c.params
    pruning = p.get("min_pruning")
    if pruning is not None and pruning < API_MIN_PRUNING:
        return True
    dangling = p.get("min_dangling_branch_length")
    if dangling is not None and dangling < API_MIN_DANGLING:
        return True
    active = p.get("active_probability_threshold")
    if active is not None and active < API_MIN_ACTIVE - 1e-12:
        return True
    hap = p.get("max_num_haplotypes_in_population")
    if hap is not None and hap > API_MAX_HAPLOTYPES:
        return True
    contam = p.get("contamination_fraction_to_filter")
    if contam is not None and contam > API_MAX_CONTAM + 1e-12:
        return True
    calling = p.get("standard_min_confidence_threshold_for_calling")
    if calling is not None and calling < API_MIN_CALLING - 1e-9:
        return True
    mx = p.get("max_assembly_region_size")
    if mx is not None and mx > API_MAX_AR:
        return True
    reads = p.get("max_reads_per_alignment_start")
    if reads is not None and reads < API_MIN_READS:
        return True
    ploidy = p.get("sample_ploidy")
    if ploidy is not None and ploidy != API_PLOIDY:
        return True
    return p.get("dont_use_soft_clipped_bases") is True


def left_the_rim(a: Config, b: Config) -> bool:
    """True when `a` is the edge family and `b` is still inside the published API box.

    A CONSERVATIVE clone of the AGGRESSIVE peak (gap 26, padding 5) is still this table.
    Family A is a measured local opt. The reserved seat is the validator-edge box.
    """
    return in_indel_family(a) and not in_indel_family(b)


ELITE_AXES = (
    "max_assembly_region_size", "phred_scaled_global_read_mismapping_rate", "base_quality_score_threshold",
    "standard_min_confidence_threshold_for_calling", "min_pruning", "min_mapping_quality_score",
    "active_probability_threshold", "min_assembly_region_size",
)
ELITE_PAIR = ("active_probability_threshold", "standard_min_confidence_threshold_for_calling")
ELITE_BIAS_SHARE = 0.40
ELITE_BAND = 1.0
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
        """Scored configs the breeder may use as parents.

        The wave-0 peak interval is measured and logged, not enforced: `_select_range` clears `_bound_pool`.
        Breeding inside it froze 15 of 19 settings on 044200Z (2026-10-16); 2,500 children were known callsets,
        generation 1 sent 32 and Step B had one seat.
        """
        bound = getattr(self, "_bound_pool", True)
        out = []
        for (key, fid), r in self.ctx.client.latest.items():
            if fid != FULL or r.score is None or r.score <= 0 or r.gate_passed is False or r.core is None:
                continue
            c = self.ctx.client.configs.get(key)
            if c is None or (bound and not self._inside(c)):
                continue
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

    async def _breathe(self) -> None:
        """Yield when this step has run for 50 ms, so the round loop can poll the fleet and write the dashboard."""
        now = time.monotonic()
        if now - getattr(self, "_yield_at", 0.0) > 0.05:
            self._yield_at = now
            await asyncio.sleep(0)

    async def far_parents(self, n: int) -> List[Point]:
        """The `n` best scores that are at least the seat distance apart from each other.

        The plain top 20 are one hill: climber grids and tied callsets fill it, and Step B then breeds children
        next to the climbers. A later score is kept only when it is far from every parent already kept. Callset
        distance when both have a sketch, otherwise the mean setting difference. Reused until the best score
        changes (v9.18). A flood of new scores at the same best must not rebuild this list: that scan used to
        freeze leftover GA for about a minute and drain the fleet.
        """
        ctx = self.ctx
        pts = sorted(self.pool(), key=lambda p: -p.score)
        if not pts or n <= 0:
            return []
        stamp = (n, round(pts[0].score, 6), bool(getattr(self, "_bound_pool", True)))
        cached = getattr(self, "_parents", None)
        if cached is not None and getattr(self, "_parent_stamp", None) == stamp:
            return cached
        dist = getattr(ctx, "dist", None)
        if dist is None:
            chosen = pts[:n]
        else:
            mc = float(tune("candidates.min_distance"))
            ms = float(tune("candidates.min_setting_distance"))
            chosen = []
            for p in pts:
                if all(dist.far(p.config, t.config, mc, ms) for t in chosen):
                    chosen.append(p)
                    if len(chosen) >= n:
                        break
                await self._breathe()
        if chosen and [p.config.key for p in chosen] != [p.config.key for p in (cached or [])]:
            ctx.say("GA parents: %d configs far apart, %.4f down to %.4f",
                    len(chosen), chosen[0].score, chosen[-1].score)
        self._parent_stamp = stamp
        self._parents = chosen
        return chosen

    async def groups(self) -> Dict[str, List[Point]]:
        """The both / core / fp parent groups, reused until the far parents change.

        Each rebuild is about 20,000 callset distances; rebuilt on every GA pass it blocked the round loop for
        seconds (one child per pass on 044200Z, 2026-10-16)."""
        pts = await self.far_parents(PARENTS)
        if not pts:
            return {"both": [], "core": [], "fp": []}
        size, band = int(tune("candidates.group_size")), float(tune("candidates.core_band"))
        stamp = (tuple(p.config.key for p in pts), size, band)
        if getattr(self, "_groups_stamp", None) == stamp:
            return self._groups
        best_core = max(p.core for p in pts)
        near = [p for p in pts if p.core >= best_core - band]
        self._groups = {"both": self.diverse(pareto_layers(pts, size), size),
                        "core": self.diverse(sorted(near, key=lambda p: (-p.core, p.fp)), size),
                        "fp": self.diverse(sorted(near, key=lambda p: (p.fp, -p.core)), size)}
        self._groups_stamp = stamp
        return self._groups

    # ---- breeding -----------------------------------------------------------------------------------------------
    def child(self, a: Config, b: Config) -> Optional[Config]:
        """One child. Each verified group comes from one parent. An unlinked setting comes from either.

        The one-rung mutation then moves only unlinked settings, so a group is not split.
        """
        from search.separate import _LINKED, _cross
        ctx = self.ctx
        crossed = _cross(ctx, a, b)
        if crossed is None:
            return None
        p = dict(crossed.params)
        frac = mutation_share(ctx)
        if ctx.rng.random() < frac:
            movable = [n for n in ctx.space.names if ctx.ladders.movable(n) and n not in _LINKED]
            k = min(int(tune("candidates.mutate_axes")), len(movable))
            for n in ctx.rng.sample(movable, k) if k else ():
                step = ctx.rng.choice((1, -1))
                v = ctx.ladders.shift(n, p.get(n), step)
                if v is None:
                    v = ctx.ladders.shift(n, p.get(n), -step)
                if v is not None:
                    p[n] = v
        child = ctx.space._assemble(p)
        if child is not None and getattr(self, "_bound_pool", True) and not self._inside(child):
            return None
        return child

    def elite_probe(self, elite: Config, mate: Config) -> Optional[Config]:
        """Copy the best config and take one to three of the mate's settings.

        This is the elite phase: find which of the best config's settings to change, then stop.
        A child that rewrites the config is a new hill, and that search does not finish before the grid.
        """
        ctx = self.ctx
        p = dict(elite.params)
        movable = [n for n in ctx.space.names if ctx.ladders.movable(n)]
        differed = [n for n in movable if mate.params.get(n) != elite.params.get(n)]
        names = differed or movable
        if not names:
            return None
        k = ctx.rng.randint(1, min(3, len(names)))
        for n in ctx.rng.sample(names, k):
            mv = mate.params.get(n)
            if n in differed and mv is not None:
                p[n] = mv
                continue
            step = ctx.rng.choice((1, -1))
            v = ctx.ladders.shift(n, p.get(n), step)
            if v is None:
                v = ctx.ladders.shift(n, p.get(n), -step)
            if v is not None:
                p[n] = v
        return ctx.space._assemble(p)

    def elite_child(self, elite: Config, mate: Config) -> Optional[Config]:
        """A child of the current best that is a different hill (at least four searched settings).

        A shortfall is filled by a whole verified group from the mate, then by unlinked settings.
        The extra one-rung move stays on an unlinked setting.
        """
        from search.separate import LINKED_GROUPS, _LINKED
        ctx = self.ctx
        c = self.child(elite, mate)
        if c is None:
            return None
        p = dict(c.params)
        names = list(ctx.space.names)
        movable = [n for n in names if ctx.ladders.movable(n)]
        diffs = lambda: sum(p.get(n) != elite.params.get(n) for n in names)
        if diffs() < ELITE_MIN_DIFFS:
            for group in LINKED_GROUPS:
                if diffs() >= ELITE_MIN_DIFFS:
                    break
                if any(mate.params.get(n) != elite.params.get(n) for n in group):
                    for n in group:
                        mv = mate.params.get(n)
                        if mv is not None:
                            p[n] = mv
            for n in movable:
                if n in _LINKED or diffs() >= ELITE_MIN_DIFFS:
                    continue
                mv, ev = mate.params.get(n), elite.params.get(n)
                if mv is not None and mv != ev and p.get(n) != mv:
                    p[n] = mv
            free = [n for n in movable if n not in _LINKED and p.get(n) == elite.params.get(n)]
            guard = 0
            while diffs() < ELITE_MIN_DIFFS and free and guard < len(free) + 2:
                guard += 1
                n = free.pop(0)
                step = ctx.rng.choice((1, -1))
                v = ctx.ladders.shift(n, p.get(n), step)
                if v is None:
                    v = ctx.ladders.shift(n, p.get(n), -step)
                if v is not None and v != elite.params.get(n):
                    p[n] = v
        extra = [n for n in movable if n not in _LINKED]
        if extra:
            n = ctx.rng.choice(extra)
            step = ctx.rng.choice((1, -1))
            v = ctx.ladders.shift(n, p.get(n), step)
            if v is None:
                v = ctx.ladders.shift(n, p.get(n), -step)
            if v is not None:
                p[n] = v
        return ctx.space._assemble(p)

    async def breed(self, gen: int, n: int, *, why: str = "", accept=None, elite: bool = False) -> List[Config]:
        from search.piles import enabled as piles_enabled
        if piles_enabled():
            return await self.breed_piles(gen, n, why=why, accept=accept)
        ctx = self.ctx
        g = await self.groups()
        if gen < 3:
            live = [x for x in (g["both"], g["core"], g["fp"]) if len(x) >= 2]
            plan = [(x, x) for x in live]
        else:
            other = g["core"] + g["fp"]
            plan = [(g["core"], g["fp"]), (g["both"], other)]
            plan = [(x, y) for x, y in plan if x and y]
        if not plan:                                           # too few scored configs: breed from the whole pool
            pool = await self.far_parents(PARENTS)
            plan = [(pool, pool)] if len(pool) >= 2 else []
        if not plan:
            return []
        frac = float(tune("candidates.elite_frac")) if elite else 0.0
        elite_cfg = None
        if frac > 0:
            top = ctx.best() if callable(getattr(ctx, "best", None)) else None
            if top is not None:
                elite_cfg = top[0]
        want_elite = int(round(n * frac)) if elite_cfg is not None else 0
        out, seen, tries, n_elite = [], set(), 0, 0
        while len(out) < n and tries < 20 * n:
            tries += 1
            x, y = plan[len(out) % len(plan)]
            a = ctx.rng.choice(x)
            b = ctx.rng.choice([q for q in y if q.config.key != a.config.key] or y)
            use_elite = elite_cfg is not None and n_elite < want_elite
            if use_elite:
                mate = b.config if b.config.key != elite_cfg.key else a.config
                if mate.key == elite_cfg.key:
                    await self._breathe()
                    continue
                c = self.elite_child(elite_cfg, mate)
            else:
                c = self.child(a.config, b.config)
            if c is None or c.key in seen or not ctx.fresh(c):
                await self._breathe()
                continue
            verdict = True if accept is None else accept(c)
            if hasattr(verdict, "__await__"):
                verdict = await verdict
            if not verdict:
                await self._breathe()
                continue
            seen.add(c.key)
            out.append(c)
            if use_elite:
                n_elite += 1
        extra = f", {n_elite} elite-anchored" if n_elite else ""
        ctx.say("%s: %d children (groups both %d, core %d, fp %d%s%s)", why or f"generation {gen}", len(out),
                len(g["both"]), len(g["core"]), len(g["fp"]), _mutate_note(ctx), extra)
        return out

    # ---- waves --------------------------------------------------------------------------------------------------
    def seeds(self) -> List[Config]:
        """Defaults, official, the edge family and the portfolio. With islands (search/islands.py) the fixed seeds go
        to island 0 and the portfolio is dealt by rank: island i of k takes ranks i, i + k, i + 2k, ..."""
        ctx = self.ctx
        island = getattr(ctx, "island", None)
        index, k = (island.index, island.k) if island is not None else (0, 1)
        out: List[Config] = []
        own = portfolio.same_round(REFERENCE_ROUND, ctx.client.round.tag)
        fixed = ({}, OFFICIAL_SEED, *(() if own else INDEL_FAMILY_SEEDS)) if index == 0 else ()
        for p in fixed:
            c = ctx.space._assemble(dict(p))
            if c is not None and all(c.key != x.key for x in out):
                out.append(c)
        if own and index == 0:
            ctx.say("seeds: the edge-family seeds are skipped; they are built on this round's best config")
        families = len(out)
        g_want = int(ctx.cap * float(tune("candidates.generalist_share")))
        n_gen = 0
        for c in portfolio.generalists(ctx.space, g_want * k)[index::k]:      # cross-round generalists first
            if all(c.key != x.key for x in out):
                out.append(c)
                n_gen += 1
        if n_gen:
            ctx.say("seeds: %d cross-round generalists first", n_gen)
        want = int(ctx.cap * float(tune("candidates.portfolio_share")))
        ranked = list(portfolio.pick(ctx.space, ctx.profile, want * k, exclude_round=ctx.client.round.tag))
        n = 0
        for c in ranked[index::k]:
            if all(c.key != x.key for x in out):
                out.append(c)
                n += 1
        if index == 0:
            ctx.say("seeds: defaults, official, %d edge-family seeds and %d from the portfolio (profile %s)",
                    max(0, families - 2), n, ctx.profile or "unknown")
        else:
            ctx.say("seeds: island %d of %d: %d from the portfolio (ranks %d, %d, %d, ...; profile %s)",
                    index, k, n, index + 1, index + 1 + k, index + 1 + 2 * k, ctx.profile or "unknown")
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

    async def wave(self, heads: List[Config], riders: int, label: str, until: float, tag: str = "A",
                   settle: bool = False, discovery: bool = False) -> None:
        """Send a wave and return once candidates.launch_back of its points are back (or at `until`)."""
        ctx = self.ctx
        futs = []
        for c in heads:
            if not ctx.fresh(c):
                continue
            batch = [c] + self.riders(c, riders)
            fs = ctx.propose(batch, tag, ["big"] * len(batch))
            futs.append(fs[0])
        top = ctx.best()
        w = {"wave": len(self.waves), "label": label, "sent": len(futs), "back": 0, "started_s": round(ctx.elapsed(), 1),
             "best_at_launch": None if top is None else round(top[1], 4), "futs": futs}
        self.waves.append(w)
        ctx.say("wave %d (%s) at t+%.1fm: %d points", w["wave"], label, ctx.elapsed() / 60, len(futs))
        if settle:
            await self._wait_opening(futs, until)
        elif discovery:
            await self._wait_discovery(futs, until)
        else:
            await ctx.wait_share(futs, float(tune("candidates.launch_back")), until)

    async def _wait_opening(self, futs, until: float) -> None:
        """End the opening once the best score has stood still, so Step B is not delayed for a marginal gain.

        At least half the heads must be back. If the best is still rising, wait until candidates.launch_back.
        """
        ctx = self.ctx
        quiet = 90.0
        floor = math.ceil(0.50 * len(futs)) if futs else 0
        cap = math.ceil(float(tune("candidates.launch_back")) * len(futs)) if futs else 0
        pending = {f for f in futs if not f.done()}
        still_since = None
        last_best = None
        while pending:
            done = len(futs) - len(pending)
            if ctx.clock() >= until or ctx.over():
                ctx.say("step A opening stops at the time limit, %d of %d back", done, len(futs))
                break
            top = ctx.best()
            score = None if top is None else round(top[1], 4)
            now = ctx.clock()
            if score != last_best:
                last_best = score
                still_since = now
            elif still_since is not None and done >= floor and now - still_since >= quiet:
                ctx.say("step A opening stops: best %.4f still for %.0fs, %d of %d back",
                        score, now - still_since, done, len(futs))
                break
            if done >= cap:
                ctx.say("step A opening stops at %.0f%% back, best %s",
                        100 * float(tune("candidates.launch_back")),
                        "unset" if score is None else f"{score:.4f}")
                break
            _, pending = await asyncio.wait(pending, timeout=1.0, return_when=asyncio.FIRST_COMPLETED)
            pending = {f for f in pending if not f.done()}
            ctx.publish()

    def _region_reps(self, reps: List[Config], best: float) -> List[Config]:
        """Promising configs at least 4 settings from every region already counted."""
        names = list(self.ctx.space.names)
        added = []
        for config, score in self.ctx.scored():
            if score < best - 0.10:
                continue
            if any(sum(config.params.get(name) != rep.params.get(name) for name in names) < ELITE_MIN_DIFFS
                   for rep in reps + added):
                continue
            added.append(config)
        return added

    async def _wait_discovery(self, futs, until: float) -> None:
        """End Generation 1 when two consecutive 2-minute windows add no new region and lift the best by under 0.10.

        A window with no finished evaluations does not count. Wave 0 has already met the coverage minimum.
        The 2-minute grain is the clock at which historical rounds have only a median of 39 genetic results,
        so one quiet window is not a handoff.
        """
        ctx = self.ctx
        window = 120.0
        ctx.say("discovery handoff: windows of %.0fs, stop after 2 dry windows; a dry window adds no region "
                "4 settings from those already counted and lifts the best by under 0.10", window)
        top = ctx.best()
        best = 0.0 if top is None else float(top[1])
        reps = self._region_reps([], best)
        ctx.say("discovery handoff: %d regions at the start of generation 1, best %.4f", len(reps), best)
        pending = {f for f in futs if not f.done()}
        window_started = ctx.clock()
        window_best = best
        window_done = len(futs) - len(pending)
        dry = 0
        while pending and ctx.clock() < until and not ctx.over():
            _, pending = await asyncio.wait(pending, timeout=1.0, return_when=asyncio.FIRST_COMPLETED)
            pending = {f for f in pending if not f.done()}
            ctx.publish()
            if ctx.clock() - window_started < window:
                continue
            done = len(futs) - len(pending)
            arrived = done - window_done
            if arrived == 0:
                ctx.say("discovery window: no finished evaluations, dry streak stays %d", dry)
            else:
                top = ctx.best()
                best = window_best if top is None else float(top[1])
                gain = best - window_best
                added = self._region_reps(reps, best)
                reps.extend(added)
                quiet = not added and gain < 0.10
                dry = dry + 1 if quiet else 0
                ctx.say("discovery window: %d finished, %d new regions, best %+.4f, dry streak %d",
                        arrived, len(added), gain, dry)
            window_started = ctx.clock()
            window_best = best
            window_done = done
            if dry >= 2:
                ctx.say("discovery handoff: yield stayed low for two windows, best %.4f, %d regions",
                        best, len(reps))
                return
        top = ctx.best()
        best = best if top is None else float(top[1])
        ctx.say("discovery handoff: generation finished on the clock, best %.4f, %d regions", best, len(reps))

    def _reference_fit(self, interval: dict) -> tuple:
        """Whether the 20261010 configuration lies in `interval`, and which settings miss it."""
        centre = self.ctx.space._assemble(dict(REFERENCE_PEAK))
        if centre is None or not interval:
            return False, ["unassembled"]
        misses = []
        for name, vals in interval.items():
            if not vals:
                continue
            value = centre.params.get(name)
            if any(value == v for v in vals):
                continue
            nums = [v for v in vals if isinstance(v, (int, float)) and not isinstance(v, bool)]
            if isinstance(value, (int, float)) and not isinstance(value, bool) and nums:
                if min(nums) - 1e-6 <= float(value) <= max(nums) + 1e-6:
                    continue
            misses.append(name)
        return not misses, misses

    def _inside(self, config: Config) -> bool:
        """True when no peak interval is set, or every setting of `config` lies on that interval."""
        interval = getattr(self, "_range", None)
        if not interval:
            return True
        for name, vals in interval.items():
            if not vals:
                continue
            value = config.params.get(name)
            if any(value == v for v in vals):
                continue
            nums = [v for v in vals if isinstance(v, (int, float)) and not isinstance(v, bool)]
            if isinstance(value, (int, float)) and not isinstance(value, bool) and nums:
                if min(nums) - 1e-6 <= float(value) <= max(nums) + 1e-6:
                    continue
            return False
        return True

    def _select_range(self) -> None:
        """Wave 0 has stopped moving the best. The interval around that best is the parent pool."""
        best, elite, used = self._elite()
        if best is None or not elite:
            self.ctx.say("wave 0 at its limit: nothing scored")
            return
        self._range = self._interval(elite)
        self._parents = None
        self._report_interval("selected", best, elite, used, self._range)
        self._bound_pool = False
        self.ctx.say("wave 0 at its limit: %d scored configs are inside the interval; generation 1 breeds from all",
                     sum(1 for c, _ in self.ctx.scored() if self._inside(c)))

    async def _open_pool_if_thin(self) -> None:
        """Use the full scored set when the peak interval cannot supply two far parents.

        Generation 1 needs a pair. Seats come from the same pool. On 20261012 and 20261014 the interval
        held one parent, generation 1 sent nothing, and the seats sat in a box that later high scores
        were outside of.
        """
        parents = await self.far_parents(2)
        if len(parents) >= 2:
            return
        self._bound_pool = False
        self._parents = None
        self.ctx.say("step A: the peak interval has %d far parent(s); generation 1 and the seats use the full pool",
                     len(parents))

    def _elite(self, margin: float = 0.05):
        """Configs within `margin` of the best score so far. Widens once if that set is a single config."""
        scored = [(c, s) for c, s in self.ctx.scored()]
        if not scored:
            return None, [], 0.0
        best = max(s for _, s in scored)
        elite = [c for c, s in scored if s >= best - margin]
        used = margin
        if len(elite) < 2:
            used = 0.10
            elite = [c for c, s in scored if s >= best - used]
        return best, elite, used

    def _interval(self, elite) -> dict:
        """Smallest ladder span that contains every elite value of each setting."""
        ctx = self.ctx
        out = {}
        for name in ctx.space.names:
            if not ctx.ladders.movable(name):
                continue
            rungs = ctx.ladders.rungs.get(name) or []
            idxs = [ctx.ladders.index(name, c.params.get(name)) for c in elite]
            if not rungs or not idxs:
                vals = []
                for c in elite:
                    v = c.params.get(name)
                    if v not in vals:
                        vals.append(v)
                out[name] = vals
                continue
            lo = max(0, int(math.floor(min(idxs) + 1e-9)))
            hi = min(len(rungs) - 1, int(math.ceil(max(idxs) - 1e-9)))
            out[name] = list(rungs[lo:hi + 1])
        return out

    def _interval_size(self, interval: dict) -> tuple:
        """Share of the ladder grid still inside the interval, and how many settings are still open."""
        volume = 1.0
        axes = 0
        open_axes = 0
        for name, vals in interval.items():
            rungs = self.ctx.ladders.rungs.get(name) or []
            full = len(rungs) if rungs else max(1, len(vals))
            span = max(1, len(vals))
            volume *= span / full
            axes += 1
            if span > 1:
                open_axes += 1
        return volume, open_axes, axes

    def _report_interval(self, label: str, best: float, elite, used: float, interval: dict) -> None:
        open_axes = [(n, vals) for n, vals in interval.items() if len(vals) > 1]
        fixed = sum(1 for vals in interval.values() if len(vals) <= 1)
        self.ctx.say("peak interval %s: best %.4f, %d configs within %.2f, %d settings fixed, %d open",
                     label, best, len(elite), used, fixed, len(open_axes))
        for name, vals in open_axes:
            self.ctx.say("peak interval %s: %s %s .. %s (%d rungs)", label, name, vals[0], vals[-1], len(vals))
        volume, n_open, n_axes = self._interval_size(interval)
        holds = any(self.ctx.score(c) == best for c in elite) if elite else False
        ref_in, misses = self._reference_fit(interval)
        self.ctx.say(
            "step A metric: size %.6g of the grid, %d of %d settings open, scored maximum %s, reference %s%s",
            volume, n_open, n_axes, "inside" if holds else "outside",
            "inside" if ref_in else "outside",
            "" if ref_in else " (" + ", ".join(misses) + ")")

    def _sample_interval(self, interval: dict, centre: Config, n: int) -> List[Config]:
        """One Latin sample inside the peak interval. A fixed setting stays at the centre's value."""
        ctx = self.ctx
        if n <= 0:
            return []
        columns = {}
        for name, vals in interval.items():
            if len(vals) <= 1:
                columns[name] = [vals[0] if vals else centre.params.get(name)] * n
                continue
            strata = list(range(n))
            ctx.rng.shuffle(strata)
            columns[name] = [vals[min(len(vals) - 1, (k * len(vals)) // n)] for k in strata]
        out, seen = [], set()
        for i in range(n):
            params = dict(centre.params)
            for name, col in columns.items():
                params[name] = col[i]
            c = ctx.space._assemble(params)
            if c is not None and c.key not in seen and ctx.fresh(c):
                seen.add(c.key)
                out.append(c)
        return out

    async def _tighten(self) -> None:
        """Second opening wave, only inside the smallest interval that holds the current peak."""
        ctx = self.ctx
        found = self._elite()
        if found[0] is None:
            ctx.say("peak interval: nothing scored")
            return
        best, elite, used = found
        interval = self._interval(elite)
        self._report_interval("opened", best, elite, used, interval)
        centre = max(elite, key=lambda c: ctx.score(c) or 0.0)
        sample = self._sample_interval(interval, centre, ctx.cap)
        if not sample:
            ctx.say("peak interval: the box is already scored")
            return
        ctx.say("peak interval: %d points inside the box", len(sample))
        await self.wave(sample, 0, "peak interval", ctx.clock() + 15 * 60)
        found = self._elite()
        if found[0] is None:
            return
        best, elite, used = found
        self._report_interval("tightened", best, elite, used, self._interval(elite))

    def _same_value(self, a, b) -> bool:
        if a == b:
            return True
        if isinstance(a, bool) or isinstance(b, bool):
            return False
        if isinstance(a, (int, float)) and isinstance(b, (int, float)):
            return abs(float(a) - float(b)) < 1e-9
        return False

    def _graft_candidates(self, margin: float) -> List[Tuple[str, object, Config]]:
        """The peak with one setting moved to a rung whose best score is within `margin` of the peak."""
        ctx = self.ctx
        scored = [(c, s) for c, s in ctx.scored()]
        if not scored:
            return []
        best = max(s for _, s in scored)
        peak = max(scored, key=lambda t: t[1])[0]
        self._graft_peak = peak
        self._graft_best = best
        out, seen = [], {peak.key}
        for name in ctx.space.names:
            if not ctx.ladders.movable(name):
                continue
            groups: List[list] = []
            for config, score in scored:
                value = config.params.get(name)
                if self._same_value(value, peak.params.get(name)):
                    continue
                hit = next((g for g in groups if self._same_value(g[0], value)), None)
                if hit is None:
                    groups.append([value, score])
                elif score > hit[1]:
                    hit[1] = score
            for value, score in groups:
                if score < best - margin:
                    continue
                params = dict(peak.params)
                params[name] = value
                graft = ctx.space._assemble(params)
                if graft is None or graft.key in seen:
                    continue
                seen.add(graft.key)
                out.append((name, graft.params.get(name), graft))
        return out

    async def _run_grafts(self, margin: float) -> None:
        """Score one-setting grafts on the peak. A rung stays when that score is within `margin`."""
        ctx = self.ctx
        grafts = self._graft_candidates(margin)
        pending = [g for g in grafts if ctx.score(g[2]) is None]
        ctx.say("grafts: %d one-setting runs on %s, %d already scored, margin %.2f",
                len(pending), self._graft_peak.key[:8], len(grafts) - len(pending), margin)
        futs = []
        for _name, _value, config in pending:
            futs.append(ctx.propose([config], "A", ["small"])[0])
        if futs:
            await ctx.wait_share(futs, 1.0, ctx.t0 + ctx.length)
        best = self._graft_best
        peak = self._graft_peak
        kept = {name: set() for name in ctx.space.names if ctx.ladders.movable(name)}
        dropped = {name: set() for name in kept}
        for name in kept:
            value = peak.params.get(name)
            kept[name].add(value)
        for name, value, config in grafts:
            score = ctx.score(config)
            if score is None:
                ctx.say("graft %s %s unscored", name, value)
                continue
            stay = score >= best - margin
            (kept if stay else dropped)[name].add(value)
            ctx.say("graft %s %s %.4f (%+.4f) %s", name, value, score, score - best,
                    "stays" if stay else "drops")
        interval = {}
        for name in kept:
            rungs = ctx.ladders.rungs.get(name) or []
            dropped_vals = dropped[name]
            if rungs:
                interval[name] = [v for v in rungs if not any(self._same_value(v, d) for d in dropped_vals)]
            else:
                interval[name] = [v for v in kept[name]]
        self._range = interval
        volume, n_open, n_axes = self._interval_size(interval)
        ref_in, misses = self._reference_fit(interval)
        ctx.say("graft domain: size %.6g of the grid, %d of %d settings open, reference %s%s",
                volume, n_open, n_axes, "inside" if ref_in else "outside",
                "" if ref_in else " (" + ", ".join(misses) + ")")
        open_because = [name for name in misses if name not in {n for n, _v, _c in grafts}]
        if ref_in:
            untested = []
            offered = ctx.space._assemble(dict(REFERENCE_PEAK))
            if offered is not None:
                for name, vals in interval.items():
                    value = offered.params.get(name)
                    tested = any(self._same_value(value, v) for v in kept[name]) or any(
                        self._same_value(value, v) for v in dropped[name])
                    if not tested:
                        untested.append(name)
            if untested:
                ctx.say("graft domain: reference is inside with %d settings untested (%s)",
                        len(untested), ", ".join(untested))

    async def _send_and_wait(self, configs: List[Config]) -> None:
        pending = [config for config in configs if config is not None and self.ctx.score(config) is None]
        if not pending:
            return
        futs = [self.ctx.propose([config], "A", ["small"])[0] for config in pending]
        await self.ctx.wait_share(futs, 1.0, self.ctx.t0 + self.ctx.length)

    async def _run_neighborhood(self) -> None:
        """One wave of ±1 neighbors around the seat, then a capped cube of the faces that beat it.

        The 0.10 span is measured against each line's own maximum. Gain against the seat is recorded
        separately. The known best configuration is scored afterwards as an audit.
        """
        import itertools
        ctx = self.ctx
        short_ladder = 12
        margin = 0.10
        cube_cap = 64
        seat = next((cfg for key, cfg in ctx.client.configs.items() if key.startswith("54d85ee7")), None)
        if seat is None or ctx.score(seat) is None:
            ctx.say("neighborhood: seat 54d85ee7 is not in the scored set")
            return
        incumbent = float(ctx.score(seat))
        ctx.say("neighborhood: seat %s %.4f", seat.key[:8], incumbent)
        ctx.say("neighborhood: span is rungs within %.2f of that line's maximum; gain versus the seat is recorded separately",
                margin)
        ctx.say("neighborhood: ladders of %d rungs or fewer are scored in full; longer ladders get ±1 rung",
                short_ladder)
        ctx.say("neighborhood: cube cap %d, pairs then triples then higher, duplicates removed", cube_cap)
        ctx.say("neighborhood: the known best configuration is an audit after the cube, not a climber proposal")
        neighbors: List[Tuple[str, object, Config]] = []
        seen = {seat.key}
        for name in ctx.space.names:
            if not ctx.ladders.movable(name):
                continue
            rungs = ctx.ladders.rungs.get(name) or []
            if len(rungs) <= 1:
                continue
            center = int(round(ctx.ladders.index(name, seat.params.get(name))))
            if len(rungs) <= short_ladder:
                indexes = [i for i in range(len(rungs)) if i != center]
            else:
                indexes = [i for i in (center - 1, center + 1) if 0 <= i < len(rungs)]
            for index in indexes:
                params = dict(seat.params)
                params[name] = rungs[index]
                config = ctx.space._assemble(params)
                if config is None or config.key in seen:
                    continue
                seen.add(config.key)
                neighbors.append((name, config.params.get(name), config))
        pending = [item for item in neighbors if ctx.score(item[2]) is None]
        ctx.say("neighborhood: %d neighbors, %d already scored", len(neighbors), len(neighbors) - len(pending))
        await self._send_and_wait([item[2] for item in pending])
        faces: List[Tuple[str, object, float]] = []
        for name in ctx.space.names:
            if not ctx.ladders.movable(name):
                continue
            points = [(ctx.ladders.index(name, seat.params.get(name)), incumbent, seat.params.get(name))]
            for axis, value, config in neighbors:
                if axis != name:
                    continue
                score = ctx.score(config)
                if score is None:
                    ctx.say("neighborhood %s %s unscored", name, value)
                    continue
                points.append((ctx.ladders.index(name, value), float(score), value))
                ctx.say("neighborhood %s %s %.4f (versus seat %+.4f)", name, value, score, score - incumbent)
            if len(points) <= 1 and not any(axis == name for axis, _value, _config in neighbors):
                continue
            line_max = max(point[1] for point in points)
            held = [point for point in points if point[1] >= line_max - margin]
            lo = max(0, int(math.floor(min(point[0] for point in held) + 1e-9)))
            hi = int(math.ceil(max(point[0] for point in held) - 1e-9))
            rungs = ctx.ladders.rungs.get(name) or []
            hi = min(len(rungs) - 1, hi) if rungs else hi
            ctx.say("neighborhood span %s line max %.4f (versus seat %+.4f), rungs %s .. %s",
                    name, line_max, line_max - incumbent,
                    rungs[lo] if rungs else lo, rungs[hi] if rungs else hi)
            others = [(value, score) for _index, score, value in points
                      if not self._same_value(value, seat.params.get(name))]
            if not others:
                continue
            best_value, best_score = max(others, key=lambda item: item[1])
            if best_score > incumbent:
                faces.append((name, best_value, best_score))
                ctx.say("neighborhood face %s %s %.4f beats the seat by %+.4f%s",
                        name, best_value, best_score, best_score - incumbent,
                        "" if best_score - incumbent >= margin else "; below 0.10, not a new region")
        ctx.say("neighborhood: %d winning faces", len(faces))
        cube: List[Config] = []
        cube_seen = {seat.key}
        for width in range(2, len(faces) + 1):
            for combo in itertools.combinations(faces, width):
                params = dict(seat.params)
                for name, value, _score in combo:
                    params[name] = value
                config = ctx.space._assemble(params)
                if config is None or config.key in cube_seen:
                    continue
                cube_seen.add(config.key)
                cube.append(config)
                if len(cube) >= cube_cap:
                    break
            if len(cube) >= cube_cap:
                break
        fresh_cube = [config for config in cube if ctx.score(config) is None]
        ctx.say("neighborhood cube: %d combinations within cap %d, %d to score",
                len(cube), cube_cap, len(fresh_cube))
        await self._send_and_wait(fresh_cube)
        climber_best = incumbent
        climber_key = seat.key
        for config in [item[2] for item in neighbors] + cube:
            score = ctx.score(config)
            if score is None:
                continue
            if config in cube:
                ctx.say("neighborhood cube %s %.4f (versus seat %+.4f)", config.key[:8], score, score - incumbent)
            if score > climber_best:
                climber_best = float(score)
                climber_key = config.key
        ctx.say("neighborhood climber best %s %.4f (versus seat %+.4f)",
                climber_key[:8], climber_best, climber_best - incumbent)
        offered = ctx.space._assemble(dict(REFERENCE_PEAK))
        if offered is not None:
            await self._send_and_wait([offered])
            audit = ctx.score(offered)
            if audit is None:
                ctx.say("neighborhood audit: known best unscored")
            else:
                ctx.say("neighborhood audit: known best %s %.4f (versus seat %+.4f); not a climber proposal",
                        offered.key[:8], audit, audit - incumbent)

    def _ingest_recent_scores(self) -> int:
        """Graft scores from the oracle run on 2026-10-14. Earlier rows in the same file are a different answer set."""
        import sqlite3
        from datetime import datetime, timezone
        from common.types import ScoreResult
        path = getattr(self.ctx.client.store, "path", None)
        if not path:
            return 0
        cut = datetime(2026, 10, 14, 1, 0, tzinfo=timezone.utc).timestamp()
        con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        rows = con.execute(
            "SELECT config, result FROM results WHERE created>=? AND score IS NOT NULL AND score>0", (cut,)).fetchall()
        con.close()
        client = self.ctx.client
        n = 0
        for raw_config, raw_result in rows:
            cfg = Config.from_dict(json.loads(raw_config))
            result = ScoreResult.from_dict(json.loads(raw_result))
            if result.score is None or (cfg.key, result.fidelity or FULL) in client.latest:
                continue
            client.configs[cfg.key] = cfg
            client.latest[(cfg.key, result.fidelity or FULL)] = result
            client._note_best(result)
            if result.gate_passed is not False:
                self.ctx.tied.note(cfg.key, cfg.params, result.score)
            n += 1
        return n

    def _sweep_indexes(self, name: str, peak_value) -> List[int]:
        """Every rung on a short ladder. Nine evenly spaced rungs, plus the peak, on a long one."""
        rungs = self.ctx.ladders.rungs.get(name) or []
        n = len(rungs)
        if n <= 1:
            return []
        if n <= 12:
            idxs = list(range(n))
        else:
            strata = 9
            idxs = [min(n - 1, int(round(i * (n - 1) / (strata - 1)))) for i in range(strata)]
        peak_i = int(round(self.ctx.ladders.index(name, peak_value)))
        idxs.append(max(0, min(n - 1, peak_i)))
        return sorted(set(idxs))

    async def _run_sweeps(self, margin: float) -> None:
        """One setting at a time across each ladder. The domain is the product of the spans that hold the maximum."""
        ctx = self.ctx
        scored = [(c, s) for c, s in ctx.scored()]
        if not scored:
            ctx.say("sweep: nothing scored")
            return
        best = max(s for _, s in scored)
        peak = max(scored, key=lambda t: t[1])[0]
        self._graft_peak = peak
        self._graft_best = best
        planned = []
        seen = {peak.key}
        for name in ctx.space.names:
            if not ctx.ladders.movable(name):
                continue
            rungs = ctx.ladders.rungs.get(name) or []
            for idx in self._sweep_indexes(name, peak.params.get(name)):
                params = dict(peak.params)
                params[name] = rungs[idx]
                config = ctx.space._assemble(params)
                if config is None or config.key in seen:
                    continue
                seen.add(config.key)
                planned.append((name, config.params.get(name), config))
        pending = [item for item in planned if ctx.score(item[2]) is None]
        ctx.say("sweep: peak %s %.4f, %d one-setting points, %d already scored",
                peak.key[:8], best, len(pending), len(planned) - len(pending))
        futs = [ctx.propose([config], "A", ["small"])[0] for _name, _value, config in pending]
        if futs:
            await ctx.wait_share(futs, 1.0, ctx.t0 + ctx.length)
        interval = {}
        offered = ctx.space._assemble(dict(REFERENCE_PEAK))
        for name in ctx.space.names:
            if not ctx.ladders.movable(name):
                continue
            rungs = ctx.ladders.rungs.get(name) or []
            points = []
            peak_value = peak.params.get(name)
            points.append((self.ctx.ladders.index(name, peak_value), best, peak_value))
            for axis, value, config in planned:
                if axis != name:
                    continue
                score = ctx.score(config)
                if score is None:
                    ctx.say("sweep %s %s unscored", name, value)
                    continue
                points.append((ctx.ladders.index(name, value), score, value))
                ctx.say("sweep %s %s %.4f (%+.4f)", name, value, score, score - best)
            if not points:
                continue
            axis_max = max(score for _i, score, _v in points)
            held = [p for p in points if p[1] >= axis_max - margin]
            lo = min(p[0] for p in held)
            hi = max(p[0] for p in held)
            lo_i = max(0, int(math.floor(lo + 1e-9)))
            hi_i = min(len(rungs) - 1, int(math.ceil(hi - 1e-9)))
            interval[name] = list(rungs[lo_i:hi_i + 1])
            holes = [p for p in points if lo - 1e-9 <= p[0] <= hi + 1e-9 and p[1] < axis_max - margin]
            ref = "" if offered is None else offered.params.get(name)
            ctx.say("sweep span %s %s .. %s (%d rungs, max %.4f, %d hold%s) reference %s",
                    name, interval[name][0], interval[name][-1], len(interval[name]), axis_max, len(held),
                    "" if not holes else f", {len(holes)} holes", ref)
        self._range = interval
        volume, n_open, n_axes = self._interval_size(interval)
        ref_in, misses = self._reference_fit(interval)
        ctx.say("sweep domain: size %.6g of the grid, %d of %d settings open, reference %s%s",
                volume, n_open, n_axes, "inside" if ref_in else "outside",
                "" if ref_in else " (" + ", ".join(misses) + ")")

    def _wave0_snapshot_path(self) -> str:
        explicit = os.environ.get("MINOS_WAVE0_SAVE")
        if explicit:
            return explicit
        tag = self.ctx.client.round.tag
        directory = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "test", "runs")
        return os.path.join(directory, f"wave0_{tag}.json")

    def _save_fixed_wave0(self) -> str:
        """The scores Wave 0 has when it stops waiting. A later test loads this and starts at Generation 1."""
        client = self.ctx.client
        rows = []
        for (key, fid), result in client.latest.items():
            if fid != FULL or result.score is None:
                continue
            cfg = client.configs.get(key)
            if cfg is None:
                continue
            rows.append({"config": cfg.to_dict(), "result": result.to_dict()})
        path = self._wave0_snapshot_path()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            json.dump({"round": client.round.tag, "n": len(rows), "rows": rows}, f)
        return path

    def _load_fixed_wave0(self, path: str) -> int:
        """Put a finished Wave 0 into the scored pool. Those configs are not sent again."""
        from common.types import ScoreResult
        blob = json.load(open(path))
        rows = blob["rows"] if isinstance(blob, dict) else blob
        client = self.ctx.client
        n = 0
        for row in rows:
            cfg = Config.from_dict(row["config"])
            result = ScoreResult.from_dict(row["result"])
            if result.score is None:
                continue
            client.configs[cfg.key] = cfg
            client.latest[(cfg.key, result.fidelity or FULL)] = result
            client._note_best(result)
            if result.gate_passed is not False:
                self.ctx.tied.note(cfg.key, cfg.params, result.score)
            n += 1
        return n

    async def run(self) -> List[Config]:
        ctx = self.ctx
        until = ctx.t0 + ctx.T(60 * float(tune("candidates.max_min")))
        from search.separate import enabled as separate_enabled, range_open
        if range_open():
            seeds = self.seeds()
            await self.wave(seeds, int(tune("candidates.riders")), "portfolio range", until)
            seen = {config.key for config in seeds}
            ctx._step_a_rest = [config for config in latin_hypercube(
                ctx.space, ctx.rng, max(0, ctx.cap - len(seeds))) if config.key not in seen]
            ctx._step_a_open = True
            from search.separate import note_range
            note_range(ctx)
            return await self.pick()
        fixed = os.environ.get("MINOS_WAVE0_FIXED")
        if fixed and os.environ.get("MINOS_NEIGHBORHOOD"):
            n = self._load_fixed_wave0(fixed)
            extra = self._ingest_recent_scores()
            ctx.say("neighborhood: %d scored configs loaded, %d recent scores added", n, extra)
            await self._run_neighborhood()
            return await self.pick()
        if fixed and os.environ.get("MINOS_RF_UCB"):
            from search.rf_ucb import propose as rf_propose
            n = self._load_fixed_wave0(fixed)
            extra = self._ingest_recent_scores()
            ctx.say("wave 0 fixed: %d scored configs loaded, %d graft scores added, wave 0 is not run", n, extra)
            batch = rf_propose(ctx, ctx.cap)
            if not batch:
                ctx.say("rf-ucb: no batch")
            else:
                ctx.say("rf-ucb batch runs, generation 1 is not bred")
                await self.wave(batch, 0, "rf-ucb", ctx.t0 + ctx.length)
            return await self.pick()
        if fixed and os.environ.get("MINOS_SWEEP"):
            n = self._load_fixed_wave0(fixed)
            extra = self._ingest_recent_scores()
            ctx.say("wave 0 fixed: %d scored configs loaded, %d graft scores added, wave 0 is not run", n, extra)
            await self._run_sweeps(0.10)
            return await self.pick()
        if fixed and os.environ.get("MINOS_GRAFT"):
            n = self._load_fixed_wave0(fixed)
            ctx.say("wave 0 fixed: %d scored configs loaded, wave 0 is not run", n)
            await self._run_grafts(0.10)
            return await self.pick()
        if fixed:
            n = self._load_fixed_wave0(fixed)
            ctx.say("wave 0 fixed: %d scored configs loaded, wave 0 is not run", n)
        else:
            wave0 = self.seeds()
            wave0 += [c for c in latin_hypercube(ctx.space, ctx.rng, max(0, ctx.cap - len(wave0)))
                      if all(c.key != x.key for x in wave0)]
            await self.wave(wave0, int(tune("candidates.riders")), "seeds + portfolio + hypercube", until,
                            settle=bool(int(tune("candidates.wave0_settle"))) or bool(os.environ.get("MINOS_FAST_HANDOFF")))
            path = self._save_fixed_wave0()
            ctx.say("wave 0 fixed for later tests: %s", path)
        if os.environ.get("MINOS_FAST_HANDOFF"):
            ctx.say("step A: the opening leaders go to step B; generation 1 is not run")
            return await self.pick()
        self._select_range()
        await self._open_pool_if_thin()
        if separate_enabled():
            from search.separate import opening_scouts, place_dummies
            await place_dummies(ctx)
            await opening_scouts(ctx)
            return await self.pick()
        # Generation 1 is waited out on its own clock. The wave-0 safety limit must not launch it and then
        # abandon it, which is what cut generation 2 on 20261008T175300Z.
        allowance = ctx.T(60 * float(tune("candidates.generation_min")))
        for gen in range(1, int(tune("candidates.generations")) + 1):
            if ctx.over() or ctx.clock() >= ctx.pick_at():
                break
            share = 0.0
            raw_share = os.environ.get("MINOS_RF_UCB")
            if raw_share:
                from search.rf_ucb import propose as rf_propose, share_of
                share = share_of(raw_share)
            n_rf = int(round(ctx.cap * share)) if gen == 1 else 0
            n_bias = int(round(ctx.cap * ELITE_BIAS_SHARE)) if gen == 1 and os.environ.get("MINOS_ELITE_BIAS") else 0
            if n_bias:
                ctx.say("elite bias: %.0f%% of generation 1 moves elite settings once the parent is within %.1f of the best; "
                        "the pair is active probability with the calling threshold; other settings stay in the genetic share",
                        100 * ELITE_BIAS_SHARE, ELITE_BAND)
            children = await self.breed(gen, max(0, ctx.cap - n_rf - n_bias))
            if n_bias:
                seen = {config.key for config in children}
                children += [config for config in self._elite_bias_children(n_bias) if config.key not in seen]
            if n_rf:
                seen = {config.key for config in children}
                children += [config for config in rf_propose(ctx, n_rf) if config.key not in seen]
            if not children:
                ctx.say("step A: generation %d has no children", gen)
                break
            ctx.say("step A: generation %d runs for %.0f min", gen, allowance / 60)
            await self.wave(children, 0, f"generation {gen}", ctx.clock() + allowance,
                            discovery=bool(os.environ.get("MINOS_DISCOVERY_HANDOFF")))
        if getattr(self, "_range", None):
            top = ctx.best()
            volume, n_open, n_axes = self._interval_size(self._range)
            inside = top is not None and self._inside(top[0])
            ref_in, misses = self._reference_fit(self._range)
            ctx.say(
                "step A metric: size %.6g of the grid, %d of %d settings open, scored maximum %s, reference %s%s",
                volume, n_open, n_axes, "inside" if inside else "outside",
                "inside" if ref_in else "outside",
                "" if ref_in else " (" + ", ".join(misses) + ")")
        return await self.pick()

    def wave_rows(self) -> List[dict]:
        return [{**{k: v for k, v in w.items() if k != "futs"}, "back": sum(1 for f in w["futs"] if f.done())}
                for w in self.waves]

    async def _apart(self, p: Point, taken: Sequence[Point], min_callset: float, min_settings: float) -> bool:
        """True when `p` is far enough from every config already taken. Yields during a long scan."""
        for t in taken:
            if p.config.key == t.config.key or not self.ctx.dist.far(p.config, t.config, min_callset, min_settings):
                return False
            # A long distance scan used to freeze the fleet and the dashboard. Yield every 50 ms of that
            # work. A fast scan does not yield, so a virtual round does not look stalled.
            await self._breathe()
        return True

    def _mutate(self, config: Config, axes: int, span: int) -> Optional[Config]:
        """Move `axes` settings by up to `span` rungs. Used to fill a thin pile and to look for a new fingerprint."""
        ctx = self.ctx
        params = dict(config.params)
        movable = [n for n in ctx.space.names if ctx.ladders.movable(n)]
        if not movable:
            return None
        for n in ctx.rng.sample(movable, min(axes, len(movable))):
            step = ctx.rng.randint(1, max(1, span)) * ctx.rng.choice((1, -1))
            v = ctx.ladders.shift(n, params.get(n), step)
            if v is None:
                v = ctx.ladders.shift(n, params.get(n), -step)
            if v is not None:
                params[n] = v
        return ctx.space._assemble(params)

    def _elite_bias_children(self, n: int) -> List[Config]:
        """Mutations of strong parents on the elite settings. The other settings are left to the genetic share."""
        ctx = self.ctx
        scored = [(config, score) for config, score in ctx.scored()]
        if n <= 0 or not scored:
            return []
        best = max(score for _config, score in scored)
        parents = [config for config, score in scored if score >= best - ELITE_BAND] or [scored[0][0]]
        axes = [name for name in ELITE_AXES if ctx.ladders.movable(name)]
        pair = [name for name in ELITE_PAIR if ctx.ladders.movable(name)]
        out, seen = [], set()
        tries = 0
        while len(out) < n and tries < n * 30:
            tries += 1
            params = dict(ctx.rng.choice(parents).params)
            if len(pair) == 2 and ctx.rng.random() < 0.5:
                moved = pair
            else:
                moved = ctx.rng.sample(axes, ctx.rng.randint(1, min(2, len(axes)))) if axes else []
            for name in moved:
                rungs = ctx.ladders.rungs.get(name) or []
                choices = [value for value in rungs if not self._same_value(value, params.get(name))]
                if choices:
                    params[name] = ctx.rng.choice(choices)
            child = ctx.space._assemble(params)
            if child is None or child.key in seen or not ctx.fresh(child):
                continue
            seen.add(child.key)
            out.append(child)
        ctx.say("elite bias: %d proposals from %d parents within %.1f of %.4f", len(out), len(parents), ELITE_BAND, best)
        return out

    async def _take(self, c: Optional[Config], seen: set, accept) -> bool:
        ctx = self.ctx
        if c is None or c.key in seen or not ctx.fresh(c):
            return False
        verdict = True if accept is None else accept(c)
        if hasattr(verdict, "__await__"):
            verdict = await verdict
        if not verdict:
            return False
        seen.add(c.key)
        return True

    async def _elite_probes(self, piles, n: int, *, why: str = "", accept=None) -> List[Config]:
        """While the GA owns the fleet, every child copies the best config and changes one to three settings."""
        from search.piles import mating_pool
        ctx = self.ctx
        top = ctx.best() if callable(getattr(ctx, "best", None)) else None
        elite = top[0] if top is not None else piles[0].config
        mates = []
        for pile in piles:
            mates.extend(p for p in mating_pool(ctx, pile) if p.key != elite.key)
        if not mates:
            mates = [elite]
        out, seen, tries = [], set(), 0
        names = list(ctx.space.names)
        while len(out) < n and tries < 40 * max(1, n):
            tries += 1
            child = self.elite_probe(elite, ctx.rng.choice(mates))
            if child is None:
                await self._breathe()
                continue
            diffs = sum(child.params.get(name) != elite.params.get(name) for name in names)
            if not 1 <= diffs <= 3:
                await self._breathe()
                continue
            if await self._take(child, seen, accept):
                out.append(child)
            else:
                await self._breathe()
        ctx.say("%s: %d elite probes of %s", why or "GA", len(out), elite.key[:8])
        return out

    async def breed_piles(self, gen: int, n: int, *, why: str = "", accept=None) -> List[Config]:
        """While the GA owns the fleet, every child is an elite probe of the best config: one to three settings.

        After that window the old split remains for a caller that still asks: half inside the piles, half
        six-setting mutations. The schedule does not ask. gen >= 3 breeds only inside the first pile.
        A pile with fewer than two parents mutates its best until 32 parents differ in at least four settings.
        """
        from search.climb import ga_hot
        from search.piles import mating_pool, refresh
        ctx = self.ctx
        piles = refresh(ctx)
        if not piles or n <= 0:
            return []
        if getattr(ctx, "ga_opening", False) or ga_hot(ctx):
            return await self._elite_probes(piles, n, why=why, accept=accept)
        if gen >= 3:
            piles = piles[:1]
            novelty_n = 0
        else:
            novelty_n = n // 2
        inside_n = n - novelty_n
        out, seen, tries = [], set(), 0
        while len(out) < inside_n and tries < 20 * max(1, inside_n):
            tries += 1
            pile = piles[len(out) % len(piles)]
            pool = mating_pool(ctx, pile)
            if len(pool) >= 2:
                a = ctx.rng.choice(pool)
                b = ctx.rng.choice([q for q in pool if q.key != a.key] or pool)
                c = self.child(a, b)
            else:
                c = self._mutate(pile.config, 4, 3)
            if await self._take(c, seen, accept):
                out.append(c)
                await self._breathe()
            else:
                await self._breathe()
        tries = 0
        while len(out) < n and tries < 20 * max(1, novelty_n):
            tries += 1
            base = ctx.rng.choice(piles).config
            c = self._mutate(base, 6, max(4, len(ctx.ladders.rungs.get("min_pruning", ())) // 4 or 4))
            if await self._take(c, seen, accept):
                out.append(c)
                await self._breathe()
            else:
                await self._breathe()
        ctx.say("%s: %d children across %d pile(s), %d looking for a new fingerprint",
                why or f"generation {gen}", len(out), len(piles), max(0, len(out) - inside_n))
        return out

    async def pick(self) -> List[Config]:
        from search.piles import enabled as piles_enabled
        if piles_enabled():
            return await self.pick_piles()
        ctx = self.ctx
        n = n_seats(ctx)
        pts = sorted(self.pool(), key=lambda p: -p.score)
        if not pts:
            ctx.say("step A: nothing scored; starting from the defaults")
            c = ctx.space._assemble({})
            return [c] if c is not None else []
        best = pts[0].score
        seated = [p for p in pts if best - p.score < SEAT_GAP]
        near = {p.config.key for p in seated}
        g = {name: [p for p in v if p.config.key in near] for name, v in (await self.groups()).items()}
        order = [seated[0]] + [p for i in range(max(len(v) for v in g.values()) if g else 0)
                               for name in ("fp", "core", "both") for p in g[name][i:i + 1]] + seated[1:]
        mc, ms = float(tune("candidates.min_distance")), float(tune("candidates.min_setting_distance"))
        names = list(ctx.space.names)
        near_n = n if n < 2 else n - 1
        taken: List[Point] = []
        for p in order:
            # Callset distance across the whole scored pool is minutes of work. Yield so the round loop
            # can poll the fleet and write the dashboard while the backup queue is built.
            if len(taken) >= near_n:
                break
            if await self._apart(p, taken, mc, ms):
                taken.append(p)
        far_seat = None
        rim_of = taken[0].config if taken else pts[0].config
        if n >= 2:
            for p in pts:
                if best - p.score > FAMILY_SEAT_FLOOR:
                    break
                if any(t.config.key == p.config.key for t in taken):
                    continue
                if not left_the_rim(p.config, rim_of):
                    continue
                if any(sum(p.config.params.get(nm) != t.config.params.get(nm) for nm in names) < ELITE_MIN_DIFFS
                       for t in taken):
                    continue
                if await self._apart(p, taken, mc, ms):
                    far_seat = p
                    break
            if far_seat is not None:
                taken.append(far_seat)
                ctx.say("step A: far seat %s %.4f (core %.4f, fp %.0f)", far_seat.config.key[:8],
                        far_seat.score, far_seat.core, far_seat.fp)
            elif len(taken) < n:
                for p in order:
                    if len(taken) >= n:
                        break
                    if await self._apart(p, taken, mc, ms):
                        taken.append(p)
        for p in order:
            if any(t.config.key == p.config.key for t in taken):
                continue
            if await self._apart(p, taken, mc, ms):
                taken.append(p)
        group_of = {t.config.key: next((name for name in ("both", "core", "fp")
                                        if any(q.config.key == t.config.key for q in g[name])), "score")
                    for t in taken}
        if far_seat is not None:
            group_of[far_seat.config.key] = "far"
        self.picked = [{"key": t.config.key, "score": round(t.score, 4), "core": round(t.core, 4), "fp": t.fp,
                        "group": group_of[t.config.key]} for t in taken[:n]]
        ctx.say("step A: %d candidate(s) climbing now, %d more when a hill is found, from %d scored configs: %s",
                min(n, len(taken)), max(0, len(taken) - n), len(pts),
                ", ".join(f"{t.config.key[:8]} {t.score:.4f} (core {t.core:.4f}, fp {t.fp:.0f})" for t in taken[:n]))
        return [t.config for t in taken]

    async def pick_piles(self) -> List[Config]:
        """One seat per callset pile, up to the climber cap. A pile with no fingerprint is not filled by a near copy."""
        from search.piles import refresh
        ctx = self.ctx
        n = n_seats(ctx)
        pts = list(self.pool())
        piles = refresh(ctx, [(p.config, p.score) for p in pts])
        ctx._piles = piles
        ctx._pile_seen = {member.key for pile in piles for member in pile.members}
        taken = piles[:n]
        self.picked = []
        for i, pile in enumerate(taken):
            point = next((p for p in pts if p.config.key == pile.config.key), None)
            self.picked.append({"key": pile.config.key, "score": round(pile.score, 4),
                                "core": round(point.core, 4) if point else 0.0,
                                "fp": point.fp if point else 0,
                                "group": pile.label(i)})
            ctx.say("step A: %s %s %.4f", pile.label(i), pile.config.key[:8], pile.score)
        ctx.say("step A: %d pile(s) climbing now, from %d scored configs", len(taken), len(pts))
        return [pile.config for pile in taken]


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
