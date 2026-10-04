"""RF-guided pattern search: forest surrogates rank axis moves (range edges included) by predicted gain per cost.

Differences from rf_ucb, each motivated by the real landscapes:
- utility is predicted IMPROVEMENT over the incumbent per unit cost, not absolute score per cost (absolute/cost
  makes any cheap variant win, so expensive axes are almost never moved);
- candidates are single-axis moves on a grid that includes both range edges, plus pair moves, random restarts,
  crossover, and joint moves of one expensive stage (reads, region, or assembly graph) with the other stages
  held still. Optima often sit at an edge (dangling 20, small padding). The forest still picks among those;
- a reserved share of full runs the forest cannot refuse. Two parents stay alive: the highest core, and the
  lowest false-positive rate among configs whose core is within a small band of it. Each reserved run copies
  one graph, region, or indel-prior setting from the low-false-positive parent onto the high-core parent.
  Read filters stay with the high-core parent, because a stricter read stage removes true indels along with
  the extra calls. Before those copies, the same share walks the calling threshold in the direction of the
  worst indel class: up when that class is mostly false positives, down when it is mostly misses. PCR model
  and read filters stay put. A flat hill raises the share from 20% to half;
- the live lineup scouts with that search, then splits the workers across promising regions. Each region keeps
  running the same search, and only on configs closer to its own center than to any other. Genotyping and
  PairHMM stay free inside every region. A center more than 1.5 below the best takes no workers, and a single
  center keeps one search on the whole space;
- warm start from previous rounds' configs as they were (cheap parameters included);
- parallel picks are kept apart (a candidate too close to a pending config is skipped);
- a cap on base-quality cutoffs: a cutoff above the BAM's highest base quality calls nothing (1000G rounds carry
  only 6/10/20/30, and on 20260928T092900Z such configs were 468 of the 469 zero scores).

Benchmark (2026-10-02): the round loop on a virtual clock against k-NN landscapes built from 8 truth-scored legacy
rounds (2.9k-117k configs each), GATK batching as LocalExecutor does it, 14 workers, 45 min, 3 seeds. Mean regret to
the round's best: this strategy 0.20 at 10 min and 0.105 final; rf_ucb + tpe 0.33 / 0.17 with the same warm start and
1.16 / 0.56 without one. Paired over 24 runs: 12 wins, 5 losses vs rf_ucb + tpe at the end (p=0.03), 21-2 vs the
lineup without a warm start. Most of the gain is the warm start; the rest is the gain-per-cost ranking and edge grid.
"""
from __future__ import annotations

import asyncio
import math
import random
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple, Union

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor

from common.params import ParamRanges
from common.scoring import DIFFICULTY_WEIGHTS, score_from_parts
from common.types import FULL, ClassCounts, Config
from search.base import Strategy
from search.client import EvalClient, Outcome
from search.family import designed_family
from search.schedule import run_streaming
from search.space import Space, _FAR_STARTS

BQ_PARAMS = ("min_base_quality_score", "base_quality_score_threshold")

# The 12 expensive settings are three GATK stages. A stage's knobs act on the same reads, the same
# region, or the same assembly graph, so a useful change usually moves more than one of them.
BLOCKS = (
    ("reads", ("min_base_quality_score", "min_mapping_quality_score", "max_reads_per_alignment_start")),
    ("region", ("active_probability_threshold", "min_assembly_region_size", "max_assembly_region_size",
                "assembly_region_padding")),
    ("graph", ("min_pruning", "min_dangling_branch_length", "recover_all_dangling_branches",
               "max_num_haplotypes_in_population")),
)
_COARSE = {
    "min_base_quality_score": (10, 15, 20, 25, 30, 40, 50),
    "min_mapping_quality_score": (0, 10, 20, 30, 40, 50, 60),
    "max_reads_per_alignment_start": (25, 50, 75, 100, 150, 200, 300),
    "active_probability_threshold": (0.001, 0.002, 0.005, 0.01, 0.02, 0.05),
    "min_assembly_region_size": (1, 25, 50, 100, 150, 200, 300),
    "max_assembly_region_size": (100, 150, 200, 300, 400, 500, 700),
    "assembly_region_padding": (0, 50, 100, 200, 300, 500),
    "min_pruning": (2, 3, 4, 5, 6, 8, 10),
    "min_dangling_branch_length": (2, 4, 8, 12, 20),
    "recover_all_dangling_branches": (False, True),
    "max_num_haplotypes_in_population": (8, 16, 32, 64, 128),
}
# Calling-confidence ladder. The live hill sits near 80, so the steps above it are the ones that
# can drop a false large indel. Steps below it are only used when the weak class is misses.
_CALL_LADDER = (30.0, 45.0, 60.0, 70.0, 80.0, 90.0, 100.0)
# One of these, copied from a nearby low-false-positive config onto the highest-core config.
# Reads are absent on purpose: stricter base quality and mapping quality take true indels with the extras.
_CORE_STEPS = (
    "min_pruning", "min_dangling_branch_length", "recover_all_dangling_branches",
    "max_num_haplotypes_in_population", "active_probability_threshold", "min_assembly_region_size",
    "max_assembly_region_size", "assembly_region_padding", "pair_hmm_gap_continuation_penalty",
    "indel_heterozygosity",
)


def bam_max_quality(bam: str, region: str, reads: int = 20000) -> Optional[int]:
    """Highest base quality on the first `reads` reads of the region, or None when the BAM cannot be read."""
    try:
        import pysam
        top = 0
        with pysam.AlignmentFile(bam) as f:
            for i, read in enumerate(f.fetch(region=region)):
                if i >= reads:
                    break
                q = read.query_qualities
                if q is not None and len(q):
                    top = max(top, max(q))
        return top or None
    except Exception:  # noqa: BLE001  (no BAM, no index, no pysam: search without the cap)
        return None


def _coarse(space: Space, name: str, cap_at: Optional[int]) -> List[object]:
    spec = space.specs[name]
    out = []
    for v in _COARSE.get(name, ()):
        if spec.type in ("int", "float") and spec.min is not None and float(v) < float(spec.min):
            continue
        if spec.type in ("int", "float") and spec.max is not None and float(v) > float(spec.max):
            continue
        if cap_at is not None and name in BQ_PARAMS and v > cap_at:
            continue
        out.append(v)
    return out


def _over_cap(c: Config, cap_at: Optional[int]) -> bool:
    return cap_at is not None and any(c.params.get(name, 0) > cap_at for name in BQ_PARAMS if name in c.params)


def _block_names(space: Space, block: Sequence[str]) -> List[str]:
    return [n for n in block if n in space.names]


def block_moves(space: Space, rng: random.Random, base: Config, n_per_block: int,
                cap_at: Optional[int] = None) -> List[Config]:
    """Joint coarse moves of one expensive stage. The other stages stay at `base`.

    Each child changes at least two settings of one stage. Nothing here is scored.
    """
    if n_per_block <= 0:
        return []
    out: List[Config] = []
    seen = {base.key}
    for _, block in BLOCKS:
        names = _block_names(space, block)
        if len(names) < 2:
            continue
        proposals: List[Dict[str, object]] = []
        for end in (0, -1):
            proposals.append({n: _coarse(space, n, cap_at)[end] for n in names if _coarse(space, n, cap_at)})
        for _ in range(n_per_block):
            proposals.append({n: rng.choice(_coarse(space, n, cap_at)) for n in names if _coarse(space, n, cap_at)})
        for chosen in proposals:
            params = dict(base.params)
            params.update(chosen)
            n_diff = sum(1 for n in names if not space._same(space.specs[n], params[n], base.params.get(n)))
            if n_diff < 2:
                continue
            c = space._assemble(params)
            if c is None or c.key in seen or _over_cap(c, cap_at):
                continue
            outside = [n for n in space.expensive if n not in names]
            if any(not space._same(space.specs[n], c.params.get(n), base.params.get(n)) for n in outside):
                continue
            seen.add(c.key)
            out.append(c)
    return out


def block_swaps(space: Space, rng: random.Random, parents: Sequence[Config], n: int,
                cap_at: Optional[int] = None) -> List[Config]:
    """Give one parent a whole expensive stage taken from another scored config."""
    if n <= 0 or len(parents) < 2:
        return []
    out: List[Config] = []
    seen = {c.key for c in parents}
    blocks = [b for _, b in BLOCKS]
    for _ in range(n):
        a, b = rng.sample(list(parents), 2)
        names = _block_names(space, rng.choice(blocks))
        if not any(not space._same(space.specs[n], a.params.get(n), b.params.get(n)) for n in names):
            continue
        params = dict(a.params)
        for name in names:
            params[name] = b.params[name]
        c = space._assemble(params)
        if c is None or c.key in seen or _over_cap(c, cap_at):
            continue
        seen.add(c.key)
        out.append(c)
    return out


def core_steps(space: Space, high_core: Config, low_fp: Config, cap_at: Optional[int] = None) -> List[Config]:
    """Copy one graph, region, or indel-prior setting from `low_fp` onto `high_core`.

    Read filters stay with `high_core`. Each child changes exactly one of those settings. Nothing here is scored.
    """
    names = [n for n in _CORE_STEPS if n in space.names]
    out: List[Config] = []
    seen = {high_core.key, low_fp.key}
    for name in names:
        spec = space.specs[name]
        if space._same(spec, high_core.params.get(name), low_fp.params.get(name)):
            continue
        params = dict(high_core.params)
        params[name] = low_fp.params[name]
        c = space._assemble(params)
        if c is None or c.key in seen or _over_cap(c, cap_at):
            continue
        changed = [n for n in space.names
                   if not space._same(space.specs[n], c.params.get(n), high_core.params.get(n))]
        if changed != [name]:
            continue
        seen.add(c.key)
        out.append(c)
    return out


def _class_f1(counts: ClassCounts) -> Tuple[float, float, float]:
    precision = counts.tp / (counts.tp + counts.fp) if counts.tp + counts.fp else 0.0
    recall = counts.tp / (counts.tp + counts.fn) if counts.tp + counts.fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return f1, precision, recall


def worst_class(counts: Dict[str, ClassCounts]) -> Optional[Tuple[float, str, float, float]]:
    """The difficulty class with the largest weight x (1 - F1): (loss, name, precision, recall)."""
    worst: Optional[Tuple[float, str, float, float]] = None
    for name, weight in DIFFICULTY_WEIGHTS.items():
        row = counts.get(name)
        if row is None:
            continue
        f1, precision, recall = _class_f1(row)
        loss = float(weight) * (1.0 - f1)
        if worst is None or loss > worst[0]:
            worst = (loss, name, precision, recall)
    return worst


def confidence_steps(space: Space, base: Config, counts: Dict[str, ClassCounts],
                     cap_at: Optional[int] = None) -> List[Config]:
    """Move only the calling threshold, in the direction the worst indel class needs.

    Precision below recall means the class is overcalled, so the threshold goes up. Recall below precision
    means it is missed, so the threshold goes down. A SNP hole returns nothing: those F1s are already high
    and a threshold move would hit the indels too. Reads and the PCR model are left as they are.
    """
    name = "standard_min_confidence_threshold_for_calling"
    if name not in space.names or not counts:
        return []
    hole = worst_class(counts)
    if hole is None or not hole[1].startswith("indel"):
        return []
    _loss, _klass, precision, recall = hole
    spec = space.specs[name]
    current = float(base.params.get(name, spec.default))
    if precision < recall:
        wanted = [v for v in _CALL_LADDER if v > current + 1e-9]
    elif recall < precision:
        wanted = [v for v in _CALL_LADDER if v < current - 1e-9]
    else:
        return []
    out: List[Config] = []
    seen = {base.key}
    for value in wanted:
        if spec.min is not None and value < float(spec.min) - 1e-9:
            continue
        if spec.max is not None and value > float(spec.max) + 1e-9:
            continue
        params = dict(base.params)
        params[name] = float(value)
        c = space._assemble(params)
        if c is None or c.key in seen or _over_cap(c, cap_at):
            continue
        changed = [n for n in space.names
                   if not space._same(space.specs[n], c.params.get(n), base.params.get(n))]
        if changed != [name]:
            continue
        seen.add(c.key)
        out.append(c)
    return out


def cross_parents(space: Space, rng: random.Random, parents: Sequence[Config], n: int,
                  cap_at: Optional[int] = None) -> List[Config]:
    """Children that mix the expensive settings of two scored configs.

    Cheap settings stay with the first parent, so a child is one full-run key and the designed family can still
    retune genotyping and PairHMM. A child is dropped when a base-quality cutoff is above `cap_at`. Nothing here
    is scored: the caller ranks the children and runs one.
    """
    if n <= 0 or len(parents) < 2:
        return []
    out: List[Config] = []
    seen = {c.key for c in parents}
    for _ in range(n):
        a, b = rng.sample(list(parents), 2)
        diffs = [name for name in space.expensive
                 if not space._same(space.specs[name], a.params.get(name), b.params.get(name))]
        if not diffs:
            continue
        params = dict(a.params)
        taken = [name for name in diffs if rng.random() < 0.5]
        if not taken:
            taken = [rng.choice(diffs)]
        for name in taken:
            params[name] = b.params[name]
        c = space._assemble(params)
        if c is None or c.key in seen:
            continue
        if cap_at is not None and any(c.params.get(name, 0) > cap_at for name in BQ_PARAMS if name in c.params):
            continue
        seen.add(c.key)
        out.append(c)
    return out


def expensive_distance(space: Space, a: Config, b: Config) -> float:
    """Distance on the expensive settings only, each axis scaled to [0, 1].

    Genotyping and PairHMM stay out of it. A region can still move the gap penalty and the base-quality
    threshold, which is where the score actually changes, while two cores do not climb the same read filter.
    """
    total = 0.0
    for name in space.expensive:
        spec = space.specs[name]
        ua = space._to_unit(spec, a.params.get(name, spec.default))
        ub = space._to_unit(spec, b.params.get(name, spec.default))
        total += (ua - ub) ** 2
    return math.sqrt(total)


def region_centers(space: Space, scored: Sequence[Tuple[float, Config]], k: int, separation: float,
                   band: float = 1.5) -> List[Config]:
    """Up to `k` configs within `band` of the best, each at least `separation` apart on the expensive axes.

    A config further than `band` below the best is not a promising region, so it does not take cores.
    """
    if not scored or k <= 0:
        return []
    best = scored[0][0]
    centers: List[Config] = []
    for score, config in scored:
        if score < best - band:
            continue
        if all(expensive_distance(space, config, other) >= separation for other in centers):
            centers.append(config)
        if len(centers) >= k:
            break
    return centers


def pattern(ranges: ParamRanges, *, tool: str = "gatk", in_flight: int = 8, family: int = 6,
            family_max_cost: float = 0.15, portfolio: int = 8, n_trees: int = 40, kappa: float = 0.5,
            kappa_end: float = 0.05, kappa_horizon: int = 150, grid: int = 7, n_pairs: int = 300,
            n_cross: int = 80, n_block: int = 8, n_block_swap: int = 24, n_random: int = 60,
            anchors: int = 3, parents: int = 8, min_sep: float = 0.03,
            late_s: float = 180.0, max_bq: Union[int, str, None] = "auto", edge_probes: bool = True,
            max_pending: int = 256, mix_share: float = 0.20, core_band: float = 0.005,
            seed: int = 0, tag: str = "pattern",
            accept: Optional[Callable[[Config], bool]] = None, slots: Optional[int] = None,
            train_inside: bool = False, anchor: Optional[Config] = None,
            stop_when: Optional[Callable[[], bool]] = None, far_starts: bool = False) -> Strategy:
    """`max_bq`: highest useful base-quality cutoff; "auto" reads it from the round BAM, None disables the cap.

    `accept` drops any config outside this run's region. `slots` is how many of this run's proposals may be in
    flight; unset, the run takes every free worker. `train_inside` fits the forest on the region only.
    `anchor` is the region's center. `stop_when` ends proposing (the scout). `far_starts` adds the two starts
    outside the default neighborhood.
    """
    space = Space(ranges, tool)
    rng = random.Random(seed)
    cap: List[Optional[int]] = [max_bq if isinstance(max_bq, int) else None]

    def capped(c: Optional[Config]) -> Optional[Config]:
        if c is None:
            return None
        if cap[0] is not None and any(c.params.get(n, space.specs[n].default) > cap[0]
                                      for n in BQ_PARAMS if n in space.specs):
            return None
        if accept is not None and not accept(c):
            return None
        return c

    def axis_values(name: str, current) -> List[object]:
        s = space.specs[name]
        if s.allowed:
            return list(s.allowed)
        if s.type == "bool":
            return [False, True]
        u = space._to_unit(s, current)
        units = {i / (grid - 1) for i in range(grid)} | {u - 0.15, u - 0.05, u + 0.05, u + 0.15}
        vals = []
        for x in sorted(units):
            if 0.0 <= x <= 1.0:
                v = space._from_unit(s, x)
                if cap[0] is not None and name in BQ_PARAMS:
                    v = min(v, cap[0])
                vals.append(v)
        return vals

    async def run(client: EvalClient) -> None:
        if max_bq == "auto":
            cap[0] = bam_max_quality(client.round.bam, client.round.region)
        failed: Set[str] = set()
        proposed: Set[str] = set()
        pending_x: Dict[str, np.ndarray] = {}
        models: List[Optional[Tuple[ExtraTreesRegressor, ExtraTreesRegressor]]] = [None]
        fitted_on = [0]
        full_sent = [0]
        reserved_sent = [0]
        step_queue: List[Config] = []
        step_key: List[Optional[Tuple]] = [None]
        wave_left: Dict[int, int] = {}
        wave_of: Dict[int, int] = {}
        wave_seq = [0]

        def rows():
            good, bad = [], set(failed)
            for (key, fid), result in client.latest.items():
                if fid != FULL or key not in client.configs or client.configs[key].tool != tool:
                    continue
                if result.score is None or result.gate_passed is False:
                    bad.add(key)
                    continue
                config = client.configs[key]
                if train_inside and accept is not None and not accept(config):
                    continue
                core = result.core if result.core is not None else result.score / 100.0
                fp = (result.region_fp / result.n_targets) if result.n_targets else 0.0
                good.append((config, float(core), float(fp), float(result.score), result.counts))
            return good, bad

        def refit():
            good, bad = rows()
            if len(good) < 4:
                return
            n = len(good) + len(bad)
            if models[0] is not None and n - fitted_on[0] < max(1, n // 20):
                return
            worst_core = min(g[1] for g in good)
            worst_fp = max(g[2] for g in good)
            configs = [g[0] for g in good] + [client.configs[k] for k in bad if k in client.configs]
            cores = [g[1] for g in good] + [worst_core] * (len(configs) - len(good))
            fps = [g[2] for g in good] + [worst_fp] * (len(configs) - len(good))
            x = np.array([space.encode(c) for c in configs])
            st = rng.randrange(1 << 30)
            models[0] = (
                ExtraTreesRegressor(n_estimators=n_trees, min_samples_leaf=2, max_features=0.5,
                                    random_state=st).fit(x, np.array(cores)),
                ExtraTreesRegressor(n_estimators=n_trees, min_samples_leaf=2, max_features=0.5,
                                    random_state=st + 1).fit(x, np.array(fps)))
            fitted_on[0] = n

        def ucb(x: np.ndarray, k: float) -> np.ndarray:
            core_m, fp_m = models[0]
            core_t = np.stack([t.predict(x) for t in core_m.estimators_])
            fp_t = np.stack([t.predict(x) for t in fp_m.estimators_])
            per_tree = np.vectorize(score_from_parts)(core_t, fp_t)
            return per_tree.mean(axis=0) + k * per_tree.std(axis=0)

        def fresh(c: Optional[Config]) -> bool:
            if c is None or (accept is not None and not accept(c)):
                return False
            return c.key not in proposed and (c.key, FULL) not in client.latest and c.key not in failed

        def late() -> bool:
            left = client.time_left()
            return left is not None and left < late_s

        def moves_from(base: Config) -> List[Config]:
            out = []
            for name in space.expensive:
                cur = base.params.get(name, space.specs[name].default)
                for v in axis_values(name, cur):
                    if v == cur:
                        continue
                    c = capped(space._assemble({**base.params, name: v}))
                    if fresh(c):
                        out.append(c)
            return out

        def candidates(ranked: List[Config]) -> List[Config]:
            cands: Dict[str, Config] = {}
            singles: List[Config] = []
            for i, a in enumerate(ranked[:anchors]):
                m = moves_from(a)
                if i == 0:
                    singles = m
                for c in m:
                    cands[c.key] = c
            inc = ranked[0]
            for _ in range(n_pairs if len(singles) > 1 else 0):
                a, b = rng.sample(singles, 2)
                diff_a = {k: v for k, v in a.params.items() if inc.params.get(k) != v}
                c = capped(space._assemble({**b.params, **diff_a}))
                if fresh(c):
                    cands[c.key] = c
            for c in cross_parents(space, rng, ranked[:parents], n_cross, cap[0]):
                if fresh(c):
                    cands[c.key] = c
            for a in ranked[:anchors]:
                for c in block_moves(space, rng, a, n_block, cap[0]):
                    if fresh(c):
                        cands[c.key] = c
            for c in block_swaps(space, rng, ranked[:parents], n_block_swap, cap[0]):
                if fresh(c):
                    cands[c.key] = c
            for _ in range(n_random):
                c = capped(space.random(rng))
                if fresh(c):
                    cands[c.key] = c
            return list(cands.values())

        def cold(base: Config) -> List[Config]:
            moves = moves_from(base)
            if not moves:
                return []
            head = moves[0]
            fam = designed_family(space, head, family,
                                  lambda c: fresh(c) and c.key != head.key
                                  and client.variant_cost(head, c) <= family_max_cost)
            return [head] + fam

        def pick() -> List[Config]:
            good, _ = rows()
            if not good or models[0] is None:
                base = good[0][0] if good else anchor
                return [] if base is None else cold(base)
            ranked = [g[0] for g in sorted(good, key=lambda t: -t[3])]
            inc = ranked[0]
            if late():
                return designed_family(space, inc, family,
                                       lambda c: fresh(c) and client.variant_cost(inc, c) <= 0.05, genotyping_only=True)
            n_obs = len(good)
            k = kappa + (kappa_end - kappa) * min(1.0, n_obs / max(kappa_horizon, 1))
            best_seen = max(g[3] for g in good)
            cands = candidates(ranked)
            cheap = designed_family(space, inc, family,
                                    lambda c: fresh(c) and client.variant_cost(inc, c) <= family_max_cost)
            pool = cands + cheap
            if not pool:
                return []
            x = np.array([space.encode(c) for c in pool])
            gain = ucb(x, k) - best_seen
            cost = np.array([max(client.variant_cost(inc, c), 0.02) if 0 < client.variant_cost(inc, c) < 1.0
                             else 1.0 for c in pool])
            util = np.where(gain > 0, gain / cost, gain)          # no predicted gain: rank by least loss, cost-blind
            if pending_x:
                px = np.stack(list(pending_x.values()))
                d = np.min(np.abs(x[:, None, :] - px[None, :, :]).max(axis=2), axis=1)
                util = np.where(d < min_sep, -np.inf, util)
            order = np.argsort(-util)
            best_i = int(order[0])
            if not np.isfinite(util[best_i]):
                return []
            chosen = pool[best_i]
            if best_i >= len(cands):        # a cheap variant of the incumbent won: send the whole designed family
                return [c for c in cheap]
            fam = designed_family(space, chosen, family,
                                  lambda c: fresh(c) and c.key != chosen.key
                                  and client.variant_cost(chosen, c) <= family_max_cost)
            return [chosen] + fam

        def flat() -> bool:
            good, _ = rows()
            if len(good) < 8:
                return False
            scores = sorted((g[3] for g in good), reverse=True)
            return scores[0] - scores[3] < 0.15

        def share() -> float:
            return 0.5 if flat() else mix_share

        def reserved_pair() -> Optional[Tuple[Config, Config]]:
            """Highest core, and the lowest false-positive rate still inside the core band that can donate a step."""
            good, _ = rows()
            if len(good) < 4:
                return None
            high = max(good, key=lambda t: (t[1], t[3]))
            floor = high[1] - core_band
            near = sorted((g for g in good if g[0].key != high[0].key and g[1] >= floor), key=lambda t: (t[2], -t[3]))
            for low in near:
                if core_steps(space, high[0], low[0], cap[0]):
                    return high[0], low[0]
            return None

        def reserved_batch() -> List[Config]:
            good, _ = rows()
            found = reserved_pair()
            if not good and found is None:
                return []
            best = max(good, key=lambda t: t[3]) if good else None
            high = found[0] if found else (best[0] if best else None)
            low = found[1] if found else None
            if high is None:
                return []
            key = (high.key, None if low is None else low.key, None if best is None else best[0].key)
            if step_key[0] != key:
                counts = best[4] if best is not None else {}
                step_queue[:] = confidence_steps(space, high if best is None else best[0], counts, cap[0])
                if low is not None:
                    step_queue.extend(core_steps(space, high, low, cap[0]))
                step_key[0] = key
            while step_queue and not fresh(step_queue[0]):
                step_queue.pop(0)
            if not step_queue:
                return []
            head = step_queue.pop(0)
            fam = designed_family(space, head, family,
                                  lambda c: fresh(c) and c.key != head.key
                                  and client.variant_cost(head, c) <= family_max_cost)
            return [head] + fam

        def note(fut):
            wid = wave_of.pop(id(fut), None)
            if wid is not None:
                wave_left[wid] -= 1
                if wave_left[wid] <= 0:
                    wave_left.pop(wid, None)
            o: Outcome = fut.result()
            if o.config is not None:
                pending_x.pop(o.config.key, None)
            if not o.ok and o.config is not None and not (o.error or "").startswith("budget"):
                failed.add(o.config.key)

        def room(pending) -> bool:
            if slots is not None:
                return len(wave_left) < slots
            limit = max(max_pending, (family + 1) * (client.capacity() or 0))
            if len(pending) >= limit:
                return False
            idle = client.idle_slots()
            return (len(pending) < in_flight) if idle is None else (idle > 0 or not pending)

        def occupy(batch: List[Config], futs: List) -> List:
            """One proposed batch is one GATK run, so it takes one core, however many cheap variants it holds."""
            if not batch:
                return []
            wid = wave_seq[0]
            wave_seq[0] += 1
            wave_left[wid] = len(futs)
            for fut in futs:
                wave_of[id(fut)] = wid
            return futs

        def send(batch: List[Config]):
            if slots is not None and len(wave_left) >= slots:
                return []
            proposed.update(c.key for c in batch)
            if batch:
                full_sent[0] += 1
                pending_x[batch[0].key] = np.array(space.encode(batch[0]))
            return occupy(batch, list(client.propose_many(batch, tag=tag)))

        def fill():
            if stop_when is not None and stop_when():
                return []
            refit()
            batch: List[Config] = []
            if not late() and reserved_sent[0] < share() * (full_sent[0] + 1):
                batch = reserved_batch()
                if batch:
                    reserved_sent[0] += 1
                    return send(batch)
            batch = pick()
            if not batch and accept is not None and not late():
                for _ in range(80):
                    c = capped(space.random(rng))
                    if fresh(c):
                        batch = [c]
                        break
            return send(batch)

        # warm start: the region center, the GATK default, previous rounds' configs as they were, and with this
        # round's cheap retune. A region keeps only the configs that fall inside it.
        seeds: List[Config] = []
        if anchor is not None:
            a = capped(anchor)
            if a is not None:
                seeds.append(a)
        d = capped(space._assemble({}))
        if d is not None:
            seeds.append(d)
        for raw in client.store.portfolio(limit=portfolio, tool=tool):
            for c in (capped(space._assemble(dict(raw.params))), capped(space.with_cheap_defaults(raw))):
                if c is not None:
                    seeds.append(c)
        n_first = len(seeds)
        # one edge probe per expensive axis from the default, so the first model sees both ends of every range
        if d is not None and edge_probes:
            for name in space.expensive:
                s = space.specs[name]
                if s.allowed or s.type == "bool":
                    continue
                cur = d.params.get(name, s.default)
                far = 1.0 if space._to_unit(s, cur) < 0.5 else 0.0
                c = capped(space._assemble({**d.params, name: space._from_unit(s, far)}))
                if c is not None:
                    seeds.append(c)
        if far_starts:
            for params in _FAR_STARTS:
                c = capped(space._assemble(dict(params)))
                if c is not None:
                    seeds.append(c)
        initial, probes = [], []
        seen = set()
        for i, c in enumerate(seeds):
            if fresh(c) and c.key not in seen:
                seen.add(c.key)
                (initial if i < n_first else probes).append(c)
        futs = []

        def launch(batch: List[Config], priority: int) -> None:
            if not batch or (slots is not None and len(wave_left) >= slots):
                return
            proposed.update(c.key for c in batch)
            futs.extend(occupy(batch, list(client.propose_many(batch, priority=priority, tag=tag))))

        if not late():
            # default and previous rounds first, each with its cheap family: the first scores arrive in one run
            for c in initial:
                if slots is not None and len(wave_left) >= slots:
                    break
                fam = designed_family(space, c, family, lambda v: fresh(v) and v.key != c.key and v.key not in seen
                                      and client.variant_cost(c, v) <= family_max_cost)
                seen.update(v.key for v in fam)
                launch([c] + fam, 1)
            if edge_probes or far_starts:
                launch(probes, 0)
        await run_streaming(client, room, fill, note, futs)

    return run


def pattern_regions(ranges: ParamRanges, *, tool: str = "gatk", in_flight: int = 8, n_regions: int = 0,
                    separation: float = 0.35, scout_s: float = 480.0, scout_min: int = 0, band: float = 1.5,
                    seed: int = 0, **kwargs) -> Strategy:
    """Scout with the pattern search, then give each promising region its own cores and its own pattern search.

    A region is the set of configs closer to its center than to any other center, measured on the expensive
    settings only. Genotyping and PairHMM stay free inside every region, so a core can still move the knobs
    that change the score. The centers stay put for the rest of the round: a core does not leave its region.
    One center, or none far enough apart, keeps a single pattern search on the whole space.
    """
    space = Space(ranges, tool)
    k = n_regions if n_regions > 0 else max(2, min(4, max(1, in_flight) // 3))
    need = scout_min if scout_min > 0 else max(16, k * 6)

    async def run(client: EvalClient) -> None:
        run.centers = []  # type: ignore[attr-defined]
        started = client.time_left()

        def scored() -> List[Tuple[float, Config]]:
            out = []
            for (key, fid), result in client.latest.items():
                if fid != FULL or key not in client.configs or client.configs[key].tool != tool:
                    continue
                if result.score is None or result.gate_passed is False:
                    continue
                out.append((float(result.score), client.configs[key]))
            out.sort(key=lambda item: -item[0])
            return out

        def scout_done() -> bool:
            centers = region_centers(space, scored(), k, separation, band)
            if len(scored()) >= need and len(centers) >= min(k, 2):
                return True
            if started is None or client.time_left() is None:
                return False
            return started - client.time_left() >= scout_s

        scout = pattern(ranges, tool=tool, in_flight=in_flight, seed=seed, tag="scout",
                        stop_when=scout_done, far_starts=True, **kwargs)
        await scout(client)
        if client.exhausted:
            return

        centers = region_centers(space, scored(), k, separation, band)
        run.centers = centers  # type: ignore[attr-defined]
        if len(centers) < 2:
            print(f"regions 1 slots {in_flight}", flush=True)
            rest = pattern(ranges, tool=tool, in_flight=in_flight, seed=seed + 1, tag="pattern",
                           edge_probes=False, **kwargs)
            await rest(client)
            return

        base, extra = divmod(in_flight, len(centers))
        slot_of = [base + (1 if i < extra else 0) for i in range(len(centers))]

        def owner(config: Config) -> int:
            distances = [expensive_distance(space, config, center) for center in centers]
            return min(range(len(distances)), key=lambda i: (distances[i], i))

        print(f"regions {len(centers)} slots {' '.join(str(s) for s in slot_of)}", flush=True)
        lanes = []
        for i, center in enumerate(centers):
            if slot_of[i] <= 0:
                continue
            lanes.append(pattern(ranges, tool=tool, in_flight=slot_of[i], slots=slot_of[i], seed=seed + 1 + i,
                                 tag=f"region{i}", accept=lambda config, i=i: owner(config) == i,
                                 train_inside=True, anchor=center, edge_probes=False, **kwargs))
        await asyncio.gather(*(lane(client) for lane in lanes))

    return run
