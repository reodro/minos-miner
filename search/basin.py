"""Basin racing for a large fleet.

Local score decides which structural neighborhood is worth a look. It does not pick a winner inside that
neighborhood. A submission changes only when a different structural signature is ahead by the band. Official
scores are not available during the round; a known official config is only the first seed.

One family is one HaplotypeCaller traversal. The other members change only genotyping parameters, so they
re-genotype in about a second and a half instead of starting another traversal.
"""
from __future__ import annotations

import math
import random
from typing import Dict, List, Optional, Sequence, Set, Tuple

from common.params import ParamRanges
from common.types import FULL, Config
from search.base import Strategy
from search.client import EvalClient
from search.schedule import run_streaming
from search.space import (CALL_CONF_GRID, CONTAMINATION_GRID, GENOTYPING_PARAMS, PINNED, Space)

# Measured on 20260928T092900Z: local 87.4 and local 85.1 were different official bands (~88.16 vs ~86.04).
# Inside 87.38–87.41 the official scores did not separate. A gain below this does not change the submission.
BAND = 0.15
# Unit-space bucket. Padding 25 and 100 land apart; a one-step confidence change does not, because
# genotyping parameters are not part of the signature.
BIN = 0.25
FINISH_S = 15 * 60
# Short intensification. Longer runs inside one neighborhood did not move the official score.
RACE_PER_BASIN = 24
RACE_BASINS = 8

# Official 88.163 on round 20260928T092900Z. A seed for the next round, not a claim about that round's basin.
OFFICIAL_SEED = {
    "min_base_quality_score": 16,
    "min_mapping_quality_score": 46,
    "base_quality_score_threshold": 22,
    "standard_min_confidence_threshold_for_calling": 80.0,
    "pcr_indel_model": "NONE",
    "min_pruning": 2,
    "max_alternate_alleles": 6,
    "min_dangling_branch_length": 20,
    "recover_all_dangling_branches": False,
    "max_num_haplotypes_in_population": 128,
    "active_probability_threshold": 0.002,
    "min_assembly_region_size": 50,
    "max_assembly_region_size": 270,
    "assembly_region_padding": 100,
    "pair_hmm_gap_continuation_penalty": 7,
    "phred_scaled_global_read_mismapping_rate": 45,
    "heterozygosity": 0.002154,
    "indel_heterozygosity": 0.000125,
    "contamination_fraction_to_filter": 0.04,
    "max_reads_per_alignment_start": 116,
}


def structural_names(space: Space) -> List[str]:
    """Parameters that change the traversal. Genotyping is cheap and stays out of the basin id."""
    skip = set(GENOTYPING_PARAMS)
    return [n for n in space.names if n not in skip]


def signature(space: Space, config: Config) -> Tuple:
    out = []
    for name in structural_names(space):
        spec = space.specs[name]
        value = config.params.get(name, spec.default)
        if spec.type == "bool" or spec.allowed:
            out.append(value)
        else:
            out.append(int(space._to_unit(spec, value) / BIN))
    return tuple(out)


def struct_distance(space: Space, a: Config, b: Config) -> float:
    names = structural_names(space)
    if not names:
        return 0.0
    total = 0.0
    for name in names:
        spec = space.specs[name]
        ua = space._to_unit(spec, a.params.get(name, spec.default))
        ub = space._to_unit(spec, b.params.get(name, spec.default))
        total += (ua - ub) ** 2
    return math.sqrt(total / len(names))


def _clamp_num(spec, value):
    if spec.type == "int":
        return int(min(spec.max, max(spec.min, round(value))))
    return float(min(spec.max, max(spec.min, value)))


def cheap_family(space: Space, base: Config, limit: int = 5) -> List[Config]:
    """Genotyping variants of `base`: confidence both ways, indel prior, contamination, max alleles.

    PCR stays where `base` put it. These share `base`'s traversal.
    """
    axes: List[Tuple[str, Sequence]] = []
    name = "standard_min_confidence_threshold_for_calling"
    spec = space.specs.get(name)
    if spec is not None and name in space.names:
        cur = float(base.params.get(name, spec.default))
        grid = [v for v in CALL_CONF_GRID if spec.min <= v <= spec.max]
        above = [v for v in grid if v > cur + 1e-6]
        below = [v for v in grid if v < cur - 1e-6]
        picked = []
        if above:
            picked.append(min(above))
        if below:
            picked.append(max(below))
        if picked:
            axes.append((name, tuple(picked)))
    het = "indel_heterozygosity"
    hspec = space.specs.get(het)
    if hspec is not None and het in space.names:
        cur = float(base.params.get(het, hspec.default))
        axes.append((het, (_clamp_num(hspec, cur * 2.0), _clamp_num(hspec, cur / 2.0))))
    cont = "contamination_fraction_to_filter"
    cspec = space.specs.get(cont)
    if cspec is not None and cont in space.names:
        cur = float(base.params.get(cont, cspec.default))
        other = 0.0 if cur > 0 else next((v for v in CONTAMINATION_GRID if v > 0), None)
        if other is not None and cspec.min <= other <= cspec.max:
            axes.append((cont, (other,)))
    alt = "max_alternate_alleles"
    aspec = space.specs.get(alt)
    if aspec is not None and alt in space.names:
        cur = int(base.params.get(alt, aspec.default))
        other = 2 if cur != 2 else min(6, int(aspec.max))
        if aspec.min <= other <= aspec.max:
            axes.append((alt, (other,)))
    out: List[Config] = []
    seen = {base.key}
    for config in space._axis_variants(base, axes):
        if config.key in seen:
            continue
        seen.add(config.key)
        out.append(config)
        if len(out) >= limit:
            break
    return out


def _reset_genotyping(space: Space, config: Config, donor: Optional[Config]) -> Optional[Config]:
    if donor is None:
        return config
    params = dict(config.params)
    for name in GENOTYPING_PARAMS:
        if name in donor.params:
            params[name] = donor.params[name]
    return space._assemble(params)


def candidate_pool(space: Space, rng: random.Random, seeds: Sequence[Config], n_pool: int) -> List[Config]:
    """Skeletons spread over the structural parameters. Genotyping is copied from the first seed."""
    donor = seeds[0] if seeds else None
    names = structural_names(space)
    pool: List[Config] = []
    seen = {c.key for c in seeds}

    def add(config: Optional[Config]) -> None:
        if config is None or config.key in seen:
            return
        seen.add(config.key)
        pool.append(config)

    target_rand = max(1, int(n_pool * 0.50))
    tries = 0
    while len(pool) < target_rand and tries < n_pool * 8:
        tries += 1
        add(_reset_genotyping(space, space.random(rng), donor))
    target_mut = target_rand + max(1, int(n_pool * 0.25))
    while seeds and len(pool) < target_mut and tries < n_pool * 12:
        tries += 1
        parent = seeds[rng.randrange(len(seeds))]
        add(space.neighbour(parent, rng, n_params=rng.randint(1, 3), scale=0.2, pool=names))
    target_cross = target_mut + max(1, int(n_pool * 0.15))
    while len(seeds) >= 2 and names and len(pool) < target_cross and tries < n_pool * 16:
        tries += 1
        a, b = rng.sample(list(seeds), 2)
        params = dict(a.params)
        for name in rng.sample(names, min(len(names), rng.randint(1, 3))):
            params[name] = b.params.get(name, space.specs[name].default)
        add(space._assemble(params))
    target_edge = n_pool
    while names and len(pool) < target_edge and tries < n_pool * 20:
        tries += 1
        parent = seeds[rng.randrange(len(seeds))] if seeds else space.random(rng)
        name = rng.choice(names)
        spec = space.specs[name]
        if spec.allowed or spec.type == "bool":
            continue
        params = dict(parent.params)
        params[name] = space._from_unit(spec, 0.0 if rng.random() < 0.5 else 1.0)
        add(_reset_genotyping(space, space._assemble(params), donor))
    return pool


def farthest(space: Space, seeds: Sequence[Config], pool: Sequence[Config], n: int) -> List[Config]:
    """Greedy farthest-point order so the first workers are not copies of one ridge."""
    chosen: List[Config] = list(seeds)
    remain = [c for c in pool if all(c.key != s.key for s in chosen)]
    while len(chosen) < n and remain:
        best_i = 0
        best_d = -1.0
        for i, config in enumerate(remain):
            dist = min(struct_distance(space, config, kept) for kept in chosen)
            if dist > best_d:
                best_d, best_i = dist, i
        if best_d < 1e-9:
            break
        chosen.append(remain.pop(best_i))
    return chosen


def choose_submission(pairs: Sequence[Tuple[Config, float]], family_only: Set[str],
                     sign_of, margin: float) -> Optional[Config]:
    """First skeleton stays. A later one replaces it only when its signature differs and it leads by `margin`.

    `pairs` is completion order. Family members are never the submission.
    """
    chosen: Optional[Config] = None
    chosen_score = 0.0
    for config, score in pairs:
        if config.key in family_only:
            continue
        if chosen is None:
            chosen, chosen_score = config, score
            continue
        if sign_of(config) == sign_of(chosen):
            continue
        if score >= chosen_score + margin:
            chosen, chosen_score = config, score
    return chosen


def basin_race(ranges: ParamRanges, *, tool: str = "gatk", in_flight: int = 8, family: int = 6,
               seed: int = 0, band: float = BAND) -> Strategy:
    """Scout structurally different skeletons, then a short race inside the top local band.

    The submission rule is installed on the client: same signature never replaces the submission, and a
    different signature replaces it only when it leads by `band`. Family members are not submissions.
    """
    pinned = dict(PINNED.get(tool, {}))
    if tool == "gatk":
        pinned["pcr_indel_model"] = "NONE"
    space = Space(ranges, tool, pinned=pinned)
    rng = random.Random(seed)

    async def run(client: EvalClient) -> None:
        family_only: Set[str] = set()
        skel_score: Dict[Tuple, Tuple[float, Config]] = {}
        race_used: Dict[Tuple, int] = {}
        margin = [band]

        def choose(pairs: List[Tuple[Config, float]]) -> Optional[Config]:
            return choose_submission(pairs, family_only, lambda config: signature(space, config), margin[0])

        client.choose_submission = choose  # type: ignore[attr-defined]

        seeds: List[Config] = []
        seen: Set[str] = set()

        def take(params: Dict[str, object]) -> None:
            config = space._assemble(dict(params))
            if config is not None and config.key not in seen:
                seen.add(config.key)
                seeds.append(config)

        take(OFFICIAL_SEED)
        take({})
        for raw in client.store.portfolio(limit=4, tool=tool):
            take(dict(raw.params))

        deadline = client.time_left()
        n_scout = 48 if deadline is None else min(5000, max(64, in_flight * 6))
        ordered = farthest(space, seeds, candidate_pool(space, rng, seeds, n_scout), n_scout)
        queue: List[Config] = [c for c in ordered if c.key not in {s.key for s in seeds}]
        # Seeds go out first, at higher priority than the GATK-default baseline.
        first: List[Config] = list(seeds)

        proposed: Set[str] = set()
        wave_left: Dict[int, int] = {}
        wave_of: Dict[int, int] = {}
        wave_seq = [0]
        finish_sent = [False]

        def fresh(config: Config) -> bool:
            return (config.key not in proposed and config.key not in family_only
                    and (config.key, FULL) not in client.latest)

        def batch_of(base: Config) -> List[Config]:
            if not fresh(base):
                return []
            kids = [c for c in cheap_family(space, base, max(0, family - 1)) if fresh(c)]
            family_only.update(c.key for c in kids)
            proposed.add(base.key)
            proposed.update(c.key for c in kids)
            return [base] + kids

        def launch(batch: List[Config], priority: int) -> List:
            if not batch or client.exhausted:
                return []
            return list(client.propose_many(batch, priority=priority, tag="basin"))

        def occupy(batch: List[Config], futs: List) -> List:
            if not futs:
                return []
            wid = wave_seq[0]
            wave_seq[0] += 1
            wave_left[wid] = len(futs)
            for fut in futs:
                wave_of[id(fut)] = wid
            return futs

        def note(fut) -> None:
            wid = wave_of.pop(id(fut), None)
            if wid is not None:
                wave_left[wid] -= 1
                if wave_left[wid] <= 0:
                    wave_left.pop(wid, None)
            outcome = fut.result()
            config = outcome.config
            if outcome.ok and config is not None and config.key not in family_only:
                skel_score[signature(space, config)] = (float(outcome.score), config)

        def room(pending) -> bool:
            if client.exhausted:
                return False
            if len(wave_left) >= max(1, in_flight):
                return False
            idle = client.idle_slots()
            return True if idle is None else (idle > 0 or not pending)

        def leaders() -> List[Config]:
            if not skel_score:
                return []
            best = max(score for score, _ in skel_score.values())
            band_reps = [cfg for score, cfg in skel_score.values() if score >= best - margin[0]]
            band_reps.sort(key=lambda cfg: -skel_score[signature(space, cfg)][0])
            picked: List[Config] = []
            for config in band_reps:
                if len(picked) >= RACE_BASINS:
                    break
                if all(struct_distance(space, config, other) > 0.15 for other in picked):
                    picked.append(config)
            return picked

        def race_batch() -> List[Config]:
            reps = [cfg for cfg in leaders() if race_used.get(signature(space, cfg), 0) < RACE_PER_BASIN]
            if not reps:
                return []
            parent = min(reps, key=lambda cfg: race_used.get(signature(space, cfg), 0))
            child = space.neighbour(parent, rng, n_params=1, scale=0.12, pool=structural_names(space))
            if child is None:
                race_used[signature(space, parent)] = RACE_PER_BASIN
                return []
            race_used[signature(space, parent)] = race_used.get(signature(space, parent), 0) + 1
            return batch_of(child)

        def finish_batch() -> List[Config]:
            if finish_sent[0] or not skel_score:
                return []
            finish_sent[0] = True
            parent = choose([(cfg, score) for score, cfg in skel_score.values()])
            if parent is None:
                parent = max(skel_score.values(), key=lambda item: item[0])[1]
            return batch_of(parent)

        step = [0]

        def fill() -> List:
            if client.exhausted:
                return []
            left = client.time_left()
            if left is not None and left < FINISH_S:
                batch = finish_batch()
                return occupy(batch, launch(batch, 1))
            step[0] += 1
            # Three scouts, then one step inside a top-band basin. The scout queue is the rest of the round.
            batch: List[Config] = []
            if queue and (step[0] % 4 != 0 or not leaders()):
                batch = batch_of(queue.pop(0))
            if not batch:
                batch = race_batch()
            if not batch and queue:
                batch = batch_of(queue.pop(0))
            return occupy(batch, launch(batch, 0))

        opening: List = []
        for config in first:
            if len(wave_left) >= max(1, in_flight):
                queue.insert(0, config)
                continue
            batch = batch_of(config)
            opening.extend(occupy(batch, launch(batch, 11)))
        await run_streaming(client, room, fill, note, opening)

    return run
