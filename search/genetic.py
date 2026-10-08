"""Recombination: the opening's genetic generations and every merge. Owner: search.

Pairing and crossover are one operation at two sizes: a child takes each setting from one parent or the other.
Pairing tries every such combination of the settings where two parents differ (exhaustive when there are few);
a genetic generation samples them and adds mutation (settings redrawn over their whole range), which reaches values
neither parent has.
"""
from __future__ import annotations

import itertools
import random
from typing import Dict, List, Optional, Sequence, Tuple

from common.types import Config
from search.space import (ACTIVITY_PARAMS, ASSEMBLY_PARAMS, FULL_RUN_PARAMS, GENOTYPING_PARAMS, PAIRHMM_PARAMS,
                          Space)
from common.tuning import tune
from search.hills import Hills

# Tunable (common/tuning.py): merges.exhaustive_max (parents differing on this many settings or fewer: every
# combination), merges.random_combos (more: this many random ones), merges.mutated (plus this many children with
# MUTATE_AXES settings redrawn over their whole range).
MUTATE_AXES = 2
MUTATE_FRAC = 0.3           # share of a generation's children that are mutated
CROSS_NICHE = 0.3           # share of a generation's children whose parents come from two different hills
STAGES = (GENOTYPING_PARAMS, PAIRHMM_PARAMS, ASSEMBLY_PARAMS, ACTIVITY_PARAMS, FULL_RUN_PARAMS)


def _same(space: Space, name: str, a, b) -> bool:
    return space._same(space.specs[name], a, b)


def differing(space: Space, names: Sequence[str], a: Config, b: Config) -> List[str]:
    out = []
    for n in names:
        s = space.specs[n]
        if not _same(space, n, a.params.get(n, s.default), b.params.get(n, s.default)):
            out.append(n)
    return out


def take(space: Space, a: Config, b: Config, settings) -> Optional[Config]:
    """`a` with the given settings taken from `b`."""
    p = dict(a.params)
    for n in settings:
        p[n] = b.params.get(n, space.specs[n].default)
    return space._assemble(p)


def mutate(space: Space, c: Config, rng: random.Random, names: Sequence[str], k: int = MUTATE_AXES) -> Optional[Config]:
    p = dict(c.params)
    for n in rng.sample(list(names), min(k, len(names))):
        p[n] = space._from_unit(space.specs[n], rng.random())       # log-scaled settings stay log-uniform
    return space._assemble(p)


def uniform_child(space: Space, a: Config, b: Config, rng: random.Random, names: Sequence[str]) -> Optional[Config]:
    return take(space, a, b, [n for n in names if rng.random() < 0.5])


def block_child(space: Space, a: Config, b: Config, rng: random.Random, names: Sequence[str]) -> Optional[Config]:
    """Whole GATK stages from one parent or the other: keeps settings that work together inside a stage together."""
    chosen = set()
    for stage in STAGES:
        if rng.random() < 0.5:
            chosen.update(n for n in stage if n in names)
    return take(space, a, b, chosen)


def merge_children(space: Space, a: Config, b: Config, rng: random.Random, names: Sequence[str],
                   exhaustive_max: Optional[int] = None, n_random: Optional[int] = None,
                   n_mutated: Optional[int] = None) -> List[Config]:
    """The children a merge of two domains evaluates: every combination of the differing settings (or `n_random`
    random ones when there are more than `exhaustive_max`), plus `n_mutated` mutated combinations."""
    exhaustive_max = tune("merges.exhaustive_max") if exhaustive_max is None else exhaustive_max
    n_random = tune("merges.random_combos") if n_random is None else n_random
    n_mutated = tune("merges.mutated") if n_mutated is None else n_mutated
    diff = differing(space, names, a, b)
    out: Dict[str, Config] = {}
    if len(diff) <= exhaustive_max:
        subsets = (s for r in range(1, len(diff)) for s in itertools.combinations(diff, r))
    else:
        subsets = (tuple(n for n in diff if rng.random() < 0.5) for _ in range(n_random))
    for s in subsets:
        c = take(space, a, b, s)
        if c is not None and c.key not in (a.key, b.key):
            out[c.key] = c
    for _ in range(n_mutated):
        base = take(space, a, b, [n for n in diff if rng.random() < 0.5]) or a
        c = mutate(space, base, rng, names)
        if c is not None and c.key not in (a.key, b.key):
            out[c.key] = c
    return list(out.values())


def across_children(space: Space, parents: Sequence[Config], rng: random.Random, names: Sequence[str],
                    n: int) -> List[Config]:
    """One generation over every parent. Each child takes each setting from any of them, so a retired domain's
    value can sit next to a live domain's in the same child. `n` children, mutated at MUTATE_FRAC."""
    if len(parents) < 2 or n <= 0:
        return []
    lead = parents[0]
    out: Dict[str, Config] = {}
    tries = 0
    while len(out) < n and tries < 8 * n:
        tries += 1
        params = {name: rng.choice(parents).params.get(name, space.specs[name].default) for name in names}
        c = space._assemble({**lead.params, **params})
        if c is not None and rng.random() < MUTATE_FRAC:
            c = mutate(space, c, rng, names) or c
        if c is not None and all(c.key != p.key for p in parents):
            out.setdefault(c.key, c)
    return list(out.values())


def merge3_children(space: Space, parents: Sequence[Config], rng: random.Random, names: Sequence[str],
                    n_generation: int = 100) -> List[Config]:
    """A three-way (or more) merge: the leader paired with each other parent, plus one generation over all of them,
    so a setting from the third parent can join a combination of the first two."""
    lead, rest = parents[0], list(parents[1:])
    out: Dict[str, Config] = {}
    for other in rest:
        for c in merge_children(space, lead, other, rng, names, n_mutated=tune("merges.mutated") // 2):
            out[c.key] = c
    for _ in range(n_generation):
        p = {n: rng.choice(parents).params.get(n, space.specs[n].default) for n in names}
        c = space._assemble({**lead.params, **p})
        if c is not None and rng.random() < MUTATE_FRAC:
            c = mutate(space, c, rng, names) or c
        if c is not None and all(c.key != x.key for x in parents):
            out[c.key] = c
    return list(out.values())


def generation(space: Space, pool: Sequence[Tuple[Config, float]], rng: random.Random, n: int, hills: Hills,
               names: Sequence[str], niches: int = 12, per_niche: int = 15) -> List[Config]:
    """One genetic generation of `n` children from `pool` (config, score), with niching so several hills survive:
    the pool is grouped into up to `niches` hills, better hills get more children, every hill gets some; most children
    mate inside one hill, CROSS_NICHE across two; MUTATE_FRAC of them are mutated. Half uniform, half stage-block
    crossover."""
    ranked = sorted(pool, key=lambda t: -t[1])[:400]
    if len(ranked) < 2 or n <= 0:
        return []
    heads = hills.cluster(ranked, niches, band=float("inf"))
    members: List[List[Tuple[Config, float]]] = [[h] for h in heads]
    for c, s in ranked:
        if any(c.key == h[0].key for h in heads):
            continue
        best_i, best_d = 0, 10 ** 9
        for i, (h, _) in enumerate(heads):
            d = hills.distance(c, h)
            if d < best_d:
                best_i, best_d = i, d
        if len(members[best_i]) < per_niche:
            members[best_i].append((c, s))
    weights = [1.0 / (1 + r) for r in range(len(heads))]           # rank-weighted: the best hill most, every hill some
    out: Dict[str, Config] = {}
    tries = 0
    while len(out) < n and tries < 6 * n:
        tries += 1
        i = rng.choices(range(len(heads)), weights)[0]
        if rng.random() < CROSS_NICHE and len(heads) > 1:
            j = rng.choice([x for x in range(len(heads)) if x != i])
            a, b = rng.choice(members[i])[0], rng.choice(members[j])[0]
        else:
            if len(members[i]) < 2:
                continue
            a, b = (m[0] for m in rng.sample(members[i], 2))
        c = (uniform_child if rng.random() < 0.5 else block_child)(space, a, b, rng, names)
        if c is not None and rng.random() < MUTATE_FRAC:
            c = mutate(space, c, rng, names) or c
        if c is not None:
            out.setdefault(c.key, c)
    return list(out.values())


def cross_generation(space: Space, core: Sequence[Config], fp: Sequence[Config], rng: random.Random, n: int,
                     names: Sequence[str]) -> List[Config]:
    """`n` children, each taking settings from one high-core parent and one high-FP parent."""
    if n <= 0 or not core or not fp:
        return []
    out: Dict[str, Config] = {}
    tries = 0
    while len(out) < n and tries < 6 * n:
        tries += 1
        a, b = rng.choice(list(core)), rng.choice(list(fp))
        if a.key == b.key:
            continue
        c = (uniform_child if rng.random() < 0.5 else block_child)(space, a, b, rng, names)
        if c is not None and rng.random() < MUTATE_FRAC:
            c = mutate(space, c, rng, names) or c
        if c is not None:
            out.setdefault(c.key, c)
    return list(out.values())


def island_generation(space: Space, pool: Sequence[Tuple[Config, float]], island_of: Dict[str, int], islands: int,
                      rng: random.Random, n: int, hills: Hills, names: Sequence[str], min_far: int,
                      per_island: int = 15) -> List[Tuple[Config, int]]:
    """One generation on separate islands: `n` children as (config, island). The first call founds the islands (the
    best configs at least `min_far` settings apart); a config of no island joins the one whose best is nearest.
    Parents always come from one island, every island gets the same number of children, and a child that lands on
    another island's best (the hill rule) is redrawn, so the islands never grow into each other. `island_of` (config
    key -> island) is updated in place."""
    ranked = sorted(pool, key=lambda t: -t[1])[:400]
    if len(ranked) < 2 or n <= 0:
        return []
    if not island_of:
        for i, (c, _) in enumerate(hills.cluster(ranked, islands, band=float("inf"), min_far=min_far)):
            island_of[c.key] = i
    best: Dict[int, Config] = {}
    for c, _ in ranked:
        i = island_of.get(c.key)
        if i is not None and i not in best:
            best[i] = c
    if not best:
        return []
    members: Dict[int, List[Config]] = {i: [] for i in best}
    for c, _ in ranked:
        i = island_of.get(c.key)
        if i is None:
            i = min(best, key=lambda j: hills.distance(c, best[j]))
            island_of[c.key] = i
        if i in members and len(members[i]) < per_island:
            members[i].append(c)
    out: Dict[str, Tuple[Config, int]] = {}
    quota = max(1, n // len(members))
    for i, ms in members.items():
        made, tries = 0, 0
        while made < quota and tries < 6 * quota:
            tries += 1
            if len(ms) >= 2:
                a, b = rng.sample(ms, 2)
                c = (uniform_child if rng.random() < 0.5 else block_child)(space, a, b, rng, names)
                if c is not None and rng.random() < MUTATE_FRAC:
                    c = mutate(space, c, rng, names) or c
            else:
                c = mutate(space, ms[0], rng, names)
            if c is None or c.key in out:
                continue
            if any(j != i and hills.same(c, best[j]) for j in best):
                continue                                         # on another island's hill: redraw
            out[c.key] = (c, i)
            made += 1
    return list(out.values())
