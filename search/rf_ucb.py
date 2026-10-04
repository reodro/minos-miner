"""Tiered RF-UCB: two random-forest surrogates, cost-normalized upper confidence, designed cheap families.

Loop: keep every worker busy. Refit when the shared result set has grown. Fit one forest on `core` and one on
`fp_per_target`, and rank candidates by the v2 formula of those predictions (mean + kappa * spread across trees)
divided by what the candidate costs. Kappa starts high and decays toward exploitation.

Candidates move the expensive parameters (full run, assembly, activity). Each pick is proposed with a designed
family of genotyping and PairHMM variants, one `propose_many`, so the executor can score them in the same GATK
run. A cheap sweep of the current best is taken first when its cost-normalized value wins, which is one batch,
not a stream of full runs. After about 100 scores, parameters the trees barely use are frozen at the incumbent.
In the last minutes only genotyping variants of the incumbent are proposed: a new full run would not finish.

The initial sample is the GATK default, one step off it on each expensive axis, and two starts outside that
neighborhood, plus expensive keys from previous rounds (configs only; their scores are not reused).

The surrogate is scikit-learn's ExtraTreesRegressor. Tuning knobs are the arguments.
"""
from __future__ import annotations

import random
from typing import List, Optional, Set, Tuple

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor

from common.params import ParamRanges
from common.scoring import score_from_parts
from common.types import FULL, Config
from search.base import Strategy
from search.client import EvalClient, Outcome
from search.family import designed_family
from search.schedule import run_streaming
from search.space import Space


def rf_ucb(ranges: ParamRanges, *, tool: str = "gatk", n_init: int = 16, in_flight: int = 8, n_candidates: int = 300,
           kappa: float = 1.0, kappa_end: float = 0.15, kappa_horizon: int = 200, top_k: int = 8,
           explore_frac: float = 0.25, family: int = 6, family_max_cost: float = 0.15, n_trees: int = 30,
           max_pending: int = 256, stable_after: int = 40, freeze_after: int = 100, late_s: float = 180.0,
           portfolio: int = 4, seed: int = 0, tag: str = "rf_ucb") -> Strategy:
    space = Space(ranges, tool)
    rng = random.Random(seed)

    async def run(client: EvalClient) -> None:
        failed: Set[str] = set()
        proposed: Set[str] = set()
        core_model: List[Optional[ExtraTreesRegressor]] = [None]
        fp_model: List[Optional[ExtraTreesRegressor]] = [None]
        fitted_on = [0]
        frozen: Set[str] = set()

        def rows() -> Tuple[List[Tuple[Config, float, float, float]], Set[str]]:
            good: List[Tuple[Config, float, float, float]] = []
            bad = set(failed)
            for (key, fid), result in client.latest.items():
                if fid != FULL or result.score is None or key not in client.configs:
                    continue
                if client.configs[key].tool != tool:
                    continue
                if result.gate_passed is False:
                    bad.add(key)
                    continue
                core = result.core if result.core is not None else result.score / 100.0
                fp = (result.region_fp / result.n_targets) if result.n_targets else 0.0
                good.append((client.configs[key], float(core), float(fp), float(result.score)))
            return good, bad

        def observations():
            good, bad = rows()
            if not good:
                return [], [], []
            worst_core = min(c for _, c, _, _ in good)
            worst_fp = max(fp for _, _, fp, _ in good)
            configs = [c for c, _, _, _ in good]
            cores = [c for _, c, _, _ in good]
            fps = [fp for _, _, fp, _ in good]
            for key in bad:
                if key in client.configs and client.configs[key].tool == tool and key not in {c.key for c in configs}:
                    configs.append(client.configs[key])
                    cores.append(worst_core)
                    fps.append(worst_fp)
            return configs, cores, fps

        def incumbent() -> Optional[Config]:
            good, _ = rows()
            if not good:
                return None
            return max(good, key=lambda t: t[3])[0]

        def kappa_now(n_obs: int) -> float:
            if kappa_horizon <= 0:
                return kappa
            t = min(1.0, n_obs / kappa_horizon)
            return kappa + (kappa_end - kappa) * t

        def refit_if_due():
            configs, cores, fps = observations()
            if len(cores) < 4:
                return
            if core_model[0] is not None and len(cores) - fitted_on[0] < max(1, len(cores) // 20):
                return
            x = np.array([space.encode(c) for c in configs])
            y_core, y_fp = np.array(cores), np.array(fps)
            state = rng.randrange(1 << 30)
            core_model[0] = ExtraTreesRegressor(
                n_estimators=n_trees, min_samples_leaf=3, max_features="sqrt", random_state=state).fit(x, y_core)
            fp_model[0] = ExtraTreesRegressor(
                n_estimators=n_trees, min_samples_leaf=3, max_features="sqrt", random_state=state + 1).fit(x, y_fp)
            fitted_on[0] = len(cores)
            good, _ = rows()
            if len(good) >= freeze_after:
                imp = 0.5 * (core_model[0].feature_importances_ + fp_model[0].feature_importances_)
                cutoff = 0.5 / max(len(space.names), 1)
                nxt = {n for n, v in zip(space.names, imp) if v < cutoff}
                frozen.clear()
                if len(nxt) < len(space.names) - 1:
                    frozen.update(nxt)

        def predict_ucb(cs: List[Config], n_obs: int) -> List[float]:
            x = np.array([space.encode(c) for c in cs])
            core_t = np.stack([t.predict(x) for t in core_model[0].estimators_])
            fp_t = np.stack([t.predict(x) for t in fp_model[0].estimators_])
            per_tree = np.vectorize(score_from_parts)(core_t, fp_t)
            k = kappa_now(n_obs)
            return list(per_tree.mean(axis=0) + k * per_tree.std(axis=0))

        def fresh(c: Optional[Config]) -> bool:
            return c is not None and c.key not in proposed and (c.key, FULL) not in client.latest and c.key not in failed

        def move_cost(anchor: Optional[Config], c: Config, explore: bool) -> float:
            if explore or anchor is None:
                return 1.0
            priced = client.variant_cost(anchor, c)
            if 0 < priced < 1.0:
                return priced
            return 1.0

        def late() -> bool:
            left = client.time_left()
            return left is not None and left < late_s

        def choose() -> List[Config]:
            inc = incumbent()
            if late():
                if inc is None:
                    return []
                return designed_family(
                    space, inc, family,
                    lambda c: fresh(c) and client.variant_cost(inc, c) <= 0.05,
                    genotyping_only=True)
            if core_model[0] is None or fp_model[0] is None:
                c = space.random(rng)
                return [c] if fresh(c) else []
            _configs, cores, _fps = observations()
            n_obs = len(cores)
            good, _ = rows()
            ranked = [c for c, _, _, _ in sorted(good, key=lambda t: -t[3])][:top_k]
            free = [n for n in space.names if n not in frozen] or list(space.names)
            ints = [n for n in space.expensive_ints if n not in frozen]
            cands: List[Config] = []
            costs: List[float] = []
            for _ in range(n_candidates):
                explore = rng.random() < explore_frac or not ranked
                if explore:
                    c = space.random(rng)
                elif inc is not None and n_obs >= stable_after and ints and rng.random() < 0.5:
                    # One expensive integer at a time, then a two-parameter step so a joint move is still possible.
                    c = space.neighbour(inc, rng, n_params=1, scale=0.12, pool=ints)
                    if c is None:
                        c = space.neighbour(rng.choice(ranked), rng, n_params=2, scale=0.15, pool=free)
                else:
                    anchor = rng.choice(ranked) if ranked else inc
                    c = space.neighbour(anchor, rng, n_params=2, scale=0.15, pool=free) if anchor is not None else None
                c = space.clamp_frozen(c, inc, frozen)
                if fresh(c):
                    cands.append(c)
                    costs.append(move_cost(inc, c, explore))
            cheap = []
            if inc is not None:
                cheap = designed_family(
                    space, inc, family,
                    lambda c: fresh(c) and client.variant_cost(inc, c) <= family_max_cost)
            if cheap:
                util_cheap = [u / max(client.variant_cost(inc, c), 1e-3) for u, c in zip(predict_ucb(cheap, n_obs), cheap)]
                costly = [(c, cost) for c, cost in zip(cands, costs) if cost > family_max_cost]
                best_cheap = max(util_cheap)
                best_costly = None
                if costly:
                    ucb = predict_ucb([c for c, _ in costly], n_obs)
                    best_costly = max(u / max(cost, 1e-3) for u, (_, cost) in zip(ucb, costly))
                if best_costly is None or best_cheap >= best_costly:
                    order = [c for _, c in sorted(zip(util_cheap, cheap), key=lambda t: -t[0])]
                    return order
            costly = [(c, cost) for c, cost in zip(cands, costs) if cost > family_max_cost]
            if not costly:
                return []
            ucb = predict_ucb([c for c, _ in costly], n_obs)
            utilities = [u / max(cost, 1e-3) for u, (_, cost) in zip(ucb, costly)]
            best = costly[max(range(len(costly)), key=utilities.__getitem__)][0]
            fam = designed_family(
                space, best, family,
                lambda c: fresh(c) and c.key != best.key and client.variant_cost(best, c) <= family_max_cost)
            return [best] + fam

        def note(fut):
            o: Outcome = fut.result()
            if not o.ok and o.config is not None and not (o.error or "").startswith("budget"):
                failed.add(o.config.key)

        def room(pending) -> bool:
            """Propose while a worker is free. A family is one job, so the cap is in configs: at least
            `max_pending`, and more on a big fleet, or a fixed cap would leave workers idle."""
            cap = max(max_pending, (family + 1) * (client.capacity() or 0))
            if len(pending) >= cap:
                return False
            idle = client.idle_slots()
            return (len(pending) < in_flight) if idle is None else (idle > 0 or not pending)

        def fill():
            refit_if_due()
            batch = choose()
            if not batch:
                return []
            proposed.update(c.key for c in batch)
            return list(client.propose_many(batch, tag=tag))

        initial = []
        if not late():
            seen = set()
            seeded = list(space.initial_design(rng, n_init))
            for raw in client.store.portfolio(limit=portfolio, tool=tool):
                warmed = space.with_cheap_defaults(raw)
                if warmed is not None:
                    seeded.append(warmed)
            for c in seeded:
                if fresh(c) and c.key not in seen:
                    seen.add(c.key)
                    initial.append(c)
            proposed.update(c.key for c in initial)
        initial_futs = list(client.propose_many(initial, tag=tag)) if initial else []
        await run_streaming(client, room, fill, note, initial_futs)

    return run
