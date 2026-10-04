"""Several trust regions, a core/FP forest inside each one, and restarts that stay outside every live radius.

The validator score is not a smooth slope. It is an integer step (one call is a few thousandths), a cliff (a
failed ti/tv or het/hom gate zeros the score), and a false-positive term that is steeper on the live estimate
of N than on the released truth. A finite-difference step usually sees a zero slope, or a cliff, and walks
into the conservative local maximum.

So a region does not estimate a gradient. It proposes one- and two-parameter steps that stay inside its
radius, and a forest on `core` and `fp_per_target` ranks them. Two candidates within half a point are ordered
by core, because a gap that small is what the live N error does to the false-positive term. A gate failure is
thrown away: it is not a score of 0 and it does not move the center. The radius shrinks only after many local
scores with no improvement, which is a real plateau, not one tied probe.

Restarts keep running for the whole round and must land outside every live radius. A restart that beats the
worst region takes that region over. The last minutes propose only genotyping variants.
"""
from __future__ import annotations

import math
import random
from typing import Dict, List, Optional, Sequence, Set, Tuple

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


class _Region:
    def __init__(self, ident: int, center: Config, score: float, radius: float):
        self.ident = ident
        self.center = center
        self.best_config = center
        self.best_score = score
        self.radius = radius
        self.stall = 0
        self.dead = False


def gradient_regions(ranges: ParamRanges, *, tool: str = "gatk", n_init: int = 16, n_regions: int = 0,
                     in_flight: int = 8, radius0: float = 0.20, radius_min: float = 0.05, shrink: float = 0.7,
                     flat_limit: int = 24, restart_frac: float = 0.30, pair_frac: float = 0.10, family: int = 6,
                     family_max_cost: float = 0.15, late_s: float = 180.0, separation: float = 0.12,
                     score_band: float = 0.5, n_trees: int = 20, kappa: float = 0.35,
                     start: Optional[Sequence[Config]] = None, max_pending: int = 256, seed: int = 0,
                     tag: str = "regions") -> Strategy:
    """`n_regions` 0 scales with `in_flight` (at least 2, at most 20). `start` replaces the default neighborhood."""
    space = Space(ranges, tool)
    rng = random.Random(seed)
    dist_dims = [n for n in space.expensive if space.specs[n].type in ("int", "float")]
    target = n_regions if n_regions > 0 else max(2, min(20, max(1, in_flight) // 4))

    async def run(client: EvalClient) -> None:
        regions: List[_Region] = []
        scored: List[Tuple[float, float, Config]] = []
        pending_far: List[Tuple[float, Config]] = []
        failed: Set[str] = set()
        proposed: Set[str] = set()
        local_of: Dict[str, int] = {}
        counts = {"local": 0, "restart": 0, "pair": 0}
        born = [0]
        next_id = [0]
        rr = [0]
        core_model: List[Optional[ExtraTreesRegressor]] = [None]
        fp_model: List[Optional[ExtraTreesRegressor]] = [None]
        fitted_on = [0]
        run.state = {"born": 0, "alive": 0, "waves": 0}  # type: ignore[attr-defined]

        def publish() -> None:
            run.state["born"] = born[0]  # type: ignore[attr-defined]
            run.state["alive"] = sum(1 for r in regions if not r.dead)  # type: ignore[attr-defined]
            run.state["waves"] = sum(r.stall for r in regions)  # type: ignore[attr-defined]

        def vec(config: Config) -> List[float]:
            return [space._to_unit(space.specs[n], config.params.get(n, space.specs[n].default)) for n in dist_dims]

        def dist(a: Config, b: Config) -> float:
            return math.sqrt(sum((x - y) ** 2 for x, y in zip(vec(a), vec(b))))

        def fresh(config: Optional[Config]) -> bool:
            return (config is not None and config.key not in proposed and config.key not in failed
                    and (config.key, FULL) not in client.latest)

        def far_from(config: Config, ignoring: Optional[_Region] = None) -> bool:
            for region in regions:
                if region.dead or region is ignoring:
                    continue
                if dist(config, region.center) <= region.radius + 1e-6:
                    return False
            return True

        def inside(config: Config, region: _Region) -> bool:
            return dist(config, region.center) <= region.radius + 1e-6

        def add_region(config: Config, score: float) -> None:
            region = _Region(next_id[0], config, score, radius0)
            next_id[0] += 1
            regions.append(region)
            pending_far[:] = [(s, c) for s, c in pending_far if c.key != config.key]
            born[0] += 1
            publish()

        def remember_far(config: Config, score: float) -> None:
            if not any(not r.dead for r in regions):
                return
            if any(config.key == c.key for _, c in pending_far):
                return
            if any(config.key == r.center.key for r in regions):
                return
            pending_far.append((score, config))
            pending_far.sort(key=lambda item: -item[0])
            del pending_far[8:]

        def live_rows() -> List[Tuple[Config, float, float, float]]:
            good = []
            for (key, fid), result in client.latest.items():
                if fid != FULL or result.score is None or result.gate_passed is False or key not in client.configs:
                    continue
                if client.configs[key].tool != tool:
                    continue
                core = result.core if result.core is not None else result.score / 100.0
                fp = (result.region_fp / result.n_targets) if result.n_targets else 0.0
                good.append((client.configs[key], float(core), float(fp), float(result.score)))
            return good

        def refit() -> None:
            good = live_rows()
            if len(good) < 8:
                return
            if core_model[0] is not None and len(good) - fitted_on[0] < max(4, len(good) // 10):
                return
            x = np.array([space.encode(c) for c, _, _, _ in good])
            state = rng.randrange(1 << 30)
            core_model[0] = ExtraTreesRegressor(
                n_estimators=n_trees, min_samples_leaf=3, max_features="sqrt", random_state=state).fit(
                x, np.array([c for _, c, _, _ in good]))
            fp_model[0] = ExtraTreesRegressor(
                n_estimators=n_trees, min_samples_leaf=3, max_features="sqrt", random_state=state + 1).fit(
                x, np.array([fp for _, _, fp, _ in good]))
            fitted_on[0] = len(good)

        def rank_key(config: Config) -> Tuple[float, float]:
            """Higher is better. Within `score_band`, the higher core wins: that gap is the live N error."""
            if core_model[0] is None or fp_model[0] is None:
                return (0.0, 0.0)
            x = np.array([space.encode(config)])
            cores = np.stack([t.predict(x) for t in core_model[0].estimators_])
            fps = np.stack([t.predict(x) for t in fp_model[0].estimators_])
            per_tree = np.vectorize(score_from_parts)(cores, fps)
            return (float(per_tree.mean() + kappa * float(per_tree.std())), float(cores.mean()))

        def better(a: Config, b: Config) -> bool:
            sa, ca = rank_key(a)
            sb, cb = rank_key(b)
            if sa > sb + score_band:
                return True
            if sb > sa + score_band:
                return False
            return ca > cb

        def shrink_region(region: _Region) -> None:
            region.stall = 0
            if region.radius <= radius_min + 1e-12:
                region.dead = True
                publish()
                return
            region.radius = max(radius_min, region.radius * shrink)
            for score, config in list(pending_far):
                if score <= region.best_score or dist(config, region.center) <= region.radius + 1e-6:
                    continue
                if not far_from(config, ignoring=region):
                    continue
                pending_far.remove((score, config))
                region.center = config
                region.best_config = config
                region.best_score = score
                region.radius = radius0
                return

        def note(fut) -> None:
            outcome: Outcome = fut.result()
            config = outcome.config
            if config is None:
                return
            gated = outcome.ok and outcome.result is not None and outcome.result.gate_passed is False
            if not outcome.ok or gated:
                if not (outcome.error or "").startswith("budget"):
                    failed.add(config.key)
                ident = local_of.pop(config.key, None)
                region = next((r for r in regions if r.ident == ident), None)
                if region is not None and not region.dead:
                    region.stall += 1
                    if region.stall >= flat_limit:
                        shrink_region(region)
                return
            core = outcome.result.core if outcome.result.core is not None else outcome.score / 100.0
            scored.append((outcome.score, float(core), config))
            ident = local_of.pop(config.key, None)
            owner = next((r for r in regions if r.ident == ident), None)
            if owner is None:
                owner = next((r for r in regions if not r.dead and inside(config, r)), None)
            if owner is not None and not owner.dead:
                if outcome.score > owner.best_score + 1e-9:
                    owner.best_score = outcome.score
                    owner.best_config = config
                    owner.center = config
                    owner.stall = 0
                elif ident is not None:
                    owner.stall += 1
                    if owner.stall >= flat_limit:
                        shrink_region(owner)
            elif far_from(config):
                remember_far(config, outcome.score)
            maybe_seed()
            publish()

        def maybe_seed() -> None:
            live = [r for r in regions if not r.dead]
            need = 1 if start is not None else min(8, n_init)
            if len(live) >= target or len(scored) < need:
                return
            pool = sorted(scored, key=lambda item: -item[0])
            if not live and pool:
                add_region(pool[0][2], pool[0][0])
                live = [r for r in regions if not r.dead]
            while len(live) < target:
                best: Optional[Tuple[float, Config]] = None
                best_d = -1.0
                for score, _core, config in pool:
                    if any(config.key == r.center.key for r in live):
                        continue
                    d = min(dist(config, r.center) for r in live)
                    if d > best_d:
                        best, best_d = (score, config), d
                if best is None or best_d < separation:
                    break
                add_region(best[1], best[0])
                live = [r for r in regions if not r.dead]

        def fill_slots() -> None:
            regions[:] = [r for r in regions if not r.dead]
            i = 0
            while len(regions) < target and i < len(pending_far):
                score, config = pending_far[i]
                if any(config.key == r.center.key for r in regions) or not far_from(config):
                    i += 1
                    continue
                pending_far.pop(i)
                add_region(config, score)

        def late() -> bool:
            left = client.time_left()
            return left is not None and left < late_s

        def archive() -> Optional[Config]:
            if not scored:
                return None
            return max(scored, key=lambda item: item[0])[2]

        def batch_for(config: Optional[Config], genotyping_only: bool = False) -> List[Config]:
            if not fresh(config):
                return []
            cap = 0.05 if genotyping_only else family_max_cost
            fam = designed_family(
                space, config, family,
                lambda c: fresh(c) and c.key != config.key and client.variant_cost(config, c) <= cap,
                genotyping_only=genotyping_only)
            return [config] + fam

        def propose_restart() -> Optional[Config]:
            best, best_d = None, -1.0
            for _ in range(40):
                config = space.random(rng)
                if not fresh(config):
                    continue
                if far_from(config):
                    return config
                d = min((dist(config, r.center) for r in regions if not r.dead), default=1.0)
                if d > best_d:
                    best, best_d = config, d
            return best

        def propose_pair() -> Optional[Config]:
            live = [r for r in regions if not r.dead]
            if not live:
                return propose_restart()
            anchor = rng.choice(live).best_config
            return space.neighbour(anchor, rng, n_params=2, scale=0.20, pool=space.expensive)

        def propose_local() -> Optional[Tuple[_Region, Config]]:
            live = [r for r in regions if not r.dead]
            if not live:
                return None
            start_i = rr[0]
            for k in range(len(live)):
                region = live[(start_i + k) % len(live)]
                cands = []
                for _ in range(16):
                    n_params = 1 if rng.random() < 0.7 else 2
                    config = space.neighbour(region.best_config, rng, n_params=n_params,
                                             scale=max(region.radius / 2.0, 0.05), pool=dist_dims or space.expensive)
                    if config is not None and fresh(config) and inside(config, region):
                        cands.append(config)
                if not cands:
                    continue
                best = cands[0]
                for config in cands[1:]:
                    if better(config, best):
                        best = config
                rr[0] = (start_i + k + 1) % len(live)
                return region, best
            return None

        def lane() -> str:
            if not any(not r.dead for r in regions):
                return "restart" if restart_frac > 0 else "local"
            total = sum(counts.values()) + 1
            targets = {"local": max(0.0, 1.0 - restart_frac - pair_frac), "restart": restart_frac, "pair": pair_frac}
            return max(targets, key=lambda name: targets[name] - counts[name] / total)

        def choose() -> List[Config]:
            fill_slots()
            refit()
            if late():
                best = archive()
                if best is None:
                    return []
                return designed_family(
                    space, best, family,
                    lambda c: fresh(c) and client.variant_cost(best, c) <= 0.05,
                    genotyping_only=True)
            proposers = ("local", "restart", "pair")
            order = [lane()]
            for name in proposers:
                if name not in order and (name != "restart" or restart_frac > 0) and (name != "pair" or pair_frac > 0):
                    order.append(name)
            for kind in order:
                for _ in range(12):
                    if kind == "local":
                        picked = propose_local()
                        if picked is None:
                            break
                        region, config = picked
                        chosen = batch_for(config)
                        if not chosen:
                            continue
                        local_of[config.key] = region.ident
                        counts[kind] += 1
                        return chosen
                    config = propose_restart() if kind == "restart" else propose_pair()
                    chosen = batch_for(config)
                    if chosen:
                        counts[kind] += 1
                        return chosen
            return []

        def room(pending) -> bool:
            cap = max(max_pending, (family + 1) * (client.capacity() or 0))
            if len(pending) >= cap:
                return False
            idle = client.idle_slots()
            return (len(pending) < in_flight) if idle is None else (idle > 0 or not pending)

        def fill():
            chosen = choose()
            if not chosen:
                return []
            proposed.update(c.key for c in chosen)
            return list(client.propose_many(chosen, tag=tag))

        initial = list(start) if start is not None else list(space.initial_design(rng, n_init))
        initial = [c for c in initial if fresh(c)]
        proposed.update(c.key for c in initial)
        initial_futs = list(client.propose_many(initial, tag=tag)) if initial else []
        await run_streaming(client, room, fill, note, initial_futs)
        publish()

    run.state = {"born": 0, "alive": 0, "waves": 0}  # type: ignore[attr-defined]
    return run
