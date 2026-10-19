"""A modest RF-UCB proposal stream. It ranks a candidate pool; it does not replace the genetic proposals.
step_b (2026-10-19): in Step B it takes climb.rf_ucb_share of the workers the grids leave to the GA.

Train on scores from this round only. Uncertainty is the standard deviation across trees, a heuristic,
not a calibrated interval. UCB is predicted score plus kappa times that deviation. The batch is then
thinned so two selected configs stay at least the seat distance apart.
"""
from __future__ import annotations

import asyncio
from typing import Awaitable, Callable, List, Sequence, Tuple

import numpy as np
from sklearn.ensemble import RandomForestRegressor

from common.types import Config
from search.axes import latin_hypercube
from search.distance import setting_distance

KAPPA = 2.0
SHARE = 0.15
MIN_TRAIN = 40
DISTANCE = 0.05
REFIT_S = 120.0            # Step B refits the forest at most this often (round clock)
TRAIN_MAX = 20000          # Step B training rows: the TOP_KEEP best and a random sample of the rest
TOP_KEEP = 5000


def share_of(raw: str) -> float:
    """`1` means the default share. A number is that share, capped at half the batch."""
    if raw == "1":
        return SHARE
    try:
        return min(0.5, max(0.0, float(raw)))
    except ValueError:
        return 0.0


def propose(ctx, n: int, kappa: float = KAPPA) -> List[Config]:
    """Up to `n` configs, highest UCB first, spread at least `DISTANCE` apart in settings."""
    scored = [(config, score) for config, score in ctx.scored()]
    if n <= 0 or len(scored) < MIN_TRAIN:
        ctx.say("rf-ucb: %d scored configs, no proposals", len(scored))
        return []
    forest = RandomForestRegressor(
        n_estimators=80, min_samples_leaf=3, max_features="sqrt", n_jobs=1, random_state=0)
    forest.fit([ctx.space.encode(config) for config, _score in scored], [score for _config, score in scored])
    pool = _pool(ctx, scored, max(400, n * 8))
    if not pool:
        ctx.say("rf-ucb: the candidate pool is empty")
        return []
    matrix = np.asarray([ctx.space.encode(config) for config in pool], dtype=float)
    predictions = np.asarray([tree.predict(matrix) for tree in forest.estimators_], dtype=float)
    mean = predictions.mean(axis=0)
    deviation = predictions.std(axis=0)
    order = np.argsort(-(mean + kappa * deviation))
    chosen: List[Config] = []
    picked_mean = picked_dev = 0.0
    for index in order:
        config = pool[int(index)]
        if not ctx.fresh(config):
            continue
        if any(setting_distance(ctx.space, config, other) < DISTANCE for other in chosen):
            continue
        chosen.append(config)
        picked_mean += float(mean[int(index)])
        picked_dev += float(deviation[int(index)])
        if len(chosen) >= n:
            break
    if chosen:
        ctx.say("rf-ucb: trained on %d, pool %d, selected %d, kappa %.1f, mean prediction %.3f, mean uncertainty %.3f",
                len(scored), len(pool), len(chosen), kappa, picked_mean / len(chosen), picked_dev / len(chosen))
    else:
        ctx.say("rf-ucb: trained on %d, pool %d, none selected", len(scored), len(pool))
    return chosen


def _fit(X, y) -> RandomForestRegressor:
    forest = RandomForestRegressor(n_estimators=80, min_samples_leaf=3, max_features="sqrt", n_jobs=4, random_state=0)
    forest.fit(X, y)
    return forest


def _refit(ctx, state: dict, scored: Sequence[Tuple[Config, float]]) -> None:
    """Start a fit off the event loop: the TOP_KEEP best scores and a sample of the rest, encodings cached."""
    rows = sorted(scored, key=lambda item: -item[1])
    if len(rows) > TRAIN_MAX:
        rows = rows[:TOP_KEEP] + ctx.rng.sample(rows[TOP_KEEP:], TRAIN_MAX - TOP_KEEP)
    enc = state["enc"]
    X = []
    for config, _score in rows:
        v = enc.get(config.key)
        if v is None:
            v = enc[config.key] = ctx.space.encode(config)
        X.append(v)
    state["task"] = asyncio.ensure_future(asyncio.to_thread(_fit, X, [score for _config, score in rows]))
    state["t"], state["n"], state["rows"] = ctx.clock(), len(scored), len(rows)


async def step_b(ctx, n: int, accept: Callable[[Config], Awaitable[bool]], kappa: float = 0.0) -> List[Config]:
    """Up to `n` Step B configs, highest predicted score + kappa x spread first, `accept`ed and DISTANCE apart.
    The forest is refit off the event loop at most every REFIT_S; until the first fit is back, nothing is proposed
    and the GA keeps those workers."""
    if n <= 0:
        return []
    state = getattr(ctx, "_rf", None)
    if state is None:
        state = ctx._rf = {"forest": None, "task": None, "t": None, "n": 0, "rows": 0, "enc": {}, "fits": 0,
                           "sent": 0}
    task = state["task"]
    if task is not None and task.done():
        state["task"] = None
        try:
            state["forest"] = task.result()
            state["fits"] += 1
            ctx.say("rf-ucb: forest %d trained on %d of %d scored configs", state["fits"], state["rows"], state["n"])
        except Exception as e:  # noqa: BLE001 - a failed fit keeps the last forest
            ctx.say("rf-ucb: fit failed (%s)", e)
    scored = list(ctx.scored())
    if len(scored) < MIN_TRAIN:
        return []
    if state["task"] is None and (state["t"] is None or
                                  (ctx.clock() - state["t"] >= REFIT_S and len(scored) > state["n"])):
        _refit(ctx, state, scored)
    forest = state["forest"]
    if forest is None:
        return []
    pool = _pool(ctx, scored, max(400, n * 8))
    if not pool:
        return []
    matrix = np.asarray([ctx.space.encode(config) for config in pool], dtype=float)
    predictions = np.asarray([tree.predict(matrix) for tree in forest.estimators_], dtype=float)
    order = np.argsort(-(predictions.mean(axis=0) + kappa * predictions.std(axis=0)))
    chosen: List[Config] = []
    for index in order:
        config = pool[int(index)]
        if any(setting_distance(ctx.space, config, other) < DISTANCE for other in chosen):
            continue
        if not await accept(config):
            continue
        chosen.append(config)
        if len(chosen) >= n:
            break
    state["sent"] += len(chosen)
    return chosen


def _pool(ctx, scored: Sequence[Tuple[Config, float]], size: int) -> List[Config]:
    """Mutations of the current best scores, plus a Latin sample of the whole grid."""
    ranked = sorted(scored, key=lambda item: -item[1])[:30]
    out: List[Config] = []
    seen = {config.key for config, _score in scored}
    for config, _score in ranked:
        for _ in range(8):
            child = _mutate(ctx, config, ctx.rng.randint(1, 3))
            if child is None or child.key in seen:
                continue
            seen.add(child.key)
            out.append(child)
            if len(out) >= size:
                return out
    for config in latin_hypercube(ctx.space, ctx.rng, max(0, size - len(out))):
        if config.key in seen:
            continue
        seen.add(config.key)
        out.append(config)
        if len(out) >= size:
            break
    return out


def _mutate(ctx, config: Config, axes: int) -> Config:
    params = dict(config.params)
    movable = [name for name in ctx.space.names if ctx.ladders.movable(name)]
    if not movable:
        return None
    for name in ctx.rng.sample(movable, min(axes, len(movable))):
        rungs = ctx.ladders.rungs.get(name) or []
        if len(rungs) < 2:
            continue
        params[name] = ctx.rng.choice([value for value in rungs if value != params.get(name)] or rungs)
    return ctx.space._assemble(params)
