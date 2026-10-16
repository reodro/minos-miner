"""Separate mountains (v9.30).

Step A places a dummy on the crest of a local quadratic at each apparent peak. One climber
sits on the highest fingerprint. A later fingerprint above that climber is valley-tested
against it: scores that already lie on the callset path, then one midpoint, then one quarter
when that midpoint sits inside the noise band. A dip of 0.10 moves the climber onto the new
fingerprint. No dip, and the callsets are farther than 0.10, leaves the lower one.

Step B compasses the highest score until the fine-tune. Up to four axes, ordered by the
one-setting score change near the centre. A face within 0.05 is not live. Two or more live
faces send their joint cube, and that cube is scored before the centre moves. Free cores look
for a fingerprint 0.10 away, then probe settings the compass has not tried. A core with
neither of those waits. The fine-tune is the same compass and does not take a new fingerprint.
"""
from __future__ import annotations

import asyncio
import itertools
import logging
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from common.types import Config
from search.engine_tag import VERSION
from search.lines import single
from search.piles import FLOOR, GAP, RADIUS, _owned, absorb_parents, mating_pool, refresh

log = logging.getLogger(__name__)

VALLEY = 0.10
WOBBLE = 0.02
SEATS = 1
FACE_MARGIN = 0.05          # a compass face this close to the centre is not a live axis
COMPASS_AXES = 4
COMPASS_SLOTS = 16          # cores left idle so the next compass starts with its own heads
RANGE_GAP = 0.50            # a candidate within this of the best score is inside Step A's range
RANGE_BATCH = 8             # points Step A may add while one range batch is still out
SCOUT_WAVES = 2
SCOUT_N = 50
TENTHS = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
NEIGHBORS = 6
QUAD_NEAR = 48
QUAD_AXES = 2
CHILDREN = 4


def enabled() -> bool:
    return (VERSION.startswith("v9.29") or VERSION.startswith("v9.30")
            or VERSION.startswith("v9.31") or VERSION.startswith("v9.32")
            or VERSION.startswith("v9.33") or VERSION.startswith("v9.34"))


def _vec(ctx, config: Config) -> List[float]:
    return list(ctx.space.encode(config))


def _dist2(a: Sequence[float], b: Sequence[float]) -> float:
    return sum((x - y) * (x - y) for x, y in zip(a, b))


def apparent_peaks(scored: Sequence[Tuple[Config, float, Sequence[float]]]) -> List[Config]:
    """A sample higher than its six nearest samples. With fewer samples, higher than every other one."""
    n = len(scored)
    if n < 2:
        return []
    mat = np.asarray([vec for _, _, vec in scored], dtype=float)
    scores = np.asarray([score for _, score, _ in scored], dtype=float)
    delta = mat[:, None, :] - mat[None, :, :]
    dist2 = np.sum(delta * delta, axis=2)
    np.fill_diagonal(dist2, np.inf)
    k = min(NEIGHBORS, n - 1)
    nearest = np.argpartition(dist2, k - 1, axis=1)[:, :k]
    return [scored[i][0] for i in range(n) if np.all(scores[i] > scores[nearest[i]])]


def _quadratic_columns(u: np.ndarray) -> np.ndarray:
    """1, the coordinates, then u_i u_j for i <= j."""
    d = int(u.shape[0])
    cols = [1.0]
    cols.extend(float(x) for x in u)
    for i in range(d):
        for j in range(i, d):
            cols.append(float(u[i] * u[j]))
    return np.asarray(cols, dtype=float)


def _quadratic_crest(origin: Sequence[float], score: float, rows: Sequence[Tuple[Sequence[float], float]]) -> Optional[np.ndarray]:
    """The crest of a downward bowl fitted in the two strongest directions. None when the bowl does not open downward
    or the crest falls outside the samples."""
    if len(rows) < 5:
        return None
    origin_v = np.asarray(origin, dtype=float)
    offsets = np.asarray([np.asarray(vec, dtype=float) - origin_v for vec, _ in rows], dtype=float)
    try:
        _u, singular, vt = np.linalg.svd(offsets, full_matrices=False)
    except np.linalg.LinAlgError:
        return None
    if singular.size == 0 or singular[0] < 1e-8:
        return None
    keep = [i for i in range(min(QUAD_AXES, int(singular.size))) if singular[i] > 1e-4 * singular[0]]
    if not keep:
        return None
    axes = vt[keep].T
    depth = int(axes.shape[1])
    needed = 1 + depth + depth * (depth + 1) // 2
    if len(rows) + 1 < needed + 2:
        return None
    local = offsets @ axes
    design = [_quadratic_columns(np.zeros(depth))]
    target = [float(score)]
    for i, (_, row_score) in enumerate(rows):
        design.append(_quadratic_columns(local[i]))
        target.append(float(row_score))
    design_m = np.asarray(design, dtype=float)
    target_v = np.asarray(target, dtype=float)
    normal = design_m.T @ design_m
    ridge = 1e-6 * (float(np.trace(normal)) / normal.shape[0])
    penalty = np.eye(normal.shape[0]) * ridge
    penalty[0, 0] = 0.0
    try:
        beta = np.linalg.solve(normal + penalty, design_m.T @ target_v)
    except np.linalg.LinAlgError:
        return None
    slope = beta[1:1 + depth]
    hess = np.zeros((depth, depth))
    cursor = 1 + depth
    for i in range(depth):
        for j in range(i, depth):
            if i == j:
                hess[i, i] = beta[cursor]
            else:
                hess[i, j] = hess[j, i] = 0.5 * beta[cursor]
            cursor += 1
    if np.any(np.linalg.eigvalsh(hess) >= -1e-8):
        return None
    try:
        spot = np.linalg.solve(2.0 * hess, -slope)
    except np.linalg.LinAlgError:
        return None
    if float(np.linalg.norm(spot)) > float(np.max(np.linalg.norm(local, axis=1))) + 1e-9:
        return None
    return np.clip(origin_v + axes @ spot, 0.0, 1.0)


def _snap_encoded(ctx, vec: np.ndarray) -> Optional[Config]:
    """Each coordinate to the legal rung nearest it in the 0–1 encoding."""
    params = {}
    for name, unit in zip(ctx.space.names, vec):
        spec = ctx.space.specs[name]
        rungs = ctx.ladders.rungs.get(name) or []
        if not rungs:
            params[name] = ctx.space._from_unit(spec, float(unit))
            continue
        units = [ctx.space._to_unit(spec, rung) for rung in rungs]
        params[name] = rungs[min(range(len(rungs)), key=lambda i: abs(units[i] - float(unit)))]
    return ctx.space._assemble(params)


def dummy_config(ctx, start: Config, samples: Sequence[Tuple[Sequence[float], float]]) -> Optional[Config]:
    """The legal config at the crest of a local quadratic. None when the bowl does not open downward or the crest
    snaps back onto this sample."""
    origin = _vec(ctx, start)
    peak_score = None
    rows = []
    for vec, score in samples:
        if _dist2(vec, origin) < 1e-18:
            peak_score = float(score) if peak_score is None else max(peak_score, float(score))
            continue
        rows.append((vec, float(score)))
    if peak_score is None:
        return None
    rows.sort(key=lambda row: _dist2(row[0], origin))
    crest = _quadratic_crest(origin, peak_score, rows[:QUAD_NEAR])
    if crest is None:
        return None
    config = _snap_encoded(ctx, crest)
    if config is None or config.key == start.key:
        return None
    return config


def dummies(ctx) -> List[Config]:
    """One dummy per apparent peak, scored configs only."""
    raw = []
    for config, score in ctx.scored():
        if score is None:
            continue
        raw.append((config, float(score), _vec(ctx, config)))
    samples = [(vec, score) for _, score, vec in raw]
    out = []
    seen = {config.key for config, _, _ in raw}
    for peak in apparent_peaks(raw):
        config = dummy_config(ctx, peak, samples)
        if config is not None and config.key not in seen and ctx.fresh(config):
            seen.add(config.key)
            out.append(config)
    return out


def line_point(ctx, higher: Config, lower: Config, t: float) -> Optional[Config]:
    """The legal config at fraction t from the higher centre toward the lower one."""
    params = dict(higher.params)
    for name in ctx.space.names:
        spec = ctx.space.specs[name]
        rungs = ctx.ladders.rungs.get(name) or []
        if not rungs:
            continue
        if spec.allowed or spec.type == "bool":
            params[name] = higher.params.get(name) if t <= 0.5 else lower.params.get(name)
            continue
        a = ctx.ladders.index(name, higher.params.get(name))
        b = ctx.ladders.index(name, lower.params.get(name))
        pos = int(round(a + t * (b - a)))
        pos = max(0, min(len(rungs) - 1, pos))
        params[name] = rungs[pos]
    return ctx.space._assemble(params)


def line_points(ctx, higher: Config, lower: Config) -> List[Config]:
    out, seen = [], set()
    for t in TENTHS:
        config = line_point(ctx, higher, lower, t)
        if config is None or config.key in seen or config.key in (higher.key, lower.key):
            continue
        seen.add(config.key)
        out.append(config)
    return out


def compass_points(ctx, centre: Config, names: Sequence[str]) -> List[Config]:
    """One rung either side of the centre on each named setting. Every other setting stays."""
    out, seen = [], {centre.key}
    for name in names:
        rungs = ctx.ladders.rungs.get(name) or []
        if not rungs:
            continue
        start = int(round(ctx.ladders.index(name, centre.params.get(name))))
        for pos in (start - 1, start + 1):
            if not 0 <= pos < len(rungs):
                continue
            config = single(ctx, centre, name, pos)
            if config is not None and config.key not in seen:
                seen.add(config.key)
                out.append(config)
    return out


def winning_faces(ctx, centre: Config, names: Sequence[str], points: Sequence[Config]) -> Dict[str, object]:
    """Setting to the ±1 value that scores above the centre. A face that does not beat the centre is left out."""
    base = ctx.score(centre)
    if base is None:
        return {}
    faces: Dict[str, Tuple[float, object]] = {}
    for config in points:
        score = ctx.score(config)
        if score is None or float(score) <= float(base):
            continue
        changed = [name for name in names if config.params.get(name) != centre.params.get(name)]
        if len(changed) != 1:
            continue
        name = changed[0]
        prev = faces.get(name)
        if prev is None or float(score) > prev[0]:
            faces[name] = (float(score), config.params.get(name))
    return {name: value for name, (_score, value) in faces.items()}


def cube_points(ctx, centre: Config, faces: Dict[str, object]) -> List[Config]:
    """Joint configs of the live faces. The one-setting faces are already scored."""
    if len(faces) < 2:
        return []
    names = list(faces)
    axes = [(centre.params.get(name), faces[name]) for name in names]
    out, seen = [], {centre.key}
    for combo in itertools.product(*axes):
        changes = dict(zip(names, combo))
        if sum(centre.params.get(name) != changes[name] for name in names) < 2:
            continue
        params = dict(centre.params)
        params.update(changes)
        config = ctx.space._assemble(params)
        if config is not None and config.key not in seen:
            seen.add(config.key)
            out.append(config)
    return out


def vision_child(ctx, centre: Config, parent: Config) -> Optional[Config]:
    """Copy the centre and take one to three settings from the parent."""
    differ = [n for n in ctx.space.names if parent.params.get(n) != centre.params.get(n)]
    if not differ:
        return None
    k = ctx.rng.randint(1, min(3, len(differ)))
    params = dict(centre.params)
    for name in ctx.rng.sample(differ, k):
        params[name] = parent.params.get(name)
    config = ctx.space._assemble(params)
    if config is None or config.key == centre.key:
        return None
    return config


def _between(ctx, a: Config, b: Config) -> Optional[float]:
    if a.key == b.key:
        return 0.0
    dist = getattr(ctx, "dist", None)
    if dist is not None:
        d, _kind = dist.between(a, b)
        return float(d)
    from search.distance import setting_distance
    return setting_distance(ctx.space, a, b)


def _pile_for(ctx, config: Config):
    piles = list(getattr(ctx, "_piles", None) or [])
    for pile in piles:
        if pile.config.key == config.key:
            return pile
        d = _between(ctx, config, pile.config)
        rad = float(getattr(pile, "radius", RADIUS))
        if d is not None and d <= rad:
            return pile
    return None


async def _wait(ctx, configs: Sequence[Config], tag: str, until: float) -> None:
    """Send `configs` first. Every worker they leave free breeds until those runs are back."""
    fresh = [c for c in configs if c is not None and ctx.fresh(c) and ctx.score(c) is None]
    if not fresh:
        return
    futs = ctx.propose(fresh, tag, ["small"] * len(fresh), priority=2)
    await _cover(ctx, futs, until)


async def _cover(ctx, futs, until: float) -> None:
    """Keep breeding on free workers until every future in `futs` is back, or the clock runs out."""
    from search.climb import free_workers

    pending = {f for f in futs if not f.done()}
    while pending and ctx.clock() < until and not ctx.over():
        idle = free_workers(ctx)
        if idle:
            try:
                climber = next((c for c in getattr(ctx, "climbers", ()) if getattr(c, "state", "") == "climbing"), None)
                await fill_free(ctx, climber, idle)
            except Exception:  # noqa: BLE001 - a failed breed must not drop the runs already out
                log.exception("filling the free workers failed")
        _, pending = await asyncio.wait(pending, timeout=2.0, return_when=asyncio.FIRST_COMPLETED)
        pending = {f for f in pending if not f.done()}
        publish = getattr(ctx, "publish", None)
        if callable(publish):
            publish()


async def place_dummies(ctx) -> None:
    """Score one dummy per apparent peak before the climbers sit down."""
    made = dummies(ctx)
    if not made:
        ctx.say("step A: the sample graph's peaks are already on scored configs")
        return
    ctx.say("step A: %d dummy peak(s)", len(made))
    await _wait(ctx, made, "A", ctx.clock() + ctx.T(60 * 15))


def _higher_pair(ctx, a: Config, b: Config) -> Tuple[Config, Config, float, float]:
    sa, sb = float(ctx.score(a)), float(ctx.score(b))
    if sa >= sb:
        return a, b, sa, sb
    return b, a, sb, sa


def _on_path(ctx, end_a: Config, end_b: Config, config: Config) -> Optional[float]:
    """Fraction of the callset path from `end_a` toward `end_b` when `config` lies between them."""
    if config.key in (end_a.key, end_b.key):
        return None
    span = _between(ctx, end_a, end_b)
    if span is None or span <= 1e-12:
        return None
    left = _between(ctx, end_a, config)
    right = _between(ctx, end_b, config)
    if left is None or right is None or left >= span or right >= span:
        return None
    if left + right > span + WOBBLE:
        return None
    return left / span


def _path_samples(ctx, end_a: Config, end_b: Config) -> List[Tuple[Config, float, float]]:
    out = []
    for config, score in ctx.scored():
        if score is None:
            continue
        place = _on_path(ctx, end_a, end_b, config)
        if place is None:
            continue
        out.append((config, float(score), place))
    return out


def _fresh_point(ctx, higher: Config, lower: Config, t: float) -> Optional[Config]:
    config = line_point(ctx, higher, lower, t)
    if config is None or config.key in (higher.key, lower.key) or ctx.score(config) is not None:
        return None
    return config


def valley_status(ctx, a: Config, b: Config) -> Tuple[str, Optional[Config]]:
    """The valley between two fingerprints, and the one config still worth scoring.

    "dip", "slope", "meet", and "undecided" are finished. "need" carries the midpoint,
    or the quarter when the midpoint landed inside the band.
    """
    if ctx.score(a) is None or ctx.score(b) is None:
        return "undecided", None
    higher, lower, _, low_score = _higher_pair(ctx, a, b)
    span = _between(ctx, higher, lower)
    if span is None:
        return "undecided", None
    if span <= GAP:
        return "meet", None
    samples = _path_samples(ctx, higher, lower)
    if any(score <= low_score - VALLEY for _, score, _ in samples):
        return "dip", None
    mid = line_point(ctx, higher, lower, 0.5)
    mid_score = ctx.score(mid) if mid is not None else None
    if mid is not None and mid_score is not None and mid.key not in (higher.key, lower.key):
        if _on_path(ctx, higher, lower, mid) is None:
            return "undecided", None
        if float(mid_score) >= low_score:
            return "slope", None
        quarter = _fresh_point(ctx, higher, lower, 0.25)
        if quarter is not None:
            return "need", quarter
        scored = line_point(ctx, higher, lower, 0.25)
        if scored is None or ctx.score(scored) is None or _on_path(ctx, higher, lower, scored) is None:
            return "undecided", None
        if float(ctx.score(scored)) <= low_score - VALLEY:
            return "dip", None
        return "slope", None
    if any(0.25 <= place <= 0.75 for _, _, place in samples):
        return "slope", None
    probe = _fresh_point(ctx, higher, lower, 0.5)
    if probe is None:
        return "undecided", None
    return "need", probe


async def _resolve(ctx, a: Config, b: Config, tag: str) -> str:
    """Finish one pair. At most the midpoint and the quarter are sent."""
    for _ in range(2):
        decision, probe = valley_status(ctx, a, b)
        if decision != "need" or probe is None:
            return decision
        await _wait(ctx, [probe], tag, ctx.clock() + ctx.T(60 * 15))
        if ctx.score(probe) is None:
            return "undecided"
    decision, _probe = valley_status(ctx, a, b)
    return "undecided" if decision == "need" else decision


def _drop_lower(kept: List[Config], higher: Config, lower: Config) -> List[Config]:
    kept = [config for config in kept if config.key != lower.key]
    if higher.key not in {config.key for config in kept}:
        kept.append(higher)
    return kept


def _apply_pair(ctx, kept: List[Config], challenger: Config, nearest: Config, decision: str) -> List[Config]:
    """Seat `challenger` against `nearest`. A slope or a meeting drops the lower fingerprint."""
    higher, lower, _, _ = _higher_pair(ctx, challenger, nearest)
    if decision == "dip":
        ctx.say("valley between %s and %s", higher.key[:8], lower.key[:8])
        if challenger.key not in {config.key for config in kept}:
            kept.append(challenger)
        return kept
    if decision == "undecided":
        ctx.say("valley %s vs %s is not decided", higher.key[:8], lower.key[:8])
        if challenger.key not in {config.key for config in kept}:
            kept.append(challenger)
        return kept
    if decision == "meet":
        hi_pile, lo_pile = _pile_for(ctx, higher), _pile_for(ctx, lower)
        if hi_pile is not None and lo_pile is not None:
            absorb_parents(ctx, hi_pile, lo_pile)
        ctx.say("met %s; %s's parents join the higher pile", higher.key[:8], lower.key[:8])
        return _drop_lower(kept, higher, lower)
    ctx.say("no valley; %s leaves before the ascent", lower.key[:8])
    return _drop_lower(kept, higher, lower)


def _nearest_seat(ctx, config: Config, seats: Sequence[Config]) -> Optional[Config]:
    best, dist = None, None
    for seat in seats:
        if seat.key == config.key:
            continue
        gap = _between(ctx, config, seat)
        if gap is None:
            continue
        if dist is None or gap < dist:
            best, dist = seat, gap
    return best


async def gate(ctx, seated: List[Config]) -> List[Config]:
    """The one climber sits on the highest fingerprint. The other piles wait."""
    refresh(ctx)
    piles = list(getattr(ctx, "_piles", None) or [])
    order = sorted((pile for pile in piles if ctx.score(pile.config) is not None), key=lambda pile: -float(pile.score))
    budget = min(SEATS, len(seated) if seated else SEATS)
    if not order:
        return list(seated)[:budget]
    kept = [order[0].config]
    ctx.say("valley: %s %.4f is the seat", kept[0].key[:8], float(order[0].score))
    for pile in order[1:]:
        if len(kept) >= budget:
            break
        nearest = _nearest_seat(ctx, pile.config, kept)
        if nearest is None:
            kept.append(pile.config)
            continue
        decision = await _resolve(ctx, pile.config, nearest, "A")
        kept = _apply_pair(ctx, kept, pile.config, nearest, decision)
    return kept[:budget]


async def claim_fingerprint(ctx, climbers) -> bool:
    """A GA fingerprint above the lowest climber takes that seat when the path to its nearest seat dips."""
    refresh(ctx)
    climbing = [cl for cl in climbers if getattr(cl, "state", "") == "climbing" and getattr(cl, "centre", None) is not None]
    if not climbing:
        return False
    lowest = min(climbing, key=lambda cl: (float(cl.score), cl.id))
    piles = list(getattr(ctx, "_piles", None) or [])
    seen = getattr(ctx, "_valley_seen", None)
    if seen is None:
        seen = set()
        ctx._valley_seen = seen
    choice = None
    for i, pile in enumerate(piles):
        if float(pile.score) <= float(lowest.score):
            continue
        if _owned(ctx, piles, i, climbing):
            continue
        nearest = min(climbing, key=lambda cl: (_between(ctx, pile.config, cl.centre) or 9e9, cl.id))
        if nearest.centre is None or nearest.centre.key == pile.config.key:
            continue
        mark = (pile.config.key, nearest.centre.key)
        if mark in seen:
            continue
        if choice is None or float(pile.score) > choice[0]:
            choice = (float(pile.score), i, pile, nearest, mark)
    if choice is None:
        return False
    _score, index, pile, nearest, mark = choice
    decision = await _resolve(ctx, pile.config, nearest.centre, "B")
    seen.add(mark)
    higher, lower, _, _ = _higher_pair(ctx, pile.config, nearest.centre)
    if decision == "dip":
        ctx.say("valley between %s and %s", higher.key[:8], lower.key[:8])
        lowest.origin = pile.label(index)
        lowest.move(pile.config, float(pile.score), "valley")
        return True
    if decision == "meet":
        hi_pile, lo_pile = _pile_for(ctx, higher), _pile_for(ctx, lower)
        if hi_pile is not None and lo_pile is not None:
            absorb_parents(ctx, hi_pile, lo_pile)
        ctx.say("met %s; %s's parents join the higher pile", higher.key[:8], lower.key[:8])
        if higher.key == pile.config.key:
            nearest.origin = pile.label(index)
            nearest.move(pile.config, float(pile.score), "met")
            return True
        return False
    if decision == "slope":
        ctx.say("no valley; %s leaves before the ascent", lower.key[:8])
        if higher.key == pile.config.key:
            nearest.origin = pile.label(index)
            nearest.move(pile.config, float(pile.score), "slope")
            return True
        return False
    ctx.say("valley %s vs %s is not decided", higher.key[:8], lower.key[:8])
    return False


def _children(ctx, centre: Config) -> List[Config]:
    pile = _pile_for(ctx, centre)
    parents = mating_pool(ctx, pile) if pile is not None else []
    parents = [p for p in parents if p.key != centre.key]
    out, seen = [], {centre.key}
    for parent in parents:
        if len(out) >= CHILDREN:
            break
        config = vision_child(ctx, centre, parent)
        if config is not None and config.key not in seen and ctx.fresh(config):
            seen.add(config.key)
            out.append(config)
    return out


def compass_axes(ctx, cl) -> List[str]:
    """Up to four axes. The largest one-setting score change near the centre comes first.

    An axis in ``cl.flat`` is skipped. Axes with no scored one-setting neighbour yet follow, read
    settings before assembly settings.
    """
    centre = cl.centre
    if centre is None:
        return []
    watched = list(getattr(cl, "read_names", ())) + list(getattr(cl, "assembly_names", ()))
    flat = set(getattr(cl, "flat", ()) or ())
    base = float(cl.score)
    gain: Dict[str, float] = {}
    for config, score in _scored_pairs(ctx):
        if score is None:
            continue
        changed = [name for name in watched if config.params.get(name) != centre.params.get(name)]
        if len(changed) != 1:
            continue
        apart = _between(ctx, config, centre)
        if apart is not None and apart > RADIUS:
            continue
        name = changed[0]
        delta = abs(float(score) - base)
        if delta > gain.get(name, -1.0):
            gain[name] = delta
    ranked = sorted((delta, name) for name, delta in gain.items() if name not in flat and delta > 0)
    ranked.sort(reverse=True)
    chosen = [name for _delta, name in ranked[:COMPASS_AXES]]
    # A neighbour that tied the centre is already measured. It does not take a slot
    # from a setting that has no neighbour yet, or the compass shrinks below four.
    fresh, tied = [], []
    for name in watched:
        if name in flat or name in chosen or (name in gain and gain[name] > 0):
            continue
        (tied if name in gain else fresh).append(name)
    for name in fresh + tied:
        if len(chosen) >= COMPASS_AXES:
            break
        chosen.append(name)
    return chosen


def _scored_pairs(ctx) -> List[Tuple[Config, float]]:
    scored = getattr(ctx, "scored", None)
    if not callable(scored):
        return []
    return [(config, score) for config, score in scored() if config is not None and score is not None]


def _apart(ctx, a: Config, b: Config) -> Tuple[Optional[float], str]:
    dist = getattr(ctx, "dist", None)
    if dist is not None and hasattr(dist, "between"):
        d, kind = dist.between(a, b)
        return float(d), str(kind)
    if a.key == b.key:
        return 0.0, "callset"
    from search.distance import setting_distance
    return setting_distance(ctx.space, a, b), "settings"


def _far_from(ctx, config: Config, others: Sequence[Config]) -> bool:
    """Callset distance at least 0.10. With no sketches, the setting-distance floor is enough to leave the pile."""
    from common.tuning import tune
    setting_floor = float(tune("candidates.min_setting_distance"))
    for other in others:
        if other is None or other.key == config.key:
            continue
        apart, kind = _apart(ctx, config, other)
        if apart is None:
            continue
        need = GAP if kind == "callset" else setting_floor
        if apart < need:
            return False
    return True


# Pearson |r| >= 0.20 on the 748 saved configs in search/portfolio/legacy.json
# (the one GVCF config left out). Reproduced to 0.01. The strongest pair takes
# both knobs. A later pair joins a free knob onto that pair. A pair whose both
# knobs are already taken stays as the stronger pairs set it, so a child never
# splits one of these three groups. Every other knob was set on its own.
LINKED_GROUPS = (
    ("min_base_quality_score", "pair_hmm_gap_continuation_penalty",
     "active_probability_threshold", "min_mapping_quality_score", "pcr_indel_model"),
    ("standard_min_confidence_threshold_for_calling",
     "phred_scaled_global_read_mismapping_rate", "base_quality_score_threshold"),
    ("max_assembly_region_size", "assembly_region_padding", "min_assembly_region_size"),
)
_LINKED = {name for group in LINKED_GROUPS for name in group}


def _shift_quarter(ctx, params: dict, name: str) -> None:
    rungs = ctx.ladders.rungs.get(name) or []
    span = max(4, len(rungs) // 4)
    step = ctx.rng.randint(max(2, span // 2), span) * ctx.rng.choice((1, -1))
    value = ctx.ladders.shift(name, params.get(name), step)
    if value is None:
        value = ctx.ladders.shift(name, params.get(name), -step)
    if value is not None:
        params[name] = value


def _far_bases(ctx) -> List[Config]:
    cached = list(getattr(ctx, "_far_bases", None) or [])
    if cached:
        return cached
    top = ctx.best() if callable(getattr(ctx, "best", None)) else None
    return [top[0]] if top is not None else []


def _cross(ctx, a: Config, b: Config) -> Optional[Config]:
    """One child. A linked group comes from one parent. An unlinked knob comes from either."""
    params = dict(a.params)
    taken = set()
    for group in LINKED_GROUPS:
        src = b if ctx.rng.random() < 0.5 else a
        for name in group:
            if name in src.params:
                params[name] = src.params[name]
            taken.add(name)
    for name in ctx.space.names:
        if name in taken:
            continue
        if ctx.rng.random() < 0.5:
            params[name] = b.params.get(name, params[name])
    return ctx.space._assemble(params)


def _pile_groups(ctx, centre: Optional[Config]) -> List[List[Config]]:
    from search.piles import mating_pool
    groups = []
    for pile in list(getattr(ctx, "_piles", None) or []):
        pool = list(mating_pool(ctx, pile) or [pile.config])
        pool = [parent for parent in pool if centre is None or parent.key != centre.key]
        seen = set()
        unique = []
        for parent in pool:
            if parent.key in seen:
                continue
            seen.add(parent.key)
            unique.append(parent)
        if unique:
            groups.append(unique)
    return groups


def hill_children(ctx, cl, n: int) -> List[Config]:
    """A cross of two piles, at least 0.10 from the climber and from children in flight.

    Linked knobs stay with one parent. With only one pile, four to six unlinked
    settings move by about a quarter of the ladder. A linked knob stays put.
    """
    if n <= 0:
        return []
    centre = getattr(cl, "centre", None) if cl is not None else None
    bases = [base for base in _far_bases(ctx) if centre is None or base.key != centre.key]
    groups = _pile_groups(ctx, centre)
    names = [name for name in ctx.space.names if ctx.ladders.movable(name)]
    if len(names) < 4 and len(groups) < 2:
        return []
    blocked = [centre] if centre is not None else []
    for pile in list(getattr(ctx, "_piles", None) or []):
        blocked.append(pile.config)
    inflight = [config for config in list(getattr(ctx, "_hill_live", None) or []) if ctx.score(config) is None]
    blocked.extend(inflight)
    out, seen = [], set()
    tries = 0
    while len(out) < n and tries < 30 * n:
        tries += 1
        if len(groups) >= 2:
            left, right = ctx.rng.sample(groups, 2)
            config = _cross(ctx, ctx.rng.choice(left), ctx.rng.choice(right))
        elif bases and len(names) >= 4:
            params = dict(bases[tries % len(bases)].params)
            free = [name for name in names if name not in _LINKED]
            if len(free) < 4:
                break
            k = ctx.rng.randint(4, min(6, len(free)))
            for name in ctx.rng.sample(free, k):
                _shift_quarter(ctx, params, name)
            config = ctx.space._assemble(params)
        else:
            break
        if config is None or config.key in seen or not ctx.fresh(config) or ctx.score(config) is not None:
            continue
        if centre is not None and sum(config.params.get(name) != centre.params.get(name) for name in names) < 4:
            continue
        if not _far_from(ctx, config, blocked):
            continue
        seen.add(config.key)
        blocked.append(config)
        out.append(config)
    return out


def population_children(ctx, n: int) -> List[Config]:
    """A cross of scored configs inside Step A's range. A linked group stays with one parent.

    Unlike the hill lane, the child may land on the climber's own callset. A higher score there
    becomes the next centre once the open compass is back. With fewer than two configs in the
    range, the cross uses the scored population.
    """
    if n <= 0:
        return []
    if range_open():
        area = [config for config, _score in candidate_range(ctx)]
    else:
        area = []
    if len(area) >= 2:
        ranked = [(config, 0.0) for config in area]
    else:
        ranked = sorted(_scored_pairs(ctx), key=lambda item: (-float(item[1]), item[0].key))
    parents: List[Config] = []
    seen_parents = set()
    for config, _score in ranked:
        if config.key in seen_parents:
            continue
        seen_parents.add(config.key)
        parents.append(config)
        if len(parents) >= 64:
            break
    if len(parents) < 2:
        return []
    out, seen = [], set()
    tries = 0
    while len(out) < n and tries < 30 * n:
        tries += 1
        left, right = ctx.rng.sample(parents, 2)
        config = _cross(ctx, left, right)
        if config is None or config.key in seen or not ctx.fresh(config) or ctx.score(config) is not None:
            continue
        seen.add(config.key)
        out.append(config)
    return out


def climb_probes(ctx, cl, n: int) -> List[Config]:
    """Copies of the centre with one to three settings the compass has not tried. They do not move the centre."""
    if cl is None or n <= 0 or getattr(cl, "centre", None) is None:
        return []
    tried = set(getattr(cl, "tried_axes", ()) or ())
    flat = set(getattr(cl, "flat", ()) or ())
    names = [name for name in list(cl.read_names) + list(cl.assembly_names) if name not in tried and name not in flat]
    if not names:
        return []
    out, seen = [], {cl.centre.key}
    tries = 0
    while len(out) < n and tries < 20 * n:
        tries += 1
        params = dict(cl.centre.params)
        k = ctx.rng.randint(1, min(3, len(names)))
        for name in ctx.rng.sample(names, k):
            step = ctx.rng.choice((1, -1))
            value = ctx.ladders.shift(name, params.get(name), step)
            if value is None:
                value = ctx.ladders.shift(name, params.get(name), -step)
            if value is not None:
                params[name] = value
        config = ctx.space._assemble(params)
        if config is None or config.key in seen or not ctx.fresh(config) or ctx.score(config) is not None:
            continue
        if config.key == cl.centre.key:
            continue
        seen.add(config.key)
        out.append(config)
    return out


def _send_lane(ctx, configs: Sequence[Config], why: str) -> None:
    if not configs or not hasattr(ctx, "propose"):
        return
    live = list(getattr(ctx, "_hill_live", None) or [])
    if why == "hill lane":
        live.extend(configs)
        ctx._hill_live = live
    futs = ctx.propose(list(configs), "G", ["big"] * len(configs), priority=0) or []
    if why == "step A":
        ctx._range_futs = [f for f in list(getattr(ctx, "_range_futs", None) or []) + list(futs)
                           if f is not None and not f.done()]
    ctx.say("%s: %d", why, len(configs))


def range_open() -> bool:
    """Step A keeps a score range and updates it while Step B pinpoints."""
    return VERSION.startswith("v9.34")


def candidate_range(ctx) -> List[Tuple[Config, float]]:
    """Scored configs close to the peak. Step A presents this area. Step B picks a point inside it."""
    pairs = [(config, float(score)) for config, score in _scored_pairs(ctx)]
    if not pairs:
        return []
    best = max(score for _config, score in pairs)
    close = [(config, score) for config, score in pairs if best - score <= RANGE_GAP]
    close.sort(key=lambda item: (-item[1], item[0].key))
    out, seen = [], set()
    for config, score in close:
        if config.key in seen:
            continue
        seen.add(config.key)
        out.append((config, score))
        if len(out) >= 32:
            break
    return out


def note_range(ctx) -> None:
    """Record the range and say so when the peak or its size changes."""
    rows = candidate_range(ctx)
    ctx.candidate_range = rows
    if not rows:
        return
    mark = (rows[0][0].key, len(rows), round(rows[0][1], 4))
    if mark == getattr(ctx, "_range_mark", None):
        return
    ctx._range_mark = mark
    ctx.say("step A: range of %d candidates, %.4f down to %.4f", len(rows), rows[0][1], rows[-1][1])


def refine_range(ctx, n: int) -> List[Config]:
    """A few points that move the range. One batch stays out at a time, so Step B keeps the fleet."""
    if not getattr(ctx, "_step_a_open", False) or n <= 0:
        return []
    pending = [f for f in getattr(ctx, "_range_futs", None) or [] if f is not None and not f.done()]
    ctx._range_futs = pending
    if pending:
        note_range(ctx)
        return []
    budget = min(RANGE_BATCH, n)
    take: List[Config] = []
    if not getattr(ctx, "_dummies_sent", False):
        ctx._dummies_sent = True
        take = [config for config in dummies(ctx)[:budget] if ctx.fresh(config) and ctx.score(config) is None]
    rest = list(getattr(ctx, "_step_a_rest", None) or [])
    while rest and len(take) < budget:
        config = rest.pop(0)
        if ctx.fresh(config) and ctx.score(config) is None:
            take.append(config)
    ctx._step_a_rest = rest
    if take:
        ctx.say("step A: range update, %d points", len(take))
    note_range(ctx)
    return take


def mating_room(n: int) -> int:
    """Free cores after a compass-sized set is left idle. A small test fleet is used whole."""
    if n <= COMPASS_SLOTS:
        return max(0, n)
    return n - COMPASS_SLOTS


def pile_ga(ctx, cl, n: int) -> List[Config]:
    """The v9.20 leftover GA, beside the compass.

    Half the children are a cross inside the piles, and that child may land on the climber.
    Half move six unlinked settings and are sent only at least 0.10 from every pile.
    A linked group stays with one parent. The compass is not held for either half.
    """
    if n <= 0:
        return []
    from search.piles import mating_pool
    centre = getattr(cl, "centre", None) if cl is not None else None
    piles = list(getattr(ctx, "_piles", None) or [])
    names = [name for name in ctx.space.names if ctx.ladders.movable(name)]
    free = [name for name in names if name not in _LINKED]
    inside_n = n - n // 2
    out, seen = [], set()
    tries = 0
    while len(out) < inside_n and tries < 20 * max(1, inside_n) and piles:
        tries += 1
        pile = piles[len(out) % len(piles)]
        pool = list(mating_pool(ctx, pile) or [pile.config])
        if len(pool) >= 2:
            left = ctx.rng.choice(pool)
            right = ctx.rng.choice([parent for parent in pool if parent.key != left.key] or pool)
            config = _cross(ctx, left, right)
        elif free:
            params = dict(pile.config.params)
            for name in ctx.rng.sample(free, min(4, len(free))):
                _shift_quarter(ctx, params, name)
            config = ctx.space._assemble(params)
        else:
            break
        if config is None or config.key in seen or not ctx.fresh(config) or ctx.score(config) is not None:
            continue
        seen.add(config.key)
        out.append(config)
    far_n = n // 2
    bases = [pile.config for pile in piles] or ([centre] if centre is not None else [])
    blocked = list(bases)
    tries = 0
    while len(out) < inside_n + far_n and tries < 30 * max(1, far_n) and bases and len(free) >= 4:
        tries += 1
        params = dict(bases[tries % len(bases)].params)
        for name in ctx.rng.sample(free, min(6, len(free))):
            _shift_quarter(ctx, params, name)
        config = ctx.space._assemble(params)
        if config is None or config.key in seen or not ctx.fresh(config) or ctx.score(config) is not None:
            continue
        if not _far_from(ctx, config, blocked):
            continue
        seen.add(config.key)
        blocked.append(config)
        out.append(config)
    return out


async def fill_free(ctx, cl, n: int) -> None:
    """The compass keeps its cores. The rest run the pile GA.

    Sixteen cores stay idle once the fleet is larger than that, so the next compass
    starts on its own heads. Step C sends climb probes only.
    """
    if n <= 0 or ctx.over() or cl is None:
        return
    room = mating_room(n)
    if (getattr(cl, "origin", "") or "") == "step C":
        _send_lane(ctx, climb_probes(ctx, cl, room), "climb probe")
        return
    area = refine_range(ctx, room) if range_open() else []
    _send_lane(ctx, area, "step A")
    _send_lane(ctx, pile_ga(ctx, cl, max(0, room - len(area))), "GA")


async def opening_scouts(ctx) -> None:
    """Two short waves from the far parents. A scout that would land on an existing pile is not sent."""
    step = getattr(ctx, "step_a", None)
    groups = getattr(step, "groups", None)
    bases: List[Config] = []
    if callable(groups):
        built = groups()
        if hasattr(built, "__await__"):
            built = await built
        for name in ("core", "fp"):
            pts = (built or {}).get(name) or []
            if pts:
                config = pts[0].config
                if all(config.key != base.key for base in bases):
                    bases.append(config)
    if not bases:
        bases = _far_bases(ctx)
    ctx._far_bases = bases
    refresh(ctx)
    for wave in range(1, SCOUT_WAVES + 1):
        if ctx.over():
            return
        kids = hill_children(ctx, None, SCOUT_N)
        if not kids:
            ctx.say("step A: the far parents have no scout 0.10 from the piles")
            return
        ctx.say("step A: far scout wave %d, %d points", wave, len(kids))
        if step is not None and hasattr(step, "wave"):
            await step.wave(kids, 0, f"far scouts {wave}", ctx.clock() + ctx.T(8 * 60))
        else:
            await _wait(ctx, kids, "A", ctx.clock() + ctx.T(8 * 60))
        refresh(ctx)


async def aim_compass(ctx, cl, valley: bool = True) -> None:
    """Point the next compass at the highest score. A new fingerprint is valley-tested first."""
    if getattr(cl, "state", "") != "climbing" or getattr(cl, "centre", None) is None:
        return
    if valley and (getattr(cl, "origin", "") or "") != "step C":
        await claim_fingerprint(ctx, [cl])
    if cl.state != "climbing" or cl.centre is None:
        return
    pile = _pile_for(ctx, cl.centre)
    rad = float(getattr(pile, "radius", RADIUS)) if pile is not None else RADIUS
    best = None
    for config, score in _scored_pairs(ctx):
        if config.key == cl.centre.key or float(score) <= float(cl.score):
            continue
        apart = _between(ctx, config, cl.centre)
        if apart is None or apart > rad:
            continue
        if best is None or float(score) > best[1]:
            best = (config, float(score))
    if best is None:
        if range_open():
            note_range(ctx)
        return
    cl.move(best[0], best[1], "best")
    ctx.say("climber %d compasses the highest %s %.4f", cl.id, best[0].key[:8], best[1])
    if range_open():
        note_range(ctx)

    if not ready:
        return
    config, score, origin = list(ready)[:SEATS][0]
    cl = Climber(ctx, len(ctx.climbers), config, score, origin=origin)
    cl.climb_until = until
    cl.compass_only = True
    ctx.climbers.append(cl)
    await aim_compass(ctx, cl)
    ctx.say("step B: compass on %s %.4f until the fine-tune", cl.centre.key[:8], cl.score)
    await cl.run(until, jumps=False)


def best_grid_point(ctx, points: Sequence[Config], centre_score: float):
    """The best scored point of a finished grid, when it beats the centre. A child that opened the view is not a point."""
    best = None
    for config in points:
        score = ctx.score(config)
        if score is None:
            continue
        if best is None or float(score) > best[1]:
            best = (config, float(score))
    if best is None or best[1] <= float(centre_score):
        return None
    return best


def _note_first(cl) -> None:
    """The basin row's first key is the centre after this climber's first view."""
    if cl.first is None:
        cl.first = (float(cl.score), cl.centre.key)
        cl.iteration += 1


def _keep_move(ctx, climber, config: Config, score: float, climbers) -> str:
    """"move" when the caller should adopt the point. "handled" when a meeting already did. "stop" when the view is done."""
    others = [cl for cl in climbers if cl is not climber and cl.state == "climbing" and cl.centre is not None]
    d_self = _between(ctx, config, climber.centre)
    pile = _pile_for(ctx, climber.centre)
    rad = float(getattr(pile, "radius", RADIUS)) if pile is not None else RADIUS
    for other in others:
        d = _between(ctx, config, other.centre)
        if d is None or d > GAP:
            continue
        hi, lo = (climber, other) if score >= other.score else (other, climber)
        if score > hi.score:
            hi.move(config, score, "met")
        hi_pile, lo_pile = _pile_for(ctx, hi.centre), _pile_for(ctx, lo.centre)
        if hi_pile is not None and lo_pile is not None:
            absorb_parents(ctx, hi_pile, lo_pile)
        lo.state = "stopped"
        ctx.say("climber %d met climber %d; the lower parents join the higher pile", hi.id, lo.id)
        return "handled"
    if d_self is None or d_self > rad:
        if all((_between(ctx, config, other.centre) or 1.0) > GAP for other in others):
            best = ctx.best() if callable(getattr(ctx, "best", None)) else None
            if best is not None and best[1] - score <= FLOOR:
                ctx.say("climber %d sees %s %.4f on another hill", climber.id, config.key[:8], score)
        return "stop"
    for other in others:
        apart = _between(ctx, climber.centre, other.centre)
        if apart is None or apart <= GAP:
            continue
        if score >= other.score:
            hi_cfg, lo_cfg, lo_score = config, other.centre, other.score
        else:
            hi_cfg, lo_cfg, lo_score = other.centre, config, score
        dips = [float(s) for point in line_points(ctx, hi_cfg, lo_cfg)
                if (s := ctx.score(point)) is not None]
        if not dips or any(s <= lo_score - VALLEY for s in dips):
            continue
        if score < other.score:
            climber.state = "stopped"
            ctx.say("climber %d has no valley against climber %d and leaves", climber.id, other.id)
            return "handled"
        other.state = "stopped"
        ctx.say("climber %d has no valley against climber %d; the lower one leaves", climber.id, other.id)
    return "move"
