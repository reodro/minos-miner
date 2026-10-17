"""Centre jumps (leaps): one expensive setting moved two or more rungs from a climber's centre.

The compass moves one rung. A jump moves the read-pass, assembly or activity setting further; every rung of the
cheap settings already rides on the centre's own run. 044200Z (2026-10-16): the 17:06 run certified 82.6168 at
min_mapping_quality_score 59. 58 and 60 scored lower and the line went straight to 0, but 46 was +0.0785.

v2 (the climbers use this): a bandit over (setting, distance class 2-3, 4-8, 9-21, 22+ rungs). Each jump comes from the
arm with the highest Thompson sample of win rate x gain per win. The priors come from four runs of this round (one-rung
peaks at 82.3+: further out on mapping quality and haplotypes won 5.6% of runs at 2 rungs, 11-13% at 4-21; smooth
settings' peaks held), and every jump's result updates its arm for the rest of the round. Inside a class, half the
picks take the rung with the best expected value (l / g over the top 15% of scores, as in TPE), half a random rung.
Before a peak is named: one jump per (setting, class, direction), then two-setting jumps on the leading arms.

v3: each arm is also a direction (mapping quality only ever won going down, max assembly region going up), a winning
jump is tried again first from the new centre and then twice as far, and the islands share one set of arms.

v1 (kept for the probe in search/probe.py): fixed sizes 2, 3, 5, 8, 13, 21, 34, half nearest first and half by
expected value x nearness. On 390 scanned lines that order needed 3.9 runs to reach a gain (nearest first 4.6, random
10.2), but from one-rung peaks the nearest jumps are the ones that win least.
"""
from __future__ import annotations

import itertools
import math
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from common.types import Config

LEAP_RUNGS = (2, 3, 5, 8, 13, 21, 34)
TOP_SHARE = 0.15            # the "good" scores of l / g
MIN_TOP = 10


def candidates(ctx, centre: Config, names: Sequence[str],
               skip: Callable[[Config], bool]) -> List[Tuple[str, int, Config]]:
    """Every leap from `centre` still worth a run, nearest first: (setting, signed rungs, config).

    A rung that assembles into the centre's own run, a config already scored and one `skip` refuses are left out."""
    out: List[Tuple[str, int, Config]] = []
    seen = {centre.key}
    ladders = ctx.ladders
    for k in LEAP_RUNGS:
        for name in names:
            if not ladders.movable(name) or k >= len(ladders.rungs.get(name, ())):
                continue
            for sign in (1, -1):
                value = ladders.shift(name, centre.params.get(name), sign * k)
                if value is None or value == centre.params.get(name):
                    continue
                config = ctx.space._assemble({**centre.params, name: value})
                if config is None or config.key in seen:
                    continue
                seen.add(config.key)
                if ctx.score(config) is not None or skip(config):
                    continue
                out.append((name, sign * k, config))
    return out


REFRESH_S = 60.0          # the expected values are rebuilt at most this often (search clock)


def _rung(ctx, name: str, value) -> Optional[int]:
    try:
        return int(round(ctx.ladders.index(name, value)))
    except Exception:  # noqa: BLE001 - a value the ladder cannot place has no expected value
        return None


def expected(ctx, names: Sequence[str]) -> Dict[str, List[float]]:
    """setting -> l / g for each ladder rung: the top 15% of this search's scores against the rest, smoothed over
    neighbouring rungs (a Gaussian of max(1, rungs / 20) rungs). Empty until there are 20 scores. Kept on the
    context for REFRESH_S: every climber's every wave reads it, and a search has tens of thousands of scores."""
    clock = getattr(ctx, "clock", None)
    now = float(clock()) if callable(clock) else 0.0
    key = tuple(names)
    hit = getattr(ctx, "_leap_values", None)
    if hit is not None and hit[0] == key and now - hit[1] < REFRESH_S:
        return hit[2]
    out = _expected(ctx, names)
    try:
        ctx._leap_values = (key, now, out)
    except AttributeError:
        pass
    return out


def _expected(ctx, names: Sequence[str]) -> Dict[str, List[float]]:
    rows = [(s, c) for c, s in ctx.scored()]
    if len(rows) < 2 * MIN_TOP:
        return {}
    rows.sort(key=lambda t: (-t[0], t[1].key))
    cut = max(MIN_TOP, int(TOP_SHARE * len(rows)))
    top, rest = rows[:cut], rows[cut:]
    out: Dict[str, List[float]] = {}
    for name in names:
        rungs = ctx.ladders.rungs.get(name) or ()
        n = len(rungs)
        if n < 2:
            continue
        at: Dict[object, Optional[int]] = {}
        for i, v in enumerate(rungs):
            try:
                at.setdefault(v, i)
            except TypeError:
                pass
        good, other = [0.0] * n, [0.0] * n
        for group, counts in ((top, good), (rest, other)):
            for _, c in group:
                v = c.params.get(name)
                try:
                    i = at[v]
                except KeyError:
                    i = at[v] = _rung(ctx, name, v)          # off the ladder (a hypercube float): placed once
                except TypeError:
                    i = _rung(ctx, name, v)
                if i is not None and 0 <= i < n:
                    counts[i] += 1.0
        h = max(1.0, n / 20.0)
        span = int(3 * h) + 1
        kern = [math.exp(-0.5 * (d / h) ** 2) for d in range(span + 1)]
        tg, to = sum(good) or 1.0, sum(other) or 1.0
        ratio = []
        for i in range(n):
            lo, hi = max(0, i - span), min(n, i + span + 1)
            lg = sum(good[j] * kern[abs(i - j)] for j in range(lo, hi)) / tg
            lr = sum(other[j] * kern[abs(i - j)] for j in range(lo, hi)) / to
            ratio.append((lg + 1e-9) / (lr + 1e-6))
        out[name] = ratio
    return out


def order(ctx, centre: Config, names: Sequence[str], skip: Callable[[Config], bool],
          ratio: Optional[Dict[str, List[float]]] = None) -> List[Config]:
    """The leaps to send, interleaved: the nearest untried one, then the best by expected value x nearness, and so on."""
    found = candidates(ctx, centre, names, skip)
    if not found:
        return []
    ratio = expected(ctx, names) if ratio is None else ratio

    def value(item: Tuple[str, int, Config]) -> float:
        name, k, config = item
        r = ratio.get(name)
        i = _rung(ctx, name, config.params.get(name)) if r else None
        if i is None or not 0 <= i < len(r):
            return 0.0
        return r[i] / (1.0 + abs(k))

    by_value = sorted(found, key=lambda item: -value(item))         # stable: equal values stay nearest first
    out: List[Config] = []
    keys = set()
    for a, b in zip(found, by_value):
        for _, _, config in (a, b):
            if config.key not in keys:
                keys.add(config.key)
                out.append(config)
    return out


# ---- v2: a bandit over (setting, distance class) --------------------------------------------------------------------
# 044200Z, four runs, from every one-rung peak on a line at 82.3 or more (both neighbours no better): a value further
# out beat it in 65-96% of the peaks on the rugged settings (mapping quality, haplotypes; and on the cheap side the
# calling threshold, gap penalty, base-quality threshold, mismapping rate) and in 0-14% on the smooth ones. Per run
# sent, on mapping quality and haplotypes: 2 rungs out won 5.6%, 3: 10.1%, 4-8: 11-13%, 9-21: 13%, 22+: 5.6%; the gain
# per run at 14-21 rungs was 17 times that at 2. So the arms below start there and learn from this round's jumps:
# Thompson sampling of win rate x gain per win (TuRBO's batch sampling, Eriksson et al. 2019; operator credit from
# the operator's own results, Fialho et al. 2010), one value per pick, sizes random inside a class (Doerr et al. 2017).
#
# v3 (2026-10-16, from the 17 jumps that moved a centre in the first two runs with jumps): an arm is a direction too.
# Mapping quality won 7 times, always down; max assembly region 4 times, always up; haplotypes 3 times, down. A winning
# move is tried again first from the new centre, then twice as far (MADS: the last successful direction first, then a
# longer step); max assembly region +3 won three times in a row that way. The islands share their arms: which way a
# setting pays did not depend on the region. 2-3 rungs won 6 of the 17, so its rugged prior is 10%, not 7%.

CLASSES = ((2, 3), (4, 8), (9, 21), (22, 10 ** 6))
# v3: the two region sizes count as rugged. The line survey called them smooth on 39 and 56 one-rung peaks, but in the
# first two runs with jumps they gave 7 of the 19 winning jumps and +0.318 of the +0.607 those jumps gained.
RUGGED = ("min_mapping_quality_score", "max_num_haplotypes_in_population", "max_assembly_region_size",
          "min_assembly_region_size")
SMOOTH = ("min_base_quality_score", "min_pruning", "assembly_region_padding", "active_probability_threshold")
WIN_PRIOR = {"rugged": (0.10, 0.12, 0.13, 0.06), "smooth": (0.02, 0.04, 0.04, 0.04), "other": (0.04, 0.08, 0.08, 0.05)}
GAIN_PRIOR = (0.017, 0.04, 0.09, 0.10)          # gain per win, by class
PRIOR_N = 10.0                                  # the win-rate prior is worth this many jumps
GAIN_N = 2.0                                    # the gain prior is worth this many wins


def kind(name: str) -> str:
    return "rugged" if name in RUGGED else ("smooth" if name in SMOOTH else "other")


def size_class(k: int) -> Optional[int]:
    k = abs(int(k))
    for i, (lo, hi) in enumerate(CLASSES):
        if lo <= k <= hi:
            return i
    return None


class Arms:
    """Win rate and gain per win of every (setting, class, direction) arm, shared by the islands for the round."""

    def __init__(self, carried: Optional[Dict[tuple, List[float]]] = None):
        self.stats: Dict[tuple, List[float]] = {}                     # arm -> [wins, tries, gain of the wins]
        self.carried: Dict[tuple, List[float]] = dict(carried or {})  # earlier rounds' counts (load_arms)

    def _post(self, arm: tuple) -> Tuple[float, float, float]:
        name, ci = arm[0], arm[1]
        wins, tries, gains = self.stats.get(arm, (0.0, 0.0, 0.0))
        cw, ct, cg = self.carried.get(arm, (0.0, 0.0, 0.0))
        wins, tries, gains = wins + cw, tries + ct, gains + cg
        p0 = WIN_PRIOR[kind(name)][ci]
        a = p0 * PRIOR_N + wins
        b = (1.0 - p0) * PRIOR_N + max(0.0, tries - wins)
        g = (GAIN_PRIOR[ci] * GAIN_N + gains) / (GAIN_N + wins)
        return a, b, g

    def sample(self, rng, arm: tuple) -> float:
        a, b, g = self._post(arm)
        return rng.betavariate(a, b) * g

    def mean(self, arm: tuple) -> float:
        a, b, g = self._post(arm)
        return a / (a + b) * g

    def learn(self, arm: tuple, gain: float, accept: float) -> bool:
        s = self.stats.setdefault(arm, [0.0, 0.0, 0.0])
        won = gain >= accept
        s[1] += 1.0
        if won:
            s[0] += 1.0
            s[2] += gain
        return won


def arms_of(ctx) -> Arms:
    arms = getattr(ctx, "_jump_arms", None)
    if arms is None:
        arms = Arms()
        try:
            ctx._jump_arms = arms
        except AttributeError:
            pass
    return arms


def _values(ctx, centre: Config, name: str, ci: int, sign: Optional[int] = None) -> List[Tuple[int, object]]:
    """(signed rungs, value) of every rung in the class from the centre, both ways unless `sign` is given."""
    lo, hi = CLASSES[ci]
    n = len(ctx.ladders.rungs.get(name) or ())
    out = []
    for k in range(lo, min(hi, n) + 1):
        for s in ((1, -1) if sign is None else (sign,)):
            v = ctx.ladders.shift(name, centre.params.get(name), s * k)
            if v is not None and v != centre.params.get(name):
                out.append((s * k, v))
    return out


def draw(ctx, centre: Config, name: str, ci: int, seen: set, skip: Callable[[Config], bool], ratio, rng,
         sign: Optional[int] = None) -> Optional[Tuple[str, int, int, Config]]:
    """One untried jump of the arm: half the time its best rung by expected value (l / g), half at random."""
    options = []
    for k, v in _values(ctx, centre, name, ci, sign):
        config = ctx.space._assemble({**centre.params, name: v})
        if config is None or config.key in seen:
            continue
        seen.add(config.key)
        if ctx.score(config) is not None or skip(config):
            continue
        options.append((k, config))
    if not options:
        return None
    r = (ratio or {}).get(name)
    if r and rng.random() < 0.5:
        def value(item):
            i = _rung(ctx, name, item[1].params.get(name))
            return r[i] if i is not None and 0 <= i < len(r) else 0.0
        top = max(value(o) for o in options)
        options = [o for o in options if value(o) >= top]
    k, config = options[rng.randrange(len(options))]
    for _k, other in options:                                       # the ones not drawn stay available
        if other.key != config.key:
            seen.discard(other.key)
    return name, ci, k, config


def pick(ctx, centre: Config, names: Sequence[str], skip: Callable[[Config], bool], n: int, rng,
         taken=()) -> List[Tuple[str, int, int, Config]]:
    """Up to n jumps from `centre`, each from the (setting, class, direction) arm with the highest Thompson sample:
    (setting, class, signed rungs, config)."""
    arms = arms_of(ctx)
    ratio = expected(ctx, names)
    live = [(name, ci, sign) for name in names if ctx.ladders.movable(name)
            for ci in range(len(CLASSES)) if CLASSES[ci][0] < len(ctx.ladders.rungs.get(name) or ())
            for sign in (1, -1)]
    out: List[Tuple[str, int, int, Config]] = []
    seen = set(taken) | {centre.key}
    while len(out) < n and live:
        arm = max(live, key=lambda a: arms.sample(rng, a))
        got = draw(ctx, centre, arm[0], arm[1], seen, skip, ratio, rng, sign=arm[2])
        if got is None:
            live.remove(arm)
            continue
        seen.add(got[3].key)
        out.append(got)
    return out


def coverage(ctx, centre: Config, names: Sequence[str], skip: Callable[[Config], bool], rng,
             taken=()) -> List[Tuple[str, int, int, Config]]:
    """The certificate: one untried jump for every (setting, class, direction) that still has one."""
    ratio = expected(ctx, names)
    out = []
    seen = set(taken) | {centre.key}
    for name in names:
        if not ctx.ladders.movable(name):
            continue
        for ci in range(len(CLASSES)):
            for sign in (1, -1):
                got = draw(ctx, centre, name, ci, seen, skip, ratio, rng, sign=sign)
                if got is not None:
                    seen.add(got[3].key)
                    out.append(got)
    return out


def pairs(ctx, centre: Config, names: Sequence[str], skip: Callable[[Config], bool], rng, top: int = 3,
          taken=()) -> List[Tuple[str, int, int, Config]]:
    """After the certificate: two settings at once, the `top` settings by posterior gain, each at its best class and
    both directions (stagnation detection raises the strength once one-setting moves are spent: Rajabi & Witt 2020)."""
    arms = arms_of(ctx)
    best = {}
    for name in names:
        if not ctx.ladders.movable(name):
            continue
        n = len(ctx.ladders.rungs.get(name) or ())
        cells = [(ci, sg) for ci in range(len(CLASSES)) if CLASSES[ci][0] < n for sg in (1, -1)]
        if cells:
            ci, sg = max(cells, key=lambda cell: arms.mean((name, cell[0], cell[1])))
            best[name] = (arms.mean((name, ci, sg)), ci)
    lead = sorted(best, key=lambda m: -best[m][0])[:top]
    out = []
    seen = set(taken) | {centre.key}
    for a_i in range(len(lead)):
        for b_i in range(a_i + 1, len(lead)):
            a, b = lead[a_i], lead[b_i]
            for sa in (1, -1):
                for sb in (1, -1):
                    va = _values(ctx, centre, a, best[a][1], sa)
                    vb = _values(ctx, centre, b, best[b][1], sb)
                    if not va or not vb:
                        continue
                    ka, xa = va[rng.randrange(len(va))]
                    kb, xb = vb[rng.randrange(len(vb))]
                    config = ctx.space._assemble({**centre.params, a: xa, b: xb})
                    if config is None or config.key in seen:
                        continue
                    seen.add(config.key)
                    if ctx.score(config) is not None or skip(config):
                        continue
                    out.append((a + "+" + b, best[a][1], ka, config))
    return out


def momentum(ctx, centre: Config, name: str, sign: int, k: int, skip: Callable[[Config], bool],
             taken=()) -> List[Tuple[str, int, int, Config]]:
    """After a winning jump: the same move again from the new centre, then twice as far (MADS polls the last
    successful direction first and lengthens the step after a success; Audet & Dennis 2006)."""
    out: List[Tuple[str, int, int, Config]] = []
    seen = set(taken) | {centre.key}
    for kk in (abs(int(k)), 2 * abs(int(k))):
        ci = size_class(kk)
        if ci is None:
            continue
        v = ctx.ladders.shift(name, centre.params.get(name), sign * kk)
        if v is None or v == centre.params.get(name):
            continue
        config = ctx.space._assemble({**centre.params, name: v})
        if config is None or config.key in seen:
            continue
        seen.add(config.key)
        if ctx.score(config) is not None or skip(config):
            continue
        out.append((name, ci, sign * kk, config))
    return out


def mutations(ctx, centre: Config, names: Sequence[str], skip: Callable[[Config], bool], n: int, rng,
              taken=(), sizes=(2, 3)) -> List[Tuple[str, int, int, Config]]:
    """Up to n configs that move two or three settings of `centre` at once, each by a value drawn from its own arm,
    the arms taken in order of a Thompson sample (an elite-anchored GA child: the centre is the only parent). Used to
    fill Step C's idle workers. The setting name is the moved settings joined by '+', which is not an arm."""
    arms = arms_of(ctx)
    live = [(name, ci, sign) for name in names if ctx.ladders.movable(name)
            for ci in range(len(CLASSES)) if CLASSES[ci][0] < len(ctx.ladders.rungs.get(name) or ())
            for sign in (1, -1)]
    out: List[Tuple[str, int, int, Config]] = []
    seen = set(taken) | {centre.key}
    tries = 0
    while len(out) < n and tries < 20 * max(1, n) and live:
        tries += 1
        want = rng.choice(sizes)
        chosen, used = [], set()
        for arm in sorted(live, key=lambda a: -arms.sample(rng, a)):
            if arm[0] in used:
                continue
            values = _values(ctx, centre, arm[0], arm[1], arm[2])
            if not values:
                continue
            k, v = values[rng.randrange(len(values))]
            chosen.append((arm, k, v))
            used.add(arm[0])
            if len(chosen) == want:
                break
        if len(chosen) < 2:
            break
        config = ctx.space._assemble({**centre.params, **{arm[0]: v for arm, _k, v in chosen}})
        if config is None or config.key in seen:
            continue
        seen.add(config.key)
        if ctx.score(config) is not None or skip(config):
            continue
        out.append(("+".join(arm[0] for arm, _k, _v in chosen), chosen[0][0][1], chosen[0][1], config))
    return out


def relink(ctx, centre: Config, score: float, skip: Callable[[Config], bool], rng, taken=(), elsewhere=(),
           targets: int = 3, per_size=(10, 10, 6), below: float = 0.10) -> List[Tuple[str, int, int, Config]]:
    """Path relinking (Glover, Laguna & Marti 2000): from a one-setting peak toward strong configs elsewhere. Targets are
    the highest scores within `below` of the centre that differ from it in three or more settings (this search's, plus
    `elsewhere`: (config, score) from the other islands); each yields the centre with 2, 3 and 4 of the target's values
    (per_size of each). 044200Z islands (2026-10-16): every single move from island 1's 82.6298 was lower, while v1's
    82.6664 in the same region was eight settings away."""
    names = list(ctx.space.names)
    pool = {}
    for c, s in list(ctx.scored()) + list(elsewhere):
        if c is None or c.key == centre.key or s is None or s < score - below:
            continue
        diff = [n for n in names if c.params.get(n) != centre.params.get(n)]
        if len(diff) >= 3 and (c.key not in pool or s > pool[c.key][0]):
            pool[c.key] = (s, c, diff)
    out: List[Tuple[str, int, int, Config]] = []
    seen = set(taken) | {centre.key}
    for s, target, diff in sorted(pool.values(), key=lambda t: (-t[0], t[1].key))[:targets]:
        for size, quota in zip((2, 3, 4), per_size):
            if size >= len(diff):
                break
            combos = list(itertools.combinations(diff, size))
            rng.shuffle(combos)
            got = 0
            for combo in combos:
                if got >= quota:
                    break
                config = ctx.space._assemble({**centre.params, **{n: target.params.get(n) for n in combo}})
                if config is None or config.key in seen:
                    continue
                seen.add(config.key)
                if ctx.score(config) is not None or skip(config):
                    continue
                out.append(("+".join(combo), 0, 0, config))
                got += 1
    return out


# ---- arms carried from round to round ------------------------------------------------------------------------------
# A round's own jumps teach its arms within the round. What every round shares (which settings are rugged, which way a
# setting tends to pay) is carried to the next round's priors: earlier counts are halved each round and capped at
# CARRY_MAX_TRIES per arm, so one round cannot set the next one's choices.
CARRY_DECAY = 0.5
CARRY_MAX_TRIES = 20.0


def _arm_key(arm: tuple) -> str:
    return "|".join(str(x) for x in arm)


def _arm_of(text: str) -> Optional[tuple]:
    parts = text.split("|")
    try:
        return (parts[0], int(parts[1]), int(parts[2])) if len(parts) == 3 else (parts[0], int(parts[1]))
    except (ValueError, IndexError):
        return None


def load_arms(path) -> Arms:
    """Arms whose priors include earlier rounds' counts (missing or unreadable file: plain priors)."""
    import json
    try:
        raw = json.load(open(path))
    except (OSError, ValueError):
        return Arms()
    carried = {}
    for text, counts in (raw.get("arms") or {}).items():
        arm = _arm_of(text)
        if arm is not None and isinstance(counts, list) and len(counts) == 3:
            carried[arm] = [float(x) for x in counts]
    return Arms(carried)


def save_arms(arms: Arms, path) -> None:
    """This round's counts added to the halved carried ones, each arm capped at CARRY_MAX_TRIES tries."""
    import json, os
    out = {}
    for arm in set(arms.stats) | set(arms.carried):
        w, t, g = arms.stats.get(arm, (0.0, 0.0, 0.0))
        cw, ct, cg = arms.carried.get(arm, (0.0, 0.0, 0.0))
        w, t, g = w + CARRY_DECAY * cw, t + CARRY_DECAY * ct, g + CARRY_DECAY * cg
        if t > CARRY_MAX_TRIES:
            f = CARRY_MAX_TRIES / t
            w, t, g = w * f, t * f, g * f
        if t > 0:
            out[_arm_key(arm)] = [round(w, 4), round(t, 4), round(g, 6)]
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = str(path) + ".tmp"
    with open(tmp, "w") as fh:
        json.dump({"arms": out}, fh)
    os.replace(tmp, path)
