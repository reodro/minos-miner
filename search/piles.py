"""Six callset piles (v9.22). Four climb. Two are candidates.

A pile's reference is its highest config within 0.10. A config joins only the nearest pile. A higher score
at distance 0.10 or less becomes that fingerprint. A pile the new fingerprint overlaps is merged into it.
The climber's centre is the highest score inside the ball. A fingerprint farther than
0.10 opens a new pile, up to six, when it is within 1.0 of the best. Four of those piles have a climber.
The other two wait, and the best one takes a seat when a climber reaches the top. A finished pile's slot
then takes the next far fingerprint. A fingerprint whose climber has already reached its peak is not climbed again.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from common.tuning import tune
from common.types import Config
from search.distance import setting_distance
from search.engine_tag import VERSION

RADIUS = 0.10
GAP = 0.10
MAX_PILES = 6
CANDIDATES = 2
POOL = 32
SETTING_DIFFS = 4
FLOOR = 1.0


class Pile:
    def __init__(self, config: Config, score: float):
        self.config, self.score = config, score
        self.members: List[Config] = [config]
        self.scores: Dict[str, float] = {config.key: score}
        self.radius = RADIUS
        self.pool_cap = POOL
        self.absorbed = 0

    def label(self, i: int) -> str:
        return f"pile {i + 1}"


def enabled() -> bool:
    return (VERSION.startswith("v9.19") or VERSION.startswith("v9.20")
            or VERSION.startswith("v9.21") or VERSION.startswith("v9.22")
            or VERSION.startswith("v9.23") or VERSION.startswith("v9.24")
            or VERSION.startswith("v9.25")             or VERSION.startswith("v9.26")
            or VERSION.startswith("v9.27") or VERSION.startswith("v9.28")
            or VERSION.startswith("v9.29") or VERSION.startswith("v9.30")
            or VERSION.startswith("v9.31") or VERSION.startswith("v9.32")
            or VERSION.startswith("v9.33") or VERSION.startswith("v9.34"))


def _dist(ctx, a: Config, b: Config) -> Optional[float]:
    if a is None or b is None:
        return None
    if a.key == b.key:
        return 0.0
    dist = getattr(ctx, "dist", None)
    if dist is not None:
        d, _kind = dist.between(a, b)
        return float(d)
    space = getattr(ctx, "space", None)
    if space is None:
        return None
    return setting_distance(space, a, b)


def _diffs(ctx, a: Config, b: Config) -> int:
    names = list(getattr(getattr(ctx, "space", None), "names", ()) or [])
    return sum(a.params.get(n) != b.params.get(n) for n in names)


def _gap_ok(ctx, config: Config, piles: Sequence[Pile], skip: int) -> bool:
    for i, pile in enumerate(piles):
        if i == skip:
            continue
        d = _dist(ctx, config, pile.config)
        if d is not None and d < GAP:
            return False
    return True


def _history(ctx) -> List[dict]:
    hist = getattr(ctx, "_fingerprints", None)
    if hist is None:
        hist = []
        ctx._fingerprints = hist
    return hist


def _say(ctx, msg: str, *args) -> None:
    say = getattr(ctx, "say", None)
    if callable(say):
        say(msg, *args)


def _trail(ctx, pile_no: int) -> str:
    parts = []
    for entry in _history(ctx):
        if entry["pile"] != pile_no:
            continue
        mark = " closed" if entry["closed"] else ""
        parts.append(f"{entry['config'].key[:8]} {entry['score']:.4f}{mark}")
    return ", ".join(parts)


def remember(ctx, pile_no: int, config: Config, score: float) -> None:
    """Append this pile's fingerprint. The same key only updates its score."""
    hist = _history(ctx)
    own = [entry for entry in hist if entry["pile"] == pile_no]
    if own and own[-1]["config"].key == config.key:
        if score > own[-1]["score"]:
            own[-1]["score"] = float(score)
        return
    hist.append({"pile": pile_no, "config": config, "score": float(score), "closed": False})
    _say(ctx, "fingerprint pile %d: %s", pile_no, _trail(ctx, pile_no))


def seal_peak(climber) -> None:
    """The climber found this fingerprint's peak. The history keeps it, so a similar pile is not opened."""
    if not enabled():
        return
    centre = getattr(climber, "centre", None)
    if centre is None:
        return
    origin = getattr(climber, "origin", "") or ""
    try:
        pile_no = int(origin.split()[1])
    except (IndexError, ValueError):
        return
    if pile_no < 1:
        return
    ctx = climber.ctx
    remember(ctx, pile_no, centre, float(getattr(climber, "score", 0.0)))
    closed = False
    for entry in _history(ctx):
        if entry["pile"] == pile_no and entry["config"].key == centre.key and not entry["closed"]:
            entry["closed"] = True
            closed = True
    if closed:
        _say(ctx, "fingerprint pile %d closed %s %.4f; history: %s",
             pile_no, centre.key[:8], float(getattr(climber, "score", 0.0)), _trail(ctx, pile_no))


def similar_fingerprint(ctx, config: Config, own_pile: Optional[int] = None):
    """The closest remembered fingerprint within 0.10, aside from this pile's own live trail.

    A closed fingerprint of this pile still counts: that pile has been explored.
    """
    best = None
    for entry in _history(ctx):
        if own_pile is not None and entry["pile"] == own_pile and not entry["closed"]:
            continue
        if entry["config"].key == config.key:
            dist = 0.0
        else:
            dist = _dist(ctx, config, entry["config"])
            if dist is None:
                continue
        if dist < GAP and (best is None or dist < best[0]):
            best = (dist, entry)
    return best


def _refuse(ctx, config: Config, score: float, hit) -> None:
    seen = getattr(ctx, "_fingerprint_refused", None)
    if seen is None:
        seen = set()
        ctx._fingerprint_refused = seen
    if config.key in seen:
        return
    seen.add(config.key)
    dist, entry = hit
    _say(ctx, "fingerprint %s %.4f is similar to pile %d %s (%.3f); it is not explored",
         config.key[:8], score, entry["pile"], entry["config"].key[:8], dist)


def history_rows(ctx) -> List[dict]:
    return [{"pile": entry["pile"], "key": entry["config"].key[:8],
             "score": round(entry["score"], 4), "closed": entry["closed"]}
            for entry in _history(ctx)]


def _owned(ctx, piles: Sequence[Pile], i: int, climbing) -> bool:
    pile = piles[i]
    for cl in climbing:
        centre = getattr(cl, "centre", None)
        if (getattr(cl, "origin", "") or "") == pile.label(i):
            return True
        if centre is None:
            continue
        d = _dist(ctx, centre, pile.config)
        if d is not None and d <= RADIUS:
            return True
    return False


def pile_rows(ctx) -> List[dict]:
    """Each pile and its three highest scores, for the live dashboard."""
    piles = list(getattr(ctx, "_piles", None) or [])
    climbing = [cl for cl in getattr(ctx, "climbers", ()) if getattr(cl, "state", "") == "climbing"]
    rows = []
    for i, pile in enumerate(piles):
        ranked = sorted(pile.scores.items(), key=lambda item: (-item[1], item[0]))[:3]
        rows.append({
            "pile": i + 1,
            "key": pile.config.key,
            "score": round(float(pile.score), 4),
            "candidate": bool(climbing) and not _owned(ctx, piles, i, climbing),
            "top": [{"key": key, "score": round(float(score), 4)} for key, score in ranked],
        })
    return rows


def note(ctx, piles: List[Pile], config: Config, score: float, best: float) -> None:
    """Fold one scored config into the piles. A higher score within 0.10 replaces that reference."""
    if not piles:
        piles.append(Pile(config, score))
        remember(ctx, 1, config, score)
        return
    dists = [_dist(ctx, config, p.config) for p in piles]
    if any(d is None for d in dists):
        return
    nearest = min(range(len(dists)), key=lambda i: dists[i])
    if dists[nearest] <= float(getattr(piles[nearest], "radius", RADIUS)):
        pile = piles[nearest]
        if score > pile.score and not _grid_owns(ctx, pile):
            claim_peak(ctx, piles, nearest, config, score)
            return
        if config.key not in pile.scores:
            pile.members.append(config)
        pile.scores[config.key] = max(score, pile.scores.get(config.key, score))
        return
    if all(d >= GAP for d in dists) and best - score <= FLOOR:
        hit = similar_fingerprint(ctx, config)
        if hit is not None:
            _refuse(ctx, config, score, hit)
            return
        if len(piles) < MAX_PILES:
            piles.append(Pile(config, score))
            remember(ctx, len(piles), config, score)
            return
        slot = _open_slot(ctx, piles, score)
        if slot is None:
            return
        piles[slot] = Pile(config, score)
        remember(ctx, slot + 1, config, score)


def refresh(ctx, rows: Optional[Sequence[Tuple[Config, float]]] = None) -> List[Pile]:
    """Piles from `rows`, or from every score the context has seen. Later calls only fold in new keys."""
    if rows is None:
        raw = list(ctx.scored()) if hasattr(ctx, "scored") else []
    else:
        raw = list(rows)
    raw.sort(key=lambda t: -t[1])
    if rows is not None:
        piles: List[Pile] = []
        best = raw[0][1] if raw else 0.0
        for config, score in raw:
            note(ctx, piles, config, score, best)
        return piles
    piles = list(getattr(ctx, "_piles", None) or [])
    seen = getattr(ctx, "_pile_seen", None)
    if seen is None:
        seen = set()
        ctx._pile_seen = seen
    fresh = [(c, s) for c, s in raw if c.key not in seen]
    best = raw[0][1] if raw else 0.0
    for config, score in fresh:
        note(ctx, piles, config, score, best)
        seen.add(config.key)
    ctx._piles = piles
    return piles


def mating_pool(ctx, pile: Pile) -> List[Config]:
    """Up to 32 parents in this pile that differ in at least four searched settings."""
    members = sorted(pile.members, key=lambda c: -pile.scores.get(c.key, pile.score))
    chosen: List[Config] = []
    for config in members:
        if all(_diffs(ctx, config, other) >= SETTING_DIFFS for other in chosen):
            chosen.append(config)
        if len(chosen) >= int(getattr(pile, "pool_cap", POOL)):
            break
    return chosen


def absorb_parents(ctx, higher: Pile, lower: Pile) -> None:
    """The higher pile takes the lower climber's parents. Its radius and parent cap grow to fit them."""
    if higher is lower:
        return
    for config in mating_pool(ctx, lower):
        score = float(lower.scores.get(config.key, lower.score))
        if config.key not in higher.scores:
            higher.members.append(config)
        higher.scores[config.key] = max(score, float(higher.scores.get(config.key, score)))
        dist = _dist(ctx, higher.config, config)
        if dist is not None:
            higher.radius = max(float(higher.radius), float(dist))
    higher.absorbed = min(3, int(higher.absorbed) + 1)
    higher.pool_cap = POOL * (1 + higher.absorbed)


def own_pile(ctx, climber) -> Optional[Pile]:
    piles = list(getattr(ctx, "_piles", None) or [])
    if not piles:
        return None
    origin = getattr(climber, "origin", "") or ""
    if origin.startswith("pile "):
        try:
            i = int(origin.split()[1]) - 1
        except (IndexError, ValueError):
            i = -1
        if 0 <= i < len(piles):
            return piles[i]
    centre = getattr(climber, "centre", None)
    if centre is None:
        return None
    dists = [(i, _dist(ctx, centre, p.config)) for i, p in enumerate(piles)]
    dists = [(i, d) for i, d in dists if d is not None]
    if not dists:
        return None
    i, d = min(dists, key=lambda t: t[1])
    return piles[i] if d <= RADIUS else None


def _grid_owns(ctx, pile: Pile) -> bool:
    """A climber is walking this pile's settings. A higher score joins the pile and leaves the fingerprint."""
    from search.climb import ga_hot
    if ga_hot(ctx):
        return False
    for cl in getattr(ctx, "climbers", ()) or ():
        if getattr(cl, "state", "") != "climbing" or not getattr(cl, "batch_live", False):
            continue
        centre = getattr(cl, "centre", None)
        if centre is None:
            continue
        d = _dist(ctx, centre, pile.config)
        if d is not None and d <= RADIUS:
            return True
    return False


def _absorb_pile(ctx, higher: Pile, lower: Pile) -> None:
    """Every member of the lower pile joins the higher one. The higher radius grows to cover them."""
    if higher is lower:
        return
    for config in list(lower.members):
        score = float(lower.scores.get(config.key, lower.score))
        if config.key not in higher.scores:
            higher.members.append(config)
        higher.scores[config.key] = max(score, float(higher.scores.get(config.key, score)))
        dist = _dist(ctx, higher.config, config)
        if dist is not None:
            higher.radius = max(float(higher.radius), float(dist))
    higher.absorbed = min(3, int(higher.absorbed) + 1)
    higher.pool_cap = POOL * (1 + higher.absorbed)


def _retarget(ctx, piles: List[Pile], winner: Pile, dropped: Sequence[int], winner_was: int) -> None:
    """Climbers follow the pile list after a merge. A climber on an absorbed pile sits on the winner."""
    winner_now = piles.index(winner)

    def new_index(old: int) -> int:
        if old in dropped:
            return winner_now
        return old - sum(1 for i in dropped if i < old)

    for cl in getattr(ctx, "climbers", ()) or ():
        origin = getattr(cl, "origin", "") or ""
        if not origin.startswith("pile "):
            continue
        try:
            old = int(origin.split()[1]) - 1
        except (IndexError, ValueError):
            continue
        if old == winner_was or old in dropped or any(i < old for i in dropped):
            cl.origin = winner.label(new_index(old)) if old in dropped else piles[new_index(old)].label(new_index(old))


def merge_overlaps(ctx, piles: List[Pile], index: int) -> int:
    """Piles within 0.10 of this fingerprint join it. The highest score keeps the fingerprint."""
    if index < 0 or index >= len(piles):
        return index
    winner = piles[index]
    dropped = []
    for i, other in enumerate(piles):
        if i == index:
            continue
        dist = _dist(ctx, winner.config, other.config)
        if dist is not None and dist < GAP:
            _absorb_pile(ctx, winner, other)
            dropped.append(i)
            _say(ctx, "pile %d %.4f absorbs pile %d %.4f (%.3f)", index + 1, float(winner.score), i + 1, float(other.score), dist)
    if not dropped:
        return index
    dropped_nos = {i + 1 for i in dropped}
    for entry in _history(ctx):
        if entry["pile"] in dropped_nos:
            entry["closed"] = True
    for i in sorted(dropped, reverse=True):
        del piles[i]
    new_index = piles.index(winner)
    old_no = index + 1
    new_no = new_index + 1
    if old_no != new_no:
        for entry in _history(ctx):
            if entry["pile"] == old_no and not entry["closed"]:
                entry["pile"] = new_no
    _retarget(ctx, piles, winner, dropped, index)
    return new_index


def claim_peak(ctx, piles: List[Pile], index: int, config: Config, score: float) -> int:
    """The highest score becomes the fingerprint. Overlapping piles merge into it."""
    pile = piles[index]
    if config.key not in pile.scores:
        pile.members.append(config)
    pile.scores[config.key] = max(float(score), float(pile.scores.get(config.key, float(score))))
    if float(score) >= float(pile.score):
        pile.config, pile.score = config, float(score)
    index = merge_overlaps(ctx, piles, index)
    pile = piles[index]
    remember(ctx, index + 1, pile.config, float(pile.score))
    return index


def best_member(pile: Pile) -> Tuple[Config, float]:
    """The highest scored config in this pile."""
    best_c, best_s = pile.config, float(pile.score)
    for config in pile.members:
        score = float(pile.scores.get(config.key, pile.score))
        if score > best_s:
            best_c, best_s = config, score
    return best_c, best_s


def follow(ctx, climber, config: Config, score: float) -> bool:
    """The climber may stand on `config` when it is within 0.10 of this pile.

    A higher score becomes the fingerprint. Any pile that fingerprint overlaps is merged into this one.
    """
    if not enabled():
        return True
    piles = getattr(ctx, "_piles", None)
    if not piles:
        return True
    own = own_pile(ctx, climber)
    if own is None:
        return True
    idx = piles.index(own)
    d_own = _dist(ctx, config, own.config)
    inside = d_own is None or d_own <= RADIUS
    if not inside and not _gap_ok(ctx, config, piles, idx):
        return False
    if score >= own.score:
        claim_peak(ctx, piles, idx, config, score)
        return True
    if config.key not in own.scores:
        own.members.append(config)
    own.scores[config.key] = max(score, own.scores.get(config.key, score))
    return True


def centre_allowed(climber, config: Config) -> bool:
    """A climber may take `config` as its centre when it is within 0.10 of this pile's fingerprint."""
    if not enabled():
        return True
    ctx = climber.ctx
    piles = list(getattr(ctx, "_piles", None) or [])
    if not piles:
        return True
    own = own_pile(ctx, climber)
    if own is None:
        return True
    d_own = _dist(ctx, config, own.config)
    if d_own is None:
        return True
    return d_own <= RADIUS


def home(ctx, config: Config) -> Optional[int]:
    """Index of the pile this fingerprint joins, or None when it is outside every ball."""
    piles = list(getattr(ctx, "_piles", None) or [])
    dists = [_dist(ctx, config, p.config) for p in piles]
    if not dists or any(d is None for d in dists):
        return None
    i = min(range(len(dists)), key=lambda k: dists[k])
    return i if dists[i] <= RADIUS else None


def already_peaked(ctx, pile: Pile) -> bool:
    """A climber already found the peak of this fingerprint.

    The same config is not climbed again. A different config inside the ball is climbed only when it beats
    that peak by climb.accept_min.
    """
    accept = float(tune("climb.accept_min"))
    best_near = None
    for cl in getattr(ctx, "climbers", ()):
        if getattr(cl, "state", "") != "peak":
            continue
        centre = getattr(cl, "centre", None)
        if centre is None:
            continue
        if centre.key == pile.config.key:
            return True
        d = _dist(ctx, centre, pile.config)
        if d is None or d > RADIUS:
            continue
        score = float(getattr(cl, "score", pile.score))
        if best_near is None or score > best_near:
            best_near = score
    return best_near is not None and pile.score - best_near < accept


def _note_peaked(ctx, i: int, pile: Pile) -> None:
    seen = getattr(ctx, "_peaked_note", None)
    if seen is None:
        seen = set()
        ctx._peaked_note = seen
    if pile.config.key in seen:
        return
    seen.add(pile.config.key)
    say = getattr(ctx, "say", None)
    if callable(say):
        say("pile %d peaked at %s %.4f; that fingerprint is not climbed again",
            i + 1, pile.config.key[:8], pile.score)


def _finished_slot(ctx, piles: Sequence[Pile]) -> Optional[int]:
    """The lowest finished pile. A new fingerprint takes its place so two candidates stay ready."""
    climbing = [cl for cl in getattr(ctx, "climbers", ()) if getattr(cl, "state", "") == "climbing"]
    finished = []
    for i, pile in enumerate(piles):
        if _owned(ctx, piles, i, climbing) or not already_peaked(ctx, pile):
            continue
        finished.append((pile.score, i))
    if not finished:
        return None
    return min(finished)[1]


def _lowest_climbing_score(ctx, piles: Sequence[Pile]) -> Optional[float]:
    climbing = [cl for cl in getattr(ctx, "climbers", ()) if getattr(cl, "state", "") == "climbing"]
    scores = []
    for cl in climbing:
        own = own_pile(ctx, cl)
        if own is None or own not in piles:
            continue
        scores.append(best_member(own)[1])
    if not scores:
        return None
    return min(scores)


def _open_slot(ctx, piles: Sequence[Pile], score: float) -> Optional[int]:
    """Where a new fingerprint sits once six piles exist.

    A peaked pile still gives its slot up. Otherwise a config above the lowest climbing pile
    replaces the lowest candidate, so the climber can move now. A peak is not required.
    """
    slot = _finished_slot(ctx, piles)
    if slot is not None:
        return slot
    low = _lowest_climbing_score(ctx, piles)
    if low is None or score <= low:
        return None
    climbing = [cl for cl in getattr(ctx, "climbers", ()) if getattr(cl, "state", "") == "climbing"]
    weak = [(float(pile.score), i) for i, pile in enumerate(piles) if not _owned(ctx, piles, i, climbing)]
    if not weak:
        return None
    cand_score, i = min(weak)
    if score <= cand_score:
        return None
    return i


def unseated(ctx, piles: Sequence[Pile]) -> List[Tuple[Config, float, str]]:
    """Piles with no climber, best score first.

    A pile whose climber already reached the peak of its fingerprint is left out. The two candidates
    are the piles that are still waiting. The best one takes a seat when a climber reaches the top.
    """
    out = []
    for cl in getattr(ctx, "climbers", ()):
        if getattr(cl, "state", "") == "peak":
            seal_peak(cl)
    climbers = [cl for cl in getattr(ctx, "climbers", ()) if getattr(cl, "state", "") == "climbing"]
    for i, pile in enumerate(piles):
        if _owned(ctx, piles, i, climbers):
            continue
        if already_peaked(ctx, pile):
            _note_peaked(ctx, i, pile)
            continue
        hit = similar_fingerprint(ctx, pile.config, own_pile=i + 1)
        if hit is not None:
            _refuse(ctx, pile.config, pile.score, hit)
            continue
        out.append((pile.config, pile.score, pile.label(i)))
    out.sort(key=lambda item: -item[1])
    return out
