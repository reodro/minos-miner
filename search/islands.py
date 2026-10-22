"""Islands: independent searches sharing one fleet.

044200Z (2026-10-16): each run settles on one region in its first 15-20 minutes and all its climbers end there
(810 workers: three climbers at 82.632-82.633; 1056: four at 82.638). The two 324-worker runs reached a region that
goes to 82.6953. The 810 run made no gain in its last 28 minutes (17,531 runs). One search per islands.workers workers
turns a bigger fleet into more regions searched, instead of the same region reached sooner.

    count     cap // islands.workers, at least 1: 810 is two islands of 405; 324 is one search, as before
    view      each island is a whole engine run (Step A, B, C, certify) on its own IslandClient. It proposes through
              the shared client and sees only the configs it proposed: latest, configs, callsets, on_result, its
              best. A config another island has already scored comes back at once; one still running is run once
              (the executor joins the second job to the first)
    workers   an island's idle_slots() leaves room for every other island that is still under its share and has
              proposed in the last ACTIVE_S seconds; while all of them are active it also stays inside its own share
    apart     the portfolio seeds are dealt by rank (search/candidates.py seeds), the defaults and the official config
              go to island 0, each island has its own random stream, and a Step A seat that differs in fewer than
              islands.min_settings settings from another island's seat is passed over while the island has another
    pick      once an island names its pick, client.pick_key is the highest score among every island's pick and best.
              The notes are the best island's, with one summary row per island and every island's climbers (their ids
              start at 100 x the island's index)
"""
from __future__ import annotations

import asyncio
import logging
import math
from typing import Callable, Dict, List, Optional, Sequence

from common.tuning import tune
from common.types import FULL, Config

log = logging.getLogger(__name__)
ACTIVE_S = 30.0          # an island that has proposed nothing for this long leaves its unused share to the others
TICK_S = 2.0


def count(cap: int) -> int:
    """How many islands `cap` workers make."""
    per = max(1, int(tune("islands.workers")))
    return max(1, int(cap or 0) // per)


def quotas(cap: int, k: int) -> List[int]:
    return [cap // k + (1 if i < cap % k else 0) for i in range(k)]


def weighted_quotas(cap: int, weights: Sequence[int]) -> List[int]:
    """cap split in proportion to weights, the remainder to the largest fractions."""
    total = sum(weights) or len(weights)
    raw = [cap * w / total for w in weights]
    out = [int(math.floor(x)) for x in raw]
    for i in sorted(range(len(raw)), key=lambda i: out[i] - raw[i])[:cap - sum(out)]:
        out[i] += 1
    return out


def climbing_in(view) -> int:
    """Climbers still climbing on one island, from its published notes."""
    eng = view.notes.get("engine") or {}
    return sum(1 for b in eng.get("basins") or [] if isinstance(b, dict) and b.get("status") == "climbing")


def rebalance(views: Sequence["IslandClient"], cap: int) -> None:
    """islands.by_climber (user, 2026-10-19): every climbing climber gets the same share of the fleet, whichever island
    it is on. An island's quota is its climbing climbers' share (one climber's worth when it has none, for its GA,
    RF-UCB and the next seat); the island's engine reads it back as its capacity (Context.publish)."""
    if not int(tune("islands.by_climber")) or len(views) < 2:
        return
    weights = [max(1, climbing_in(v)) for v in views]
    for v, q in zip(views, weighted_quotas(cap, weights)):
        v.quota = q


class Island:
    """What an island's engine run knows about the islands: its index, how many there are, and their seats."""

    def __init__(self, index: int, k: int, seats: Dict[int, List[Config]], arms=None, views=None):
        self.index, self.k, self.seats = index, k, seats
        self.arms = arms                                              # the centre-jump arms every island shares
        self.views = views or []

    def elsewhere(self) -> List[tuple]:
        """(config, score) of every other island's best: path-relinking targets at a peak (search/leaps.py relink)."""
        out = []
        for v in self.views:
            if v.index == self.index:
                continue
            key = _best_key(v)
            r = v.latest.get((key, FULL)) if key else None
            c = v.configs.get(key) if key else None
            if r is not None and c is not None and r.score is not None:
                out.append((c, float(r.score)))
        return out

    def apart(self, ctx, starts: Sequence[Config]) -> List[Config]:
        """Step A's starts without the ones that differ from another island's seat in fewer than islands.min_settings
        settings, while at least one start remains. The starts kept become this island's seats."""
        starts = list(starts)
        limit = int(tune("islands.min_settings"))
        others = [c for i, cs in self.seats.items() if i != self.index for c in cs]
        if others and limit > 0 and starts:
            names = list(ctx.space.names)

            def apart_from_all(c: Config) -> bool:
                return all(sum(c.params.get(n) != o.params.get(n) for n in names) >= limit for o in others)

            # The island's best start always stays (044200Z islands, 2026-10-16: island 1 passed over its 82.1127 and
            # opened on one weaker seat); only the later ones must be apart.
            kept = starts[:1] + [c for c in starts[1:] if apart_from_all(c)]
            if len(kept) < len(starts):
                ctx.say("islands: %d of %d later step A seat(s) differ from another island's seat in fewer than %d "
                        "settings; passed over (the best start stays)", len(starts) - len(kept), len(starts), limit)
            starts = kept
        self.seats[self.index] = list(starts)
        return starts


class IslandClient:
    """One island's view of the shared EvalClient: its own proposals, results, workers, notes and pick."""

    def __init__(self, base, index: int, quota: int, peers: List["IslandClient"]):
        self.base, self.index, self.quota, self.peers = base, index, int(quota), peers
        self.round, self.ranges, self.store, self.executor = base.round, base.ranges, base.store, base.executor
        self.clock = base.clock
        self.budget, self.deadline_ts = base.budget, base.deadline_ts
        self.own: set = set()                                         # keys this island proposed
        self.latest: Dict[tuple, object] = {}
        self.callsets: Dict[tuple, object] = {}
        self.configs: Dict[str, Config] = {}
        self.on_result: List[Callable] = []
        self.on_rescore: List[Callable] = []
        self.notes: Dict[str, object] = {}
        self._pick: Optional[str] = None
        self._open: List[tuple] = []                                  # (futures of one propose call, its cost in runs)
        self._futs: Dict[asyncio.Future, str] = {}                    # this island's futures still out -> config key
        self.last_propose = self.clock()
        base.on_result.append(self._heard)
        base.on_rescore.append(self._rescored)

    # ---- what the engine reads --------------------------------------------------------------------------------------
    @property
    def answer_set(self):
        return self.base.answer_set

    @property
    def spent(self) -> float:
        return self.base.spent

    @property
    def exhausted(self) -> bool:
        return self.base.exhausted

    def time_left(self) -> Optional[float]:
        return self.base.time_left()

    def capacity(self) -> Optional[int]:
        return self.quota

    def variant_cost(self, main: Config, variant: Config, fidelity: str = FULL) -> float:
        return self.base.variant_cost(main, variant, fidelity)

    @property
    def in_flight(self) -> int:
        self._prune()
        return len(self._futs)

    def running(self) -> set:
        self._prune()
        return set(self._futs.values())

    def best(self, fidelity: str = FULL):
        top = None
        for (key, fid), r in self.latest.items():
            if fid == fidelity and r.score is not None and r.gate_passed is not False and (top is None or r.score > top.score):
                top = r
        return top

    def picked(self):
        return self.latest.get((self._pick, FULL)) if self._pick else None

    def _note_best(self, result) -> None:
        """The shared client keeps the bests; an island's best is read from its own results."""

    @property
    def pick_key(self) -> Optional[str]:
        return self._pick

    @pick_key.setter
    def pick_key(self, key: Optional[str]) -> None:
        self._pick = key
        name_pick(self.base, self.peers)

    # ---- proposing ----------------------------------------------------------------------------------------------------
    def propose(self, config: Config, *, fidelity: str = FULL, priority: int = 0, tag: str = "",
                batch: str = "big") -> asyncio.Future:
        return self.propose_many([config], fidelity=fidelity, priority=priority, tag=tag, batch=batch)[0]

    def propose_many(self, configs, *, fidelity: str = FULL, priority: int = 0, tag: str = "", batch: str = "big",
                     batches: Optional[Sequence[str]] = None) -> List[asyncio.Future]:
        configs = list(configs)
        for c in configs:
            self.own.add(c.key)
        before = self.base.spent
        futs = self.base.propose_many(configs, fidelity=fidelity, priority=priority, tag=tag, batch=batch,
                                      batches=batches)
        cost = max(0.0, float(self.base.spent) - float(before))
        self.last_propose = self.clock()
        out = [f for f in futs if not f.done()]
        if out:
            self._open.append((out, cost))
        for c, f in zip(configs, futs):
            if not f.done():
                self._futs[f] = c.key
                continue
            outcome = f.result()                                       # already scored (by another island) or invalid
            if getattr(outcome, "result", None) is not None:
                self._copy(c.key, fidelity)
                for cb in list(self.on_result):
                    cb(outcome)
        return futs

    async def propose_all(self, configs, *, quorum: float = 1.0, fidelity: str = FULL, priority: int = 0,
                          tag: str = "") -> list:
        futs = self.propose_many(configs, fidelity=fidelity, priority=priority, tag=tag)
        need = math.ceil(quorum * len(futs))
        pending = set(futs)
        while len(futs) - len(pending) < need:
            _, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
        return [f.result() if f.done() else None for f in futs]

    # ---- workers ------------------------------------------------------------------------------------------------------
    def outstanding(self) -> float:
        """Runs this island has out (a rider counts as its share of a run)."""
        self._open = [(fs, cost) for fs, cost in self._open if any(not f.done() for f in fs)]
        return sum(cost for _, cost in self._open)

    def idle_slots(self) -> Optional[int]:
        free = self.base.idle_slots()
        if free is None:
            return None
        now = self.clock()
        others, everyone_active = 0.0, True
        for p in self.peers:
            if p is self:
                continue
            if now - p.last_propose <= ACTIVE_S:
                others += max(0.0, p.quota - p.outstanding())
            else:
                everyone_active = False                               # its unused share is free for the others
        room = free - others
        if everyone_active:
            room = min(room, self.quota - self.outstanding())
        return max(0, int(math.floor(room)))

    # ---- results ------------------------------------------------------------------------------------------------------
    def _prune(self) -> None:
        for f in [f for f in self._futs if f.done()]:
            del self._futs[f]

    def _copy(self, key: str, fidelity: str = FULL) -> None:
        c = self.base.configs.get(key)
        if c is not None:
            self.configs[key] = c
        for fid in {FULL, fidelity}:
            r = self.base.latest.get((key, fid))
            if r is not None:
                self.latest[(key, fid)] = r
            cs = self.base.callsets.get((key, fid))
            if cs is not None:
                self.callsets[(key, fid)] = cs

    def _heard(self, outcome) -> None:
        c = getattr(outcome, "config", None)
        if c is None or c.key not in self.own:
            return
        result = getattr(outcome, "result", None)
        self._copy(c.key, getattr(result, "fidelity", None) or FULL)
        for cb in list(self.on_result):
            cb(outcome)

    def _rescored(self) -> None:
        for (key, fid) in list(self.latest):
            r = self.base.latest.get((key, fid))
            if r is not None:
                self.latest[(key, fid)] = r
        for cb in list(self.on_rescore):
            cb()

    def close(self) -> None:
        for lst, cb in ((self.base.on_result, self._heard), (self.base.on_rescore, self._rescored)):
            if cb in lst:
                lst.remove(cb)


def _score(base, key: Optional[str]) -> Optional[float]:
    if not key:
        return None
    r = base.latest.get((key, FULL))
    if r is None or r.score is None or r.gate_passed is False:
        return None
    return float(r.score)


def _best_key(view: IslandClient) -> Optional[str]:
    top = None
    for (key, fid), r in view.latest.items():
        if fid == FULL and r.score is not None and r.gate_passed is not False and (top is None or r.score > top[1]):
            top = (key, r.score)
    return None if top is None else top[0]


def name_pick(base, views: Sequence[IslandClient]) -> None:
    """client.pick_key, once any island has named its pick: the highest score among every island's pick and best.

    The pick is the highest local score. An island whose clock jumped past its pick time (a simulated job) has not
    named one, and its best still counts."""
    if not any(v.pick_key for v in views):
        return
    best = None
    for v in views:
        for key in (v.pick_key, _best_key(v)):
            s = _score(base, key)
            if s is not None and (best is None or s > best[0]):
                best = (s, key)
    if best is not None and base.pick_key != best[1]:
        base.pick_key = best[1]


def island_distances(rows: Sequence[dict]) -> List[dict]:
    """Between every two islands' best configs: the unit-cube distance and how many searched settings differ."""
    out = []
    for i, a in enumerate(rows):
        for b in rows[i + 1:]:
            va, vb = a.get("best_vec"), b.get("best_vec")
            if not va or not vb or len(va) != len(vb):
                continue
            out.append({"a": a["island"], "b": b["island"],
                        "distance": round(sum((x - y) ** 2 for x, y in zip(va, vb)) ** 0.5, 4),
                        "settings": sum(1 for x, y in zip(va, vb) if abs(x - y) > 1e-9)})
    return out


def race_cull(base, views: Sequence[IslandClient], made: Dict[int, "Island"]) -> None:
    """islands.race_k (generalized Step A, 2026-10-22): when the race ends, the `keep` islands with the highest best
    score climb on; the others stop (their climbers end, their workers go to the winners through rebalance)."""
    keep = max(1, int(getattr(next(iter(made.values())), "keep", 2)))
    ranked = sorted(views, key=lambda v: -(_score(base, _best_key(v)) or float("-inf")))
    winners = {v.index for v in ranked[:keep]}
    for v in views:
        if v.index not in winners and v.index in made:
            made[v.index].stopped = True
    log.info("islands: race over: %s climb on; %s stop", sorted(winners), sorted(set(made) - winners))


def merge_notes(base, views: Sequence[IslandClient], cap: int) -> None:
    """client.notes["engine"]: the best island's notes, every island's summary and every island's climbers."""
    rows, lead, lead_best = [], None, None
    for v in views:
        eng = v.notes.get("engine")
        if not isinstance(eng, dict):
            continue
        b = eng.get("best")
        if lead is None or (b is not None and (lead_best is None or b > lead_best)):
            lead, lead_best = eng, b
        rows.append(eng)
    if lead is None:
        return
    merged = dict(lead)
    merged["capacity"] = cap
    merged["islands"] = [{"island": v.index, "capacity": v.quota,
                          "best": (v.notes.get("engine") or {}).get("best"),
                          "phase": (v.notes.get("engine") or {}).get("phase"),
                          "pick": (v.notes.get("engine") or {}).get("pick"),
                          "climbers": len((v.notes.get("engine") or {}).get("basins") or []),
                          "best_vec": (v.notes.get("engine") or {}).get("best_vec"),
                          "sources": (v.notes.get("engine") or {}).get("sources"),
                          "cloud": (v.notes.get("engine") or {}).get("cloud") or []} for v in views]
    merged["island_distances"] = island_distances(merged["islands"])
    basins, events = [], []
    for v in views:
        eng = v.notes.get("engine") or {}
        basins += [dict(row, island=v.index) for row in eng.get("basins") or [] if isinstance(row, dict)]
        events += list(eng.get("events") or [])
    merged["basins"] = basins
    merged["events"] = sorted(events, key=lambda e: e[0] if isinstance(e, (list, tuple)) and e else 0)[-200:]
    base.notes["engine"] = merged


async def _guarded(coro, index: int, board: Optional[dict] = None) -> None:
    try:
        await coro
    except asyncio.CancelledError:
        raise
    except Exception:  # noqa: BLE001 - one island's failure must not end the others
        log.exception("engine: island %d failed", index)
    finally:
        if board is not None:                           # far.step_c_top: the others stop waiting for this island
            board.setdefault("done", set()).add(index)


async def run(base, cap: int, k: int, start: Callable[[IslandClient, Island], object],
              attrs: Optional[Sequence[dict]] = None) -> List[IslandClient]:
    """Run k islands to the end of the round. `start(view, island)` is one whole engine run. `attrs[i]`: extra
    attributes of island i (offline test T1: its move_k)."""
    views: List[IslandClient] = []
    for i, q in enumerate(quotas(cap, k)):
        views.append(IslandClient(base, i, q, views))
    seats: Dict[int, List[Config]] = {}
    from search.engine import carry_arms, round_arms
    arms = round_arms(base)
    board: dict = {}                                  # shared by the islands (far.step_c_top: the round's top climbers)

    made: Dict[int, Island] = {}

    def island_of(v):
        isl = Island(v.index, k, seats, arms, views)
        made[v.index] = isl
        isl.board = board
        for name, value in ((attrs[v.index] if attrs and v.index < len(attrs) else None) or {}).items():
            setattr(isl, name, value)
        return isl
    tasks = [asyncio.ensure_future(_guarded(start(v, island_of(v)), v.index, board)) for v in views]
    virtual = hasattr(base.executor, "now")
    try:
        culled = False
        while True:
            race_end = next((getattr(i, "race_end", None) for i in made.values()), None)
            if not culled and race_end is not None and base.clock() >= race_end:
                culled = True
                race_cull(base, views, made)
            rebalance(views, cap)
            merge_notes(base, views, cap)
            name_pick(base, views)
            if all(t.done() for t in tasks):
                break
            if virtual:
                await asyncio.sleep(0)
                await asyncio.wait(tasks, timeout=0.05)
            else:
                await asyncio.wait(tasks, timeout=TICK_S)
    finally:
        for t in tasks:
            if not t.done():
                t.cancel()
        for v in views:
            v.close()
    merge_notes(base, views, cap)
    name_pick(base, views)
    base.on_result.append(lambda outcome: _late_pick(base, outcome))
    carry_arms(base, arms)
    return views


def _late_pick(base, outcome) -> None:
    """A run that comes back after every island has stopped still becomes the pick when it scores higher: the pick is
    the highest local score."""
    if not base.pick_key:
        return
    r = getattr(outcome, "result", None)
    c = getattr(outcome, "config", None)
    if r is None or c is None or r.score is None or r.gate_passed is False:
        return
    cur = _score(base, base.pick_key)
    if cur is None or r.score > cur:
        base.pick_key = c.key
