"""Engine v3.2: isolated iterated local search with racing, under a global search that never stops. Owner: search.

The objective is discrete and expensive, the fleet is wide and a round allows ~15-20 dependent steps, so:

  opening     two global waves, 0-OPENING_S: wave 1 = seeds + portfolio + a Latin hypercube over every axis (each
              axis stratified across its whole range, every pair of axes jointly covered), 8 genotyping riders each;
              wave 2 (from WAVE2_AT, or once wave 1 is mostly back) = points chosen by the landscape model, half where
              it predicts well and half where it is most uncertain. A clear wave-1 leader (EARLY_LEAD above every far
              rival) founds its domain at once; the rest are founded when the opening ends: the DIVISION axes (the 3
              live axes with the widest spread among the good configs) and up to N_DOMAINS domains at the best configs
              pairwise far on them, above the hill margin
  scout       SCOUT_SHARE of the workers for the whole round: candidates (fresh hypercube points, variations of the
              best configs) ranked by the model's upper confidence bound - predicted score plus uncertainty - so it is
              pulled to the regions nobody has sampled; never inside a standing domain's hill. A scout find above the
              hill margin in unowned territory is the candidate for the next free domain slot
  domains     each an isolated iterated local search (domain.py): coarse and fine sweeps, best-of-sweep acceptance
              with merges, sideways moves, shakes only when exhausted, batches with a bracket, siblings on structural
              moves. A config belongs to the domain whose centre is nearest on the division axes; a domain never
              proposes into another's territory; nothing crosses between live domains
  merges      one hill by distance (fewer than MIN_FAR live axes apart) or by a midpoint run when the scores are
              within MERGE_CLOSE: the survivor inherits transplants, the loser's tested set and finer steps
  funnel      at FUNNEL times keep the holder of the best unconditionally, then the largest recent gains: 6 -> 3 -> 1;
              a freed slot (freeze, merge, funnel) is refilled from the scout's finds (rolling domains)
  final       the last FINAL_S on the single best point among live and frozen domains: halving steps and all pairs,
              no filler; the scout takes the spare workers

Riders only in the global waves and the scout (each cell's threshold sweep). A climbing head is exactly the move
its domain asked for, with two cheap siblings; every move goes in the throughput lane so moves sharing a traversal
pack. No barrier anywhere. Main submits the local best on every improvement; the engine never submits.
"""
from __future__ import annotations

import asyncio
import logging
import random
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Set, Tuple

from common.params import ParamRanges
from common.types import FULL, Config
from search.axes import AxisBook, far_axes, latin_hypercube
from search.base import Strategy
from search.client import EvalClient, Outcome
from search.domain import ACCEPT_MIN as ACCEPT_MIN_MERGE, Domain, Move
from search.engine import FAMILY_SEED, OFFICIAL_SEED, settled_capacity
from search.families import genotyping_riders
from search.ranker import Ranker
from search.space import GENOTYPING_PARAMS, Space, bam_caps

log = logging.getLogger(__name__)

OPENING_S = 10 * 60
WAVE2_AT = 5 * 60                  # the second global wave starts here (or once wave 1 is WAVE_BACK back)
WAVE_BACK = 0.8                    # a wave counts as back at this share: the slow tail must not idle the fleet
WAVE2_IDLE = 0.25                  # ... or wave 2 goes out once this share of the workers is idle with half of wave 1 back
TICK_S = 2.0                       # the loop wakes at least this often: timed steps must not wait for a result
OPENING_BACK = 0.6                 # at OPENING_S the domains open once this share of wave 2 is back
EARLY_LEAD = 0.1                   # a wave-1 leader this far above every far rival founds its domain at once
N_DOMAINS = 6
DIVISION_AXES = 3
GOOD_BAND = 1.0                    # opening configs within this of the best define the division axes
WORTH_A_HILL = 0.5                 # a domain is never founded this far below the best
MIN_FAR = 4                        # domains closer than this (live axes > 10% apart) are one hill: merge
MERGE_CLOSE = 0.02                 # ... or this close in score, if one midpoint run shows no valley between them
MIDPOINT_RETEST_S = 10 * 60        # a pair found distinct is not retested for this long
BATCH_MAX, BATCH_MIN = 300, 30
FUNNEL: Tuple[Tuple[float, int], ...] = ((28 * 60, 3), (60 * 60, 1))
FINAL_S = 10 * 60
LATE_BIRTH_S = FINAL_S + 10 * 60   # no new domain this close to the end: it could not be climbed in time
SCOUT_SHARE = 0.15
SCOUT_KAPPA = 1.0                  # the scout's upper confidence bound: mean + kappa x spread
RIDERS, RIDERS_SMALL = 8, 6
SMALL_FLEET = 64


@dataclass
class Unit:
    role: str                        # opening | scout | move | midpoint | final
    head: Config
    domain: Optional[int] = None
    move: Optional[Move] = None
    wave: int = 0


def searcher(ranges: ParamRanges, workers: int, *, tool: str = "gatk", seed: int = 0, tag: str = "engine3") -> Strategy:
    rng = random.Random(seed)

    async def run(client: EvalClient) -> None:
        r = client.round
        caps = bam_caps(r.bam, r.chrom, r.start, r.end) if r.bam else {}
        space = Space(ranges, tool, caps=caps)
        if caps:
            log.info("engine3: caps from the round's BAM: %s", caps)
        book = AxisBook(space)
        ranker = Ranker(space, seed=seed)
        clock = client.clock
        t0 = clock()
        cap = await settled_capacity(client, workers)
        small = cap < SMALL_FLEET
        riders = RIDERS_SMALL if small else RIDERS
        n_domains_max = max(1, min(N_DOMAINS, cap // 8))
        virtual = hasattr(client.executor, "now")
        pending: Dict[asyncio.Future, Tuple[Unit, bool]] = {}
        proposed: Set[str] = set()
        domains: Dict[int, Domain] = {}
        next_id = [0]
        wave_heads = [0, 0]                                      # heads sent per wave
        wave_back = [0, 0]
        wave2_sent = [False]
        early_done = [False]
        opened = [False]
        candidates: Dict[str, Tuple[Config, float]] = {}
        was_centre: Set[str] = set()
        division: List[str] = []
        batch_tests = [BATCH_MAX]
        funnel_done: Set[float] = set()
        final_sent = [False]
        scouts_inflight = [0]
        scout_queue: List[Config] = []
        fitting: List[Optional[asyncio.Task]] = [None]
        midpoints: Dict[str, Tuple[int, int]] = {}               # midpoint key -> the pair under test
        distinct_until: Dict[Tuple[int, int], float] = {}
        last_log = [t0]
        wake = asyncio.Event()
        client.on_result.append(lambda _o: wake.set())

        def elapsed() -> float:
            return clock() - t0

        def left() -> Optional[float]:
            return client.time_left()

        def fresh(c: Optional[Config]) -> bool:
            return c is not None and c.key not in proposed and (c.key, FULL) not in client.latest

        def live() -> List[Domain]:
            return [d for d in domains.values() if not d.retired and not d.frozen]

        def standing() -> List[Domain]:
            return [d for d in domains.values() if not d.retired]

        genotyping = [a for a in GENOTYPING_PARAMS if a in space.names]

        def live_names() -> List[str]:
            return book.live(clock()) + genotyping          # the genotyping axes are never parked: steps on them are free

        def owner(c: Config) -> int:
            """The live domain whose centre is nearest on the division axes (-1: none)."""
            ld = live()
            if not ld or not division:
                return -1
            def dist(d: Domain) -> float:
                return sum((space._to_unit(space.specs[a], c.params.get(a, space.specs[a].default))
                            - space._to_unit(space.specs[a], d.centre.params.get(a, space.specs[a].default))) ** 2
                           for a in division if not (space.specs[a].allowed or space.specs[a].type == "bool")) \
                    + sum(float(c.params.get(a) != d.centre.params.get(a))
                          for a in division if (space.specs[a].allowed or space.specs[a].type == "bool"))
            return min(ld, key=dist).id

        def unowned(c: Config) -> bool:
            """Outside every standing domain's hill: far from each centre on MIN_FAR live axes."""
            names = live_names()
            return all(far_axes(space, names, c, d.centre) >= MIN_FAR for d in standing())

        def keep_now() -> int:
            k = n_domains_max
            for at, keep in FUNNEL:
                if elapsed() >= at:
                    k = min(k, keep)
            return k

        # ---- proposing -------------------------------------------------------------------------------------------
        def propose(unit: Unit, rider_configs: Sequence[Config], priority: int, lane: str) -> int:
            batch = [unit.head] + [c for c in rider_configs if fresh(c)]
            if not fresh(unit.head) and unit.role != "move":
                batch = batch[1:]
            if not batch:
                return 0
            lanes = [lane] + ["throughput"] * (len(batch) - 1)
            kind = unit.move.kind if unit.move is not None else unit.role
            futs = client.propose_many(batch, priority=priority, tag=f"{tag}:{kind}", lanes=lanes)
            for c, f in zip(batch, futs):
                proposed.add(c.key)
                pending[f] = (unit, c is unit.head)
            return len(futs)

        def propose_head(role: str, head: Config, priority: int, lane: str, domain: Optional[Domain] = None,
                         move: Optional[Move] = None, wave: int = 0) -> bool:
            """Global points (the waves, the scout) carry the threshold riders; a climbing head is exactly the move
            the domain asked for, plus its two cheap siblings."""
            unit = Unit(role, head, domain.id if domain else None, move, wave)
            if role in ("opening", "scout"):
                kids = genotyping_riders(space, head, riders, fresh)
            elif domain is not None and move is not None:
                kids = domain.siblings(move, rng)
            else:
                kids = []
            return propose(unit, kids, priority, lane) > 0

        def seeds() -> List[Config]:
            out: List[Config] = []
            for params in ({}, OFFICIAL_SEED, FAMILY_SEED):
                c = space._assemble(dict(params))
                if c is not None and all(c.key != s.key for s in out):
                    out.append(c)
            for raw in client.store.portfolio(limit=max(4, cap // 20), tool=tool):
                c = space._assemble(dict(raw.params))
                if c is not None and all(c.key != s.key for s in out):
                    out.append(c)
            return out

        def wave1() -> None:
            sd = seeds()
            cells = [c for c in latin_hypercube(space, rng, max(0, cap - len(sd))) if all(c.key != s.key for s in sd)]
            for i, c in enumerate(sd + cells):
                if propose_head("opening", c, 5 if i < len(sd) else 4, "throughput", wave=0):
                    wave_heads[0] += 1
            batch_tests[0] = max(BATCH_MIN, min(BATCH_MAX, round(2 * cap / n_domains_max)))
            log.info("engine3: capacity %d -> up to %d domains, %d riders, batch %d tests, scout share %.0f%%; "
                     "wave 1: %d seeds + %d hypercube cells", cap, n_domains_max, riders, batch_tests[0],
                     SCOUT_SHARE * 100, len(sd), len(cells))

        def model_pick(n: int, kappa: float) -> List[Config]:
            """`n` fresh global candidates: hypercube points and variations of the best configs, ranked by the model's
            upper confidence bound (unranked while the model is not fitted)."""
            if n <= 0:
                return []
            pool = [c for c in latin_hypercube(space, rng, 3 * n) if fresh(c) and unowned(c)]
            bests = [c for c, _ in sorted(candidates.values(), key=lambda t: -t[1])[:8]]
            tries = 0
            while bests and len(pool) < 4 * n and tries < 4 * n:
                tries += 1
                c = space.neighbour(rng.choice(bests), rng, n_params=rng.randint(2, 4), scale=0.3)
                if c is not None and fresh(c) and unowned(c) and all(c.key != p.key for p in pool):
                    pool.append(c)
            if not pool:
                return []
            if ranker.ready:
                mean, spread = ranker.predict(pool)
                top = max((s for _, s in candidates.values()), default=0.0)
                order = sorted(range(len(pool)), key=lambda i: -(mean[i] - top + kappa * spread[i]))
                pool = [pool[i] for i in order]
            return pool[:n]

        def wave2() -> None:
            wave2_sent[0] = True
            t_a = time.monotonic()
            idle = client.idle_slots()
            n = max(1, min(cap, idle if idle is not None else cap))
            if ranker.ready:
                half = n // 2
                picks = model_pick(half, kappa=0.0) + model_pick(n - half, kappa=SCOUT_KAPPA * 3)
            else:
                picks = model_pick(n, kappa=SCOUT_KAPPA)
            t_b = time.monotonic()
            seen = set()
            sent = 0
            for c in picks:
                if c.key in seen:
                    continue
                seen.add(c.key)
                if propose_head("opening", c, 4, "throughput", wave=1):
                    wave_heads[1] += 1
                    sent += 1
            log.info("engine3: wave 2 at t+%.1fm: %d points (%s) with wave 1 %d/%d back; picking %.1f s, sending %.1f s",
                     elapsed() / 60, sent, "model-ranked" if ranker.ready else "hypercube", wave_back[0], wave_heads[0],
                     t_b - t_a, time.monotonic() - t_b)

        def choose_division() -> List[str]:
            """The live axes with the widest spread among the good opening configs."""
            ranked = sorted(candidates.values(), key=lambda t: -t[1])
            if not ranked:
                return []
            top = ranked[0][1]
            good = [c for c, s in ranked if s >= top - GOOD_BAND] or [ranked[0][0]]
            spread = {}
            for a in live_names():
                s = space.specs[a]
                if s.allowed or s.type == "bool":
                    spread[a] = len({c.params.get(a) for c in good}) / max(2, len(s.allowed or (True, False)))
                else:
                    us = [space._to_unit(s, c.params.get(a, s.default)) for c in good]
                    spread[a] = max(us) - min(us)
            return sorted(spread, key=lambda a: -spread[a])[:DIVISION_AXES]

        def early_leader() -> None:
            """Wave 1 mostly back: a config EARLY_LEAD above every far rival climbs now instead of at OPENING_S."""
            ranked = sorted(candidates.values(), key=lambda t: -t[1])
            if not ranked or live():
                return
            c, s = ranked[0]
            names = live_names()
            rivals = [sx for x, sx in ranked[1:] if far_axes(space, names, c, x) >= MIN_FAR]
            if not rivals or s - max(rivals) >= EARLY_LEAD:
                division[:] = choose_division()
                birth(c, s, clock(), "a clear wave-1 leader")

        def open_domains() -> None:
            opened[0] = True
            division[:] = choose_division()
            fill_slots("the opening")
            log.info("engine3: opening done at t+%.1fm: waves %d/%d and %d/%d back, division axes %s, %d domains",
                     elapsed() / 60, wave_back[0], wave_heads[0], wave_back[1], wave_heads[1], division, len(live()))

        def fill_slots(why: str) -> int:
            """Free domain slots (keep_now minus live) filled from the best candidates in unowned territory above the
            hill margin: the opening's starts, and later the scout's finds (rolling domains)."""
            ranked = sorted(candidates.values(), key=lambda t: -t[1])
            if not ranked:
                return 0
            t = left()
            if t is not None and t < LATE_BIRTH_S:
                return 0
            top = ranked[0][1]
            now = clock()
            n = 0
            for c, s in ranked:
                if len(live()) >= keep_now() or s < top - WORTH_A_HILL:
                    break
                if c.key in was_centre or not unowned(c):
                    continue
                if division and any(far_axes(space, division, c, d.centre) < 1 for d in live()):
                    continue
                birth(c, s, now, why)
                n += 1
            return n

        def birth(c: Config, s: float, now: float, why: str) -> Domain:
            d = Domain(next_id[0], space, book, c, s, now, owner)
            next_id[0] += 1
            domains[d.id] = d
            was_centre.add(c.key)
            log.info("engine3: domain %d born at %.4f (%s) t+%.1fm", d.id, s, why, elapsed() / 60)
            return d

        def final_round() -> None:
            final_sent[0] = True
            st = standing()
            ld = max(st, key=lambda d: d.best) if st else None      # the best point, live or frozen
            for d in domains.values():
                if d is not ld:
                    d.retired = True
            if ld is None:
                return
            n = ld.start_final(live_names(), rng, cap)
            log.info("engine3: final at t+%.1fm on domain %d (%.4f): %d moves", elapsed() / 60, ld.id, ld.best, n)

        def scout_target() -> int:
            return max(1, int(cap * SCOUT_SHARE)) if (opened[0] or live()) else 0

        def send_scout(idle: Optional[int]) -> bool:
            if not scout_queue:
                scout_queue.extend(model_pick(max(1, min(64, idle or 64)), SCOUT_KAPPA))
            while scout_queue:
                c = scout_queue.pop(0)
                if fresh(c) and propose_head("scout", c, 1, "throughput"):
                    scouts_inflight[0] += 1
                    return True
            return False

        def fill() -> int:
            n = 0
            while not client.exhausted:
                idle = client.idle_slots()
                if idle is not None and idle <= 0:
                    break
                if idle is None and len(pending) >= cap:
                    break
                t = left()
                if t is not None and t < FINAL_S and opened[0] and not final_sent[0]:
                    final_round()
                    n += 1
                    continue
                ld = live()
                if not wave2_sent[0] and not ld:
                    break                                        # wave 1 is out; wave 2 comes from its results
                if ld and scouts_inflight[0] < scout_target() and send_scout(idle):
                    n += 1
                    continue
                if ld:
                    target = max(1.0, cap * (1 - SCOUT_SHARE) / len(ld))
                    for d in ld:
                        if not d.pool and d.inflight < target:
                            d.refill(live_names(), rng, int(target))
                    ready = [d for d in ld if d.pool]
                    under = [d for d in ready if d.inflight < target]
                    picked = min(under, key=lambda d: d.inflight / target) if under else \
                        (min(ready, key=lambda d: d.inflight) if ready else None)
                    if picked is not None:
                        move = picked.take()
                        if move is not None and propose_head("move", move.config, 3, "throughput", picked, move):
                            n += 1
                        elif move is not None:
                            picked.inflight = max(0, picked.inflight - 1)
                        continue
                if wave2_sent[0] and send_scout(idle):           # idle capacity never waits: the scout takes it
                    n += 1
                    continue
                break
            return n

        # ---- results -------------------------------------------------------------------------------------------
        def note_result(unit: Unit, is_head: bool, o: Outcome) -> None:
            now = clock()
            score = o.score if o.ok else None
            failed = (not o.ok) or score is None or score <= 0
            rr = o.result
            ranker.add(o.config, rr.core if (rr and not failed) else None,
                       (rr.region_fp / rr.n_targets) if (rr and not failed and rr.n_targets) else None, failed)
            if not failed:
                k = o.config.key
                if k not in candidates or score > candidates[k][1]:
                    candidates[k] = (o.config, score)
                if len(candidates) > 600:
                    for key in sorted(candidates, key=lambda x: -candidates[x][1])[300:]:
                        del candidates[key]
            if unit.role in ("opening", "scout"):
                if is_head:
                    if unit.role == "opening":
                        wave_back[unit.wave] += 1
                    else:
                        scouts_inflight[0] = max(0, scouts_inflight[0] - 1)
                if opened[0] and not failed:
                    adopt_if_better(o.config, score, now)
                return
            if unit.role == "midpoint":
                judge_midpoint(unit.head.key, score if not failed else None)
                return
            d = domains.get(unit.domain) if unit.domain is not None else None
            if d is None or unit.move is None:
                return
            if d.retired:
                if not failed:
                    adopt_if_better(o.config, score, now)
                return
            was_frozen = d.frozen
            moved = d.consider(unit.move, is_head, o.config, score, failed, now)
            if moved:
                log.info("engine3: domain %d -> %.4f via %s %s t+%.1fm", d.id, d.best, unit.move.kind,
                         ",".join(unit.move.axes), elapsed() / 60)
            if d.bracket is not None and d.bracket[1] >= 0.8 * d.bracket[0]:
                gain = d.close_batch()
                log.info("engine3: domain %d batch %d closed at %.4f (gain %+.4f, %s) t+%.1fm", d.id, d.batch, d.best,
                         gain, d.mode, elapsed() / 60)
                merges()
            elif d.bracket is None and d.completed >= batch_tests[0] and d.phase != "final":
                sent = d.start_bracket(live_names())
                log.info("engine3: domain %d bracket of %d after %d tests t+%.1fm", d.id, sent, d.completed, elapsed() / 60)
            if not d.frozen and d.phase != "final" and d.stalled(now):
                d.frozen = True
            if d.frozen and not was_frozen and d.phase != "final":
                log.info("engine3: domain %d frozen at %.4f (%s) t+%.1fm", d.id, d.best,
                         "exhausted at the floor" if d.exhausted_at_floor else "no gain above the resolution for 12 min",
                         elapsed() / 60)

        def adopt_if_better(c: Config, s: float, now: float) -> None:
            """A global point (a wave cell, a scout) that beats the domain owning it, on its hill."""
            oid = owner(c)
            d = domains.get(oid)
            if d is not None and not d.retired and s > d.best + 0.005 and far_axes(space, live_names(), c, d.centre) < MIN_FAR:
                d._recentre(c, s, now, ())
                log.info("engine3: domain %d -> %.4f, a global point on its hill t+%.1fm", d.id, s, elapsed() / 60)

        def merge_into(survivor: Domain, loser: Domain, why: str) -> None:
            loser.retired = True
            n = survivor.inherit(loser, live_names())
            log.info("engine3: domain %d (%.4f) merged into domain %d (%.4f), %s: %d transplants t+%.1fm", loser.id,
                     loser.best, survivor.id, survivor.best, why, n, elapsed() / 60)

        def merges() -> None:
            """Converged domains merge: centres fewer than MIN_FAR live axes apart at once; centres within
            MERGE_CLOSE in score after one midpoint run shows no valley. Sharing happens only here."""
            ld = live()
            for a, b in [(a, b) for a in ld for b in ld if a.id < b.id]:
                if a.retired or b.retired or a.frozen or b.frozen:
                    continue
                worse, better = (a, b) if a.best < b.best else (b, a)
                if far_axes(space, live_names(), a.centre, b.centre) < MIN_FAR:
                    merge_into(better, worse, "one hill by distance")
                    continue
                pair = (a.id, b.id)
                if abs(a.best - b.best) <= MERGE_CLOSE and distinct_until.get(pair, 0.0) <= clock() \
                        and pair not in midpoints.values():
                    test_midpoint(a, b)

        def test_midpoint(a: Domain, b: Domain) -> None:
            params = dict(a.centre.params if a.best >= b.best else b.centre.params)
            for name in live_names():
                s = space.specs[name]
                if s.allowed or s.type == "bool":
                    continue
                u = (space._to_unit(s, a.centre.params.get(name, s.default)) + space._to_unit(s, b.centre.params.get(name, s.default))) / 2
                params[name] = space._from_unit(s, u)
            mid = space._assemble(params)
            if mid is None or mid.key in midpoints:
                return
            midpoints[mid.key] = (a.id, b.id)
            if not propose_head("midpoint", mid, 4, "throughput"):
                midpoints.pop(mid.key, None)

        def judge_midpoint(key: str, score: Optional[float]) -> None:
            pair = midpoints.pop(key, None)
            if pair is None:
                return
            a, b = domains.get(pair[0]), domains.get(pair[1])
            if a is None or b is None or a.retired or b.retired:
                return
            if score is not None and score >= min(a.best, b.best) - ACCEPT_MIN_MERGE:
                worse, better = (a, b) if a.best < b.best else (b, a)
                merge_into(better, worse, f"one hill by midpoint {score:.4f}")
            else:
                distinct_until[pair] = clock() + MIDPOINT_RETEST_S

        def funnel() -> None:
            for at, keep in FUNNEL:
                if at in funnel_done or elapsed() < at:
                    continue
                funnel_done.add(at)
                merges()
                ld = live()
                if not ld:
                    continue
                holder = max(ld, key=lambda d: d.best)             # the holder of the best is always kept
                rest = sorted((d for d in ld if d is not holder), key=lambda d: (-d.gain_recent(), -d.best))
                kept = [holder] + rest[:max(0, keep - 1)]
                for d in rest[max(0, keep - 1):]:
                    d.retired = True
                    log.info("engine3: domain %d (%.4f, gain %+.4f) dropped by the funnel at t+%.1fm", d.id, d.best,
                             d.gain_recent(), elapsed() / 60)
                log.info("engine3: funnel at t+%.1fm keeps %s", elapsed() / 60,
                         " ".join(f"{d.id}:{d.best:.4f}/{d.gain_recent():+.4f}" for d in kept))

        def maybe_fit() -> None:
            if ranker.stale() and (fitting[0] is None or fitting[0].done()):
                fitting[0] = asyncio.create_task(asyncio.to_thread(ranker.fit))

        # ---- the loop --------------------------------------------------------------------------------------------
        wave1()
        while True:
            if pending:
                done, _ = await asyncio.wait(list(pending), timeout=TICK_S, return_when=asyncio.FIRST_COMPLETED)
                t_r = time.monotonic()
                for fut in done:
                    unit, is_head = pending.pop(fut)
                    try:
                        note_result(unit, is_head, fut.result())
                    except Exception:  # noqa: BLE001
                        log.exception("engine3: result handling failed")
                maybe_fit()
                if time.monotonic() - t_r > 5:
                    log.info("engine3: slow: %d results handled in %.1f s at t+%.1fm", len(done), time.monotonic() - t_r, elapsed() / 60)
            w1_back = bool(wave_heads[0]) and wave_back[0] >= WAVE_BACK * wave_heads[0]
            idle_now = client.idle_slots()
            w1_idle = bool(wave_heads[0]) and wave_back[0] >= 0.5 * wave_heads[0] and idle_now is not None \
                and idle_now >= WAVE2_IDLE * cap                 # the slow tail is holding the fleet: go now
            if not wave2_sent[0] and (elapsed() >= WAVE2_AT or w1_back or w1_idle):
                if not early_done[0] and w1_back:
                    early_done[0] = True
                    early_leader()
                wave2()
            if not opened[0] and wave2_sent[0] and (
                    wave_back[1] >= WAVE_BACK * max(1, wave_heads[1])
                    or (elapsed() >= OPENING_S and wave_back[1] >= OPENING_BACK * max(1, wave_heads[1]))
                    or not pending):
                open_domains()
            if opened[0]:
                funnel()
                if len(live()) < keep_now():
                    fill_slots("a scout find")
                t = left()
                if not live() and not final_sent[0] and standing() and (t is None or t < LATE_BIRTH_S):
                    final_round()                                # nothing left to climb this late: the final, early
            t_f = time.monotonic()
            n = fill()
            if time.monotonic() - t_f > 5:
                log.info("engine3: slow: fill of %d proposals took %.1f s at t+%.1fm", n, time.monotonic() - t_f, elapsed() / 60)
            if clock() - last_log[0] > 120:
                last_log[0] = clock()
                ld = live()
                frozen = [d for d in domains.values() if d.frozen and not d.retired]
                log.info("engine3: t+%.1fm domains %s | frozen %s | scouts %d | candidates %d | pending %d", elapsed() / 60,
                         " ".join(f"{d.id}:{d.best:.3f}/{d.mode if d.phase != 'final' else 'final'}/{d.inflight}" for d in ld),
                         " ".join(f"{d.id}:{d.best:.3f}" for d in frozen), scouts_inflight[0], len(candidates), len(pending))
                client.notes["engine"] = {
                    "basins": [{"id": d.id, "best": round(d.best, 4), "generation": d.generation, "inflight": d.inflight,
                                "phase": d.phase, "batch": d.batch, "born_min": round((d.born - t0) / 60, 1)} for d in ld],
                    "division": list(division), "scouts": scouts_inflight[0], "candidates": len(candidates),
                    "pending": len(pending), "caps": caps, "dead_axes": sorted(set(book.names) - set(book.live(clock())))}
            if client.exhausted and not pending:
                return
            if not pending and n == 0:
                t = left()
                if t is not None and t <= 0:
                    return
                if virtual:
                    await asyncio.sleep(0)
                else:
                    wake.clear()
                    try:
                        await asyncio.wait_for(wake.wait(), timeout=1.0)
                    except asyncio.TimeoutError:
                        pass

    return run
