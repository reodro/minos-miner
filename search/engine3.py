"""Engine v3: isolated iterated local search with racing. Owner: search.

The objective is discrete and expensive, the fleet is wide and a round allows ~15-20 dependent steps, so:

  opening     0-8 min: seeds + a designed wave + random cells over the whole space, 8 genotyping riders each.
              From it: the DIVISION axes (the 3 live axes with the widest spread among the good configs) and up to
              N_DOMAINS domains, one at each best config pairwise far on the division axes and above the hill margin;
              dead axes are parked by the axis book as evidence comes in
  domains     each an isolated iterated local search (domain.py): neighbourhood 1 at effective steps, neighbourhood 2
              pairs, sideways moves on plateaus, shakes of 2/3/5/8 axes when exhausted. A config belongs to the
              domain whose centre is nearest on the division axes; a domain never proposes into another's territory
  batches     per domain, after BATCH completed tests: the bracket (4 values across every live axis + all pairs of the
              5 strongest axes), then merges of domains that converged (fewer than MIN_FAR live axes apart)
  funnel      at FUNNEL times keep the domains with the largest gain over their last two batches: 6 -> 3 -> 1
  final       the last FINAL_S on the single best point: halving steps, all pairs of the top axes, small shakes

Riders only in the opening (each cell's threshold sweep). A climbing head is exactly the move its domain asked for;
the genotyping axes are ordinary axes of the neighbourhoods, and the evaluator packs the genotyping steps of one
centre into one traversal by itself.

No barrier anywhere: a domain closes its own batch when its own count is reached. Main submits the local best on
every improvement; the engine never submits.
"""
from __future__ import annotations

import asyncio
import logging
import random
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Set, Tuple

from common.params import ParamRanges
from common.types import FULL, Config
from search.axes import AxisBook, design_wave, far_axes
from search.base import Strategy
from search.client import EvalClient, Outcome
from search.domain import ACCEPT_MIN as ACCEPT_MIN_MERGE, Domain, Move
from search.engine import FAMILY_SEED, OFFICIAL_SEED, settled_capacity
from search.families import genotyping_riders
from search.space import GENOTYPING_PARAMS
from search.space import Space, bam_caps

log = logging.getLogger(__name__)

OPENING_S = 8 * 60
OPENING_BACK = 0.6                 # ... and at least this share of the opening wave back
OPENING_EARLY = 0.8                # this share back ends the opening early: the slow tail must not idle the fleet
N_DOMAINS = 6
DIVISION_AXES = 3
GOOD_BAND = 1.0                    # opening configs within this of the best define the division axes
WORTH_A_HILL = 0.5                 # a domain is never founded this far below the best
MIN_FAR = 4                        # domains closer than this (live axes > 10% apart) are one hill: merge
MERGE_CLOSE = 0.02                 # ... or this close in score, if one midpoint run shows no valley between them
MIDPOINT_RETEST_S = 10 * 60        # a pair found distinct is not retested for this long
BATCH_MAX, BATCH_MIN = 300, 30
FUNNEL: Tuple[Tuple[float, int], ...] = ((24 * 60, 3), (60 * 60, 1))
FINAL_S = 10 * 60
RIDERS, RIDERS_SMALL = 8, 6
SMALL_FLEET = 64
RANDOM_SHARE = 0.2                 # of the opening wave: uniform random cells beside the designed ones


@dataclass
class Unit:
    role: str                        # opening | move | final
    head: Config
    domain: Optional[int] = None
    move: Optional[Move] = None


def searcher(ranges: ParamRanges, workers: int, *, tool: str = "gatk", seed: int = 0, tag: str = "engine3") -> Strategy:
    rng = random.Random(seed)

    async def run(client: EvalClient) -> None:
        r = client.round
        caps = bam_caps(r.bam, r.chrom, r.start, r.end) if r.bam else {}
        space = Space(ranges, tool, caps=caps)
        if caps:
            log.info("engine3: caps from the round's BAM: %s", caps)
        book = AxisBook(space)
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
        opening_heads = [0]
        opening_back = [0]
        opened = [False]
        candidates: Dict[str, Tuple[Config, float]] = {}
        division: List[str] = []
        batch_tests = [BATCH_MAX]
        funnel_done: Set[float] = set()
        final_sent = [False]
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

        midpoints: Dict[str, Tuple[int, int]] = {}                 # midpoint key -> the pair under test
        distinct_until: Dict[Tuple[int, int], float] = {}

        def leader() -> Optional[Domain]:
            ld = live()
            return max(ld, key=lambda d: d.best) if ld else None

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
                         move: Optional[Move] = None) -> bool:
            """Only the opening wave carries riders (each cell's threshold sweep); a climbing head is exactly the
            move the domain asked for - its genotyping steps are moves of their own, packed by the evaluator."""
            unit = Unit(role, head, domain.id if domain else None, move)
            if role == "opening":
                kids = genotyping_riders(space, head, riders, fresh)
            elif domain is not None and move is not None:
                kids = domain.siblings(move, rng)                # one genotyping and one PairHMM variation of the move
            else:
                kids = []
            return propose(unit, kids, priority, lane) > 0

        def opening() -> None:
            seeds: List[Config] = []
            for params in ({}, OFFICIAL_SEED, FAMILY_SEED):
                c = space._assemble(dict(params))
                if c is not None and all(c.key != s.key for s in seeds):
                    seeds.append(c)
            for raw in client.store.portfolio(limit=max(4, cap // 20), tool=tool):
                c = space._assemble(dict(raw.params))
                if c is not None and all(c.key != s.key for s in seeds):
                    seeds.append(c)
            n_wave = max(0, cap - len(seeds))
            n_random = int(n_wave * RANDOM_SHARE)
            cells = design_wave(space, rng, n_wave - n_random, seeds[:3])
            seen = {c.key for c in seeds} | {c.key for c in cells}
            randoms: List[Config] = []
            while len(randoms) < n_random:
                c = space.random(rng)
                if c.key not in seen:
                    seen.add(c.key)
                    randoms.append(c)
            heads = seeds + cells + randoms
            for i, c in enumerate(heads):
                if propose_head("opening", c, 5 if i < len(seeds) else 4, "throughput"):
                    opening_heads[0] += 1
            batch_tests[0] = max(BATCH_MIN, min(BATCH_MAX, round(2 * cap / n_domains_max)))
            log.info("engine3: capacity %d -> up to %d domains, %d riders, batch %d tests; opening %d seeds + %d design + %d random",
                     cap, n_domains_max, riders, batch_tests[0], len(seeds), len(cells), len(randoms))

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

        def open_domains() -> None:
            opened[0] = True
            division[:] = choose_division()
            ranked = sorted(candidates.values(), key=lambda t: -t[1])
            if not ranked:
                return
            top = ranked[0][1]
            now = clock()
            for c, s in ranked:
                if len(live()) >= n_domains_max or s < top - WORTH_A_HILL:
                    break
                if all(far_axes(space, division, c, d.centre) >= 1 for d in live()):
                    birth(c, s, now)
            log.info("engine3: opening done at t+%.1fm: %d heads back of %d, division axes %s, %d domains",
                     elapsed() / 60, opening_back[0], opening_heads[0], division, len(live()))

        def birth(c: Config, s: float, now: float) -> Domain:
            d = Domain(next_id[0], space, book, c, s, now, owner)
            next_id[0] += 1
            domains[d.id] = d
            log.info("engine3: domain %d born at %.4f t+%.1fm", d.id, s, elapsed() / 60)
            return d

        def final_round() -> None:
            final_sent[0] = True
            ld = leader()
            if ld is None:                                       # every domain frozen: the best of them finishes
                frozen = [d for d in domains.values() if not d.retired]
                ld = max(frozen, key=lambda d: d.best) if frozen else None
            for d in domains.values():
                if d is not ld:
                    d.retired = True
            if ld is None:
                return
            n = ld.start_final(live_names(), rng, cap)
            log.info("engine3: final at t+%.1fm on domain %d (%.4f): %d moves", elapsed() / 60, ld.id, ld.best, n)

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
                if not ld:
                    if opened[0] and not final_sent[0] and any(not d.retired for d in domains.values()):
                        final_round()                            # every domain frozen: the final on the best, early
                        n += 1
                        continue
                    break
                target = max(1.0, cap / len(ld))
                for d in ld:
                    if not d.pool and d.inflight < target:
                        d.refill(live_names(), rng, int(target))
                ready = [d for d in ld if d.pool]
                under = [d for d in ready if d.inflight < target]
                picked = min(under, key=lambda d: d.inflight / target) if under else (min(ready, key=lambda d: d.inflight) if ready else None)
                if picked is None:
                    break
                move = picked.take()
                if move is None:
                    continue
                # every move in the throughput lane: moves of one centre that share a traversal (genotyping, PairHMM,
                # assembly, activity steps) then pack into one run; a latency job never rides in another's batch
                if propose_head("move", move.config, 3, "throughput", picked, move):
                    n += 1
                else:
                    picked.inflight = max(0, picked.inflight - 1)
            return n

        # ---- results -------------------------------------------------------------------------------------------
        def note_result(unit: Unit, is_head: bool, o: Outcome) -> None:
            now = clock()
            score = o.score if o.ok else None
            failed = (not o.ok) or score is None or score <= 0
            if not failed:
                k = o.config.key
                if k not in candidates or score > candidates[k][1]:
                    candidates[k] = (o.config, score)
                if len(candidates) > 600:
                    for key in sorted(candidates, key=lambda x: -candidates[x][1])[300:]:
                        del candidates[key]
            if unit.role == "opening":
                if is_head:
                    opening_back[0] += 1
                if opened[0] and not failed:
                    adopt_if_better(o.config, score, now)
                return
            if unit.role == "midpoint":
                judge_midpoint(unit.head.key, score if not failed else None)
                return
            d = domains.get(unit.domain) if unit.domain is not None else None
            if d is None or unit.move is None:
                if d is not None and not d.retired and not failed and score > d.best + 0.005:
                    d._recentre(o.config, score, now, ())          # a dense final rider on top
                    log.info("engine3: domain %d -> %.4f via final rider t+%.1fm", d.id, score, elapsed() / 60)
                return
            if d.retired:
                if not failed:
                    adopt_if_better(o.config, score, now)
                return
            moved = d.consider(unit.move, is_head, o.config, score, failed, now)
            if moved:
                log.info("engine3: domain %d -> %.4f via %s %s t+%.1fm", d.id, d.best, unit.move.kind,
                         ",".join(unit.move.axes), elapsed() / 60)
            if d.bracket is not None and d.bracket[1] >= 0.8 * d.bracket[0]:
                gain = d.close_batch()
                log.info("engine3: domain %d batch %d closed at %.4f (gain %+.4f, %s) t+%.1fm", d.id, d.batch, d.best,
                         gain, d.mode, elapsed() / 60)
                if d.frozen:
                    log.info("engine3: domain %d frozen at %.4f after %d batches without a gain t+%.1fm", d.id, d.best,
                             len(d.gains), elapsed() / 60)
                merges()
            elif d.bracket is None and d.completed >= batch_tests[0] and d.phase != "final":
                sent = d.start_bracket(live_names())
                log.info("engine3: domain %d bracket of %d after %d tests t+%.1fm", d.id, sent, d.completed, elapsed() / 60)

        def adopt_if_better(c: Config, s: float, now: float) -> None:
            """A late opening result (or a retired domain's straggler) that beats the domain owning it."""
            oid = owner(c)
            d = domains.get(oid)
            if d is not None and not d.retired and s > d.best + 0.005 and far_axes(space, live_names(), c, d.centre) < MIN_FAR:
                d._recentre(c, s, now, ())
                log.info("engine3: domain %d -> %.4f, an opening config on its hill t+%.1fm", d.id, s, elapsed() / 60)

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
                ld = sorted(live(), key=lambda d: (-d.gain_recent(), -d.best))
                for d in ld[keep:]:
                    d.retired = True
                    log.info("engine3: domain %d (%.4f, gain %+.4f) dropped by the funnel at t+%.1fm", d.id, d.best,
                             d.gain_recent(), elapsed() / 60)
                log.info("engine3: funnel at t+%.1fm keeps %s", elapsed() / 60,
                         " ".join(f"{d.id}:{d.best:.4f}/{d.gain_recent():+.4f}" for d in ld[:keep]))

        # ---- the loop --------------------------------------------------------------------------------------------
        opening()
        while True:
            if pending:
                done, _ = await asyncio.wait(list(pending), return_when=asyncio.FIRST_COMPLETED)
                for fut in done:
                    unit, is_head = pending.pop(fut)
                    try:
                        note_result(unit, is_head, fut.result())
                    except Exception:  # noqa: BLE001
                        log.exception("engine3: result handling failed")
            if not opened[0] and opening_heads[0] and (
                    opening_back[0] >= OPENING_EARLY * opening_heads[0]
                    or (elapsed() >= OPENING_S and opening_back[0] >= OPENING_BACK * opening_heads[0])
                    or (not pending)):
                open_domains()
            if opened[0]:
                funnel()
            n = fill()
            if clock() - last_log[0] > 120:
                last_log[0] = clock()
                ld = live()
                frozen = [d for d in domains.values() if d.frozen and not d.retired]
                log.info("engine3: t+%.1fm domains %s | frozen %s | candidates %d | pending %d", elapsed() / 60,
                         " ".join(f"{d.id}:{d.best:.3f}/{d.mode if d.phase != 'final' else 'final'}/{d.inflight}" for d in ld),
                         " ".join(f"{d.id}:{d.best:.3f}" for d in frozen), len(candidates), len(pending))
                client.notes["engine"] = {
                    "basins": [{"id": d.id, "best": round(d.best, 4), "generation": d.generation, "inflight": d.inflight,
                                "phase": d.phase, "batch": d.batch, "born_min": round((d.born - t0) / 60, 1)} for d in ld],
                    "division": list(division), "candidates": len(candidates), "pending": len(pending),
                    "caps": caps, "dead_axes": sorted(set(book.names) - set(book.live(clock())))}
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
