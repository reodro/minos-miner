"""The search engine: asynchronous basin climbing. Owner: search.

One streaming loop with no barriers. A round allows only ~15-20 dependent steps (a GATK traversal is minutes under
load) and as much width as the fleet has, so every step is as wide as the free workers allow and an improvement is
used the moment it arrives:

  opening     seeds (GATK default, the official seed, earlier rounds' bests) + a designed wave over the dominant
              axes (axes.design_wave), every skeleton with a genotyping family as riders
  basins      k regions (MIN_FAR live axes apart, a valley between them confirmed by one midpoint run; k from the
              fast nodes' capacity), each climbed asynchronously; born freely in the opening, later only by a leader
              (basins.py): a move is dispatched whenever this basin has a free worker, in the latency lane; an
              accept re-centres at once
  scouts      a share of the width keeps exploring (crossover and mutations of the basins' bests, unrun design
              cells, a few random configs), ranked by the (core, FP) forest; a scout that beats the weakest basin
              by a margin births a basin, if it is a region of its own
  finish      the last minutes: dense genotyping + PairHMM families on the best skeletons, half steps only

Every count scales with the executor's capacity, so a 16-worker replay and a 2,000-worker fleet run the same code.
Main submits the local best; the engine never submits.
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
from search.axes import AxisBook, design_wave, far_axes, rms_distance
from search.base import Strategy
from search.basins import Basin, Move
from search.client import EvalClient, Outcome
from search.families import genotyping_riders, pairhmm_riders, screen_moves, structural_riders
from search.ranker import Ranker
from search.space import Space, bam_caps

log = logging.getLogger(__name__)
TRACE = False                 # print every fill decision (debugging on the virtual clock)

RADIUS = 0.15                 # basins this far apart (RMS unit distance over structural axes) are distinct regions
MIN_FAR = 4                   # ... and now: on at least this many live axes far apart (> 10% of the range)
OPENING_S = 12 * 60           # basins are born freely only this long; later a newcomer must lead the best
LEAD_MARGIN = 0.05            # ... by this much (and sit on a confirmed hill of its own)
WORTH_A_HILL = 0.5            # a candidate this far below the best is not worth a valley test (the legacy margin)
FAST_FACTOR = 1.5             # nodes within this factor of the fastest are the climbing lane's capacity
BIRTH_UNTIL_S = 15 * 60       # no new basin in the last 15 min (it could not be climbed in time)
FINISH_S = 8 * 60             # the finish phase
RETIRE_BEHIND = 0.10          # a basin this far behind the leader after its patience is retired
RETIRE_PATIENCE_S = 12 * 60
BIRTH_MARGIN = 0.02           # a candidate must beat the weakest basin by this to replace it
REPLACE_AGE_S = 5 * 60        # ... and that basin must have had its chance: this old, or stalled for STALL_S
STALL_S = 4 * 60
MOVES_PER_BASIN = 200         # moves in flight per basin on a big fleet: sets how many basins the fleet races
MAX_BASINS = 12
SCOUT_SHARE, SCOUT_SHARE_LATE = 0.25, 0.15
RIDERS, RIDERS_SCOUT, RIDERS_SMALL = 8, 4, 6      # genotyping riders: the threshold sweep, the indel prior, contamination
STRUCTURAL_RIDERS = 4                             # a scout's family: structurally different configs in its traversal
SCREEN_S = FINISH_S + 10 * 60                     # the last stretch: multi-axis screens around the peak take the width
SMALL_FLEET = 64
CAPACITY_WAIT_S, CAPACITY_SETTLE_S, CAPACITY_POLL_S = 120.0, 20.0, 5.0   # a fleet's nodes join in the first minute

# A known strong config (official 88.163 on 2026-09-28 09:29): a seed, never a claim about the next round.
OFFICIAL_SEED = {
    "min_base_quality_score": 16, "min_mapping_quality_score": 46, "base_quality_score_threshold": 22,
    "standard_min_confidence_threshold_for_calling": 80.0, "pcr_indel_model": "NONE", "min_pruning": 2,
    "max_alternate_alleles": 6, "min_dangling_branch_length": 20, "recover_all_dangling_branches": False,
    "max_num_haplotypes_in_population": 128, "active_probability_threshold": 0.002, "min_assembly_region_size": 50,
    "max_assembly_region_size": 270, "assembly_region_padding": 100, "pair_hmm_gap_continuation_penalty": 7,
    "phred_scaled_global_read_mismapping_rate": 45, "heterozygosity": 0.002154, "indel_heterozygosity": 0.000125,
    "contamination_fraction_to_filter": 0.04, "max_reads_per_alignment_start": 116,
}
# The family shape that won our 2026-10-02 test round locally (its basin, not its exact values).
FAMILY_SEED = {
    "min_base_quality_score": 35, "min_mapping_quality_score": 8, "base_quality_score_threshold": 37,
    "standard_min_confidence_threshold_for_calling": 74.0, "pcr_indel_model": "AGGRESSIVE", "min_pruning": 2,
    "min_dangling_branch_length": 20, "recover_all_dangling_branches": False, "max_num_haplotypes_in_population": 17,
    "active_probability_threshold": 0.005, "min_assembly_region_size": 18, "max_assembly_region_size": 100,
    "assembly_region_padding": 13, "pair_hmm_gap_continuation_penalty": 13,
    "phred_scaled_global_read_mismapping_rate": 48, "heterozygosity": 0.008, "max_reads_per_alignment_start": 64,
}


@dataclass
class Unit:
    """One proposal: a head skeleton and its riders."""
    role: str                                   # seed | design | scout | move | finish
    head: Config
    basin: Optional[int] = None
    move: Optional[Move] = None
    best: Optional[float] = None
    best_config: Optional[Config] = None


async def settled_capacity(client: EvalClient, workers: int) -> int:
    """The workers the round has. A fleet's nodes join during the first minute (self-test, bundle push), so the
    first number seen would size a 1,000-worker fleet as a one-basin search: wait until the count has stood still
    for CAPACITY_SETTLE_S (at most CAPACITY_WAIT_S). Local and virtual executors answer at once."""
    cap = client.capacity()
    nodes = getattr(client.executor, "nodes", None)
    if nodes is None:                                          # not a fleet: its workers, now
        return cap or workers or 8
    t0 = time.time()
    stable_since = t0
    while time.time() - t0 < CAPACITY_WAIT_S:
        await asyncio.sleep(CAPACITY_POLL_S)
        c = client.capacity() or 0
        if c != cap:
            cap, stable_since = c, time.time()
            continue
        joining = [n for n in nodes() if not n["ready"] and not n["rejected"]]   # still in its self-test or push
        if cap and not joining and time.time() - stable_since >= CAPACITY_SETTLE_S:
            break
    log.info("engine: fleet capacity settled at %d workers after %.0f s (%d nodes ready)", cap or 0,
             time.time() - t0, sum(1 for n in nodes() if n["ready"]))
    return cap or workers or 8


def searcher(ranges: ParamRanges, workers: int, *, tool: str = "gatk", seed: int = 0, tag: str = "engine") -> Strategy:
    rng = random.Random(seed)

    async def run(client: EvalClient) -> None:
        r = client.round
        caps = bam_caps(r.bam, r.chrom, r.start, r.end) if r.bam else {}
        space = Space(ranges, tool, caps=caps)                 # the round's space: its data bounds the axes
        if caps:
            log.info("engine: caps from the round's BAM: %s", caps)
        book = AxisBook(space)
        ranker = Ranker(space, seed=seed)
        clock = client.clock
        t0 = clock()
        cap = await settled_capacity(client, workers)
        small = cap < SMALL_FLEET
        riders = RIDERS_SMALL if small else RIDERS
        riders_scout = max(2, RIDERS_SMALL // 2) if small else RIDERS_SCOUT

        def fast_capacity() -> int:
            """Workers on the nodes within FAST_FACTOR of the fastest: a sweep is as slow as its slowest move, so
            the climbing lane is sized to them (a lone machine, or a virtual executor: everything)."""
            nodes = getattr(client.executor, "nodes", None)
            rows = [n for n in (nodes() if nodes else []) if n.get("ready") and n.get("workers")]
            if not rows:
                return cap
            fastest = min(n["solo_s"] for n in rows)
            return sum(n["workers"] for n in rows if n["solo_s"] <= FAST_FACTOR * fastest) or cap

        fast_cap = fast_capacity()
        k_max = 1 if small else max(1, min(MAX_BASINS, round(fast_cap * (1 - SCOUT_SHARE) / MOVES_PER_BASIN)))
        virtual = hasattr(client.executor, "now")
        pending: Dict[asyncio.Future, Tuple[Unit, bool]] = {}      # future -> (unit, is_head)
        proposed: Set[str] = set()
        basins: Dict[int, Basin] = {}
        next_basin = [0]
        candidates: Dict[str, Tuple[Config, float]] = {}            # basin candidates: key -> (config, family best)
        retired_centres: List[Tuple[Config, float]] = []           # where retired basins ended, and their bests
        was_centre: Set[str] = set()                                # candidate keys that already seeded a basin
        valley: Dict[str, Tuple[Config, float, int]] = {}          # midpoint key -> (candidate, its score, basin id)
        valley_done: Set[str] = set()                               # candidate keys already tested
        inflight_scouts = [0]
        scout_queue: List[Config] = []                              # ranked scouts waiting for a slot
        design_left: List[Config] = []
        fitting: List[Optional[asyncio.Task]] = [None]
        finish_sent = [False]
        last_log = [t0]
        wake = asyncio.Event()
        client.on_result.append(lambda _o: wake.set())

        def elapsed() -> float:
            return clock() - t0

        def left() -> Optional[float]:
            return client.time_left()

        def fresh(c: Optional[Config]) -> bool:
            return c is not None and c.key not in proposed and (c.key, FULL) not in client.latest

        def live() -> List[Basin]:
            return [b for b in basins.values() if not b.retired]

        def leader() -> Optional[Basin]:
            lb = live()
            return max(lb, key=lambda b: b.best) if lb else None

        def best_known() -> float:
            return max((s for _, s in candidates.values()), default=0.0)

        def scout_target() -> int:
            t = left()
            share = SCOUT_SHARE_LATE if (t is not None and t < FINISH_S + 10 * 60) else SCOUT_SHARE
            return max(1, int(cap * share)) if live() else cap

        # ---- proposing ---------------------------------------------------------------------------------------
        def propose(unit: Unit, rider_configs: Sequence[Config], priority: int, lane: str) -> int:
            batch = [unit.head] + [c for c in rider_configs if fresh(c)]
            batch = [c for c in batch if fresh(c) or c is unit.head]
            if not batch:
                return 0
            lanes = [lane] + ["throughput"] * (len(batch) - 1)            # riders pack into the head's batch
            futs = client.propose_many(batch, priority=priority, tag=f"{tag}:{unit.role}", lanes=lanes)
            for c, f in zip(batch, futs):
                proposed.add(c.key)
                if f.done():                                   # cached by the client: handle it on the next pass
                    pending[f] = (unit, c is unit.head)
                else:
                    pending[f] = (unit, c is unit.head)
            return len(futs)

        def propose_skeleton(role: str, head: Config, n_riders: int, priority: int, lane: str,
                             basin: Optional[Basin] = None, move: Optional[Move] = None) -> bool:
            unit = Unit(role, head, basin.id if basin else None, move)
            kids = genotyping_riders(space, head, n_riders, fresh)
            if role in ("seed", "design", "scout"):              # width: a family of structural variants per traversal
                kids += structural_riders(space, head, STRUCTURAL_RIDERS, rng, fresh)
            sent = propose(unit, kids, priority, lane)
            if TRACE:
                print(f"   propose {role}: {sent} futures, {len(kids)} riders, head fresh={fresh(head)} idle_after={client.idle_slots()}")
            if sent == 0:
                return False
            if basin is not None:
                pass                                         # basin.inflight was counted by take()
            elif role != "finish":
                inflight_scouts[0] += 1
            return True

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
            n_design = max(0, cap - len(seeds))
            cells = design_wave(space, rng, n_design, seeds[:3])
            for c in seeds:
                propose_skeleton("seed", c, riders, 5, "throughput")
            design_left.extend(cells)
            log.info("engine: capacity %d (%d on the fast nodes) -> %d basins max, %d riders; %d seeds, %d design cells",
                     cap, fast_cap, k_max, riders, len(seeds), len(cells))

        def scout_candidates(n: int) -> List[Config]:
            out: List[Config] = []
            while design_left and len(out) < n:
                c = design_left.pop(0)
                if fresh(c):
                    out.append(c)
            if len(out) >= n:
                return out
            bests = [b.centre for b in live()]
            if not bests:
                bests = [c for c, _ in sorted(candidates.values(), key=lambda t: -t[1])[:4]]
            pool: List[Config] = []
            names = book.names
            t_left = left()
            if bests and t_left is not None and t_left < SCREEN_S:   # the last stretch: screens around the peaks
                ld = leader()
                tops = [ld.centre] if ld is not None else []
                tops += [c for c, _ in sorted(candidates.values(), key=lambda t: -t[1])[:3] if all(c.key != x.key for x in tops)]
                for top in tops:
                    pool += screen_moves(space, top, rng, max(1, n // max(1, len(tops))), book.live(clock()), fresh)
            tries = 0
            while len(pool) < 4 * n and tries < 20 * n:
                tries += 1
                r = rng.random()
                if r < 0.10 or not bests:
                    c = space.random(rng)
                elif r < 0.45 and len(bests) >= 2:
                    a, b = rng.sample(bests, 2)
                    params = dict(a.params)
                    for name in names:
                        if rng.random() < 0.5:
                            params[name] = b.params.get(name, space.specs[name].default)
                    c = space._assemble(params)
                    if c is not None and rng.random() < 0.3:
                        c = space.neighbour(c, rng, n_params=2, scale=0.2, pool=names) or c
                elif r < 0.85:
                    c = space.neighbour(rng.choice(bests), rng, n_params=rng.randint(1, 3), scale=0.2, pool=names)
                else:
                    base = rng.choice(bests)
                    name = rng.choice(names)
                    spec = space.specs[name]
                    if spec.allowed or spec.type == "bool":
                        continue
                    params = dict(base.params)
                    params[name] = space._from_unit(spec, 0.0 if rng.random() < 0.5 else 1.0)
                    c = space._assemble(params)
                if c is None or not fresh(c) or any(c.key == p.key for p in pool):
                    continue
                if not all(book.allowed(a, c.params.get(a)) for a in names):
                    continue
                pool.append(c)
            if ranker.ready and pool:
                lb = leader()
                order = ranker.rank(pool, best=lb.best if lb else 0.0)
                pool = [pool[i] for i in order]
            out.extend(pool[:n - len(out)])
            return out

        def finish_round() -> None:
            """Dense families on the best skeletons of every surviving basin, and the top candidates."""
            finish_sent[0] = True
            tops: List[Config] = [b.centre for b in sorted(live(), key=lambda b: -b.best)]
            for c, _ in sorted(candidates.values(), key=lambda t: -t[1])[:max(1, cap // 40)]:
                if all(rms_distance(space, book.names, c, t) > 1e-9 for t in tops):
                    tops.append(c)
            n = 0
            for c in tops[:max(1, cap // 40)]:
                kids = genotyping_riders(space, c, riders, fresh) + pairhmm_riders(space, c, 6, fresh)
                unit = Unit("finish", c)
                n += propose(unit, kids, 2, "throughput")
            for a in book.axes.values():                    # the last steps are half steps
                a.step = max(0.02, a.step * 0.5)
            log.info("engine: finish at t+%.1fm: %d rider configs on %d skeletons", elapsed() / 60, n, len(tops))

        def fill() -> int:
            n = 0
            while not client.exhausted:
                idle = client.idle_slots()
                if idle is not None and idle <= 0:
                    break
                if idle is None and len(pending) >= cap:
                    break
                t = left()
                if t is not None and t < FINISH_S:
                    if not finish_sent[0]:
                        finish_round()
                        n += 1
                        continue
                lb = live()
                want_scouts = scout_target()
                picked = None
                if TRACE:
                    print(f"fill t+{elapsed()/60:.1f}m idle={idle} scouts={inflight_scouts[0]}/{want_scouts} basins="
                          + " ".join(f"{b.id}:pool{len(b.pool)}/inflight{b.inflight}/gen{b.generation}/refilled{b.refilled_generation}" for b in lb))
                if lb and inflight_scouts[0] >= want_scouts:
                    ready = [b for b in lb if b.pool or (b.wants_refill()
                                                         and b.refill(rng, [o.centre for o in lb if o is not b], clock(), ranker))]
                    if ready:
                        target = max(1.0, fast_cap * (1 - SCOUT_SHARE) / max(1, len(lb)))   # a sweep per basin, on fast nodes
                        under = [b for b in ready if b.inflight < target]
                        picked = min(under, key=lambda b: b.inflight / target) if under else None
                if picked is not None:
                    move = picked.take()
                    if move is not None and propose_skeleton("move", move.config, riders, 3, "latency", picked, move):
                        n += 1
                        continue
                    picked.inflight = max(0, picked.inflight - 1)
                if t is None or t >= FINISH_S:
                    if not scout_queue:
                        scout_queue.extend(scout_candidates(max(1, min(64, want_scouts - inflight_scouts[0], idle or 64))))
                    while scout_queue:
                        c = scout_queue.pop(0)
                        if fresh(c):
                            break
                    else:
                        c = None
                    if c is not None and propose_skeleton("scout", c, riders_scout, 1, "throughput"):
                        n += 1
                        continue
                if picked is None and lb:                    # scouts full (or off) and no basin picked: climb anyway
                    ready = [b for b in lb if b.pool or (b.wants_refill()
                                                         and b.refill(rng, [o.centre for o in lb if o is not b], clock(), ranker))]
                    if ready:
                        b = min(ready, key=lambda b: b.inflight)
                        move = b.take()
                        if move is not None and propose_skeleton("move", move.config, riders, 3, "latency", b, move):
                            n += 1
                            continue
                        b.inflight = max(0, b.inflight - 1)
                # nothing left to climb and the scout share is full, yet workers are idle: idle capacity never
                # waits - width is cheap, and a basin's sweep in flight cannot use the slot anyway
                if not scout_queue:
                    scout_queue.extend(scout_candidates(max(1, min(32, idle or 32))))
                c = next((c for c in scout_queue if fresh(c)), None)
                if c is not None:
                    scout_queue.remove(c)
                    if propose_skeleton("scout", c, riders_scout, 1, "throughput"):
                        n += 1
                        continue
                break
            return n

        # ---- results -----------------------------------------------------------------------------------------
        def note_result(unit: Unit, is_head: bool, o: Outcome) -> None:
            now = clock()
            score = o.score if o.ok else None
            failed = (not o.ok) or score is None or score <= 0
            r = o.result
            ranker.add(o.config, r.core if (r and not failed) else None,
                       (r.region_fp / r.n_targets) if (r and not failed and r.n_targets) else None, failed)
            if not failed and (unit.best is None or score > unit.best):
                unit.best, unit.best_config = score, o.config
            b = basins.get(unit.basin) if unit.basin is not None else None
            if b is not None:
                if is_head:
                    b.inflight = max(0, b.inflight - 1)
                if not b.retired:
                    moved = b.consider(unit.move if is_head else None, o.config, score, failed, now)
                    if moved:
                        log.info("engine: basin %d -> %.4f via %s %s at t+%.1fm", b.id, b.best,
                                 unit.move.kind if unit.move else "rider", ",".join(unit.move.axes) if unit.move else "",
                                 elapsed() / 60)
                    return
            if is_head and unit.role in ("seed", "design", "scout", "valley"):
                inflight_scouts[0] = max(0, inflight_scouts[0] - 1)
            if is_head and unit.role == "valley":
                judge_valley(unit.head.key, score if not failed else None, now)
            if not failed:
                k = unit.head.key
                if k not in candidates or score > candidates[k][1]:
                    candidates[k] = (o.config, score)
                if len(candidates) > 600:                       # keep the pool small: the top 300 by score
                    for key in sorted(candidates, key=lambda x: -candidates[x][1])[300:]:
                        del candidates[key]

        def live_names() -> List[str]:
            return book.live(clock())

        def same_region(a: Config, b: Config) -> bool:
            return far_axes(space, live_names(), a, b) < MIN_FAR

        def test_valley(config: Config, score: float, basin: Basin, now: float) -> bool:
            """One run at the midpoint between a candidate and the nearest basin's centre (on the live axes): a
            midpoint as good as both ends means one hill - the candidate becomes a transplant source there; a dip
            means a valley - the candidate founds a basin. Returns True when the test was sent."""
            names = live_names()
            params = dict(config.params)
            for name in names:
                s = space.specs[name]
                if s.allowed or s.type == "bool":
                    continue
                u = (space._to_unit(s, config.params.get(name, s.default)) + space._to_unit(s, basin.centre.params.get(name, s.default))) / 2
                params[name] = space._from_unit(s, u)
            mid = space._assemble(params)
            valley_done.add(config.key)
            if mid is None or not fresh(mid):
                return False
            valley[mid.key] = (config, score, basin.id)
            return propose_skeleton("valley", mid, riders_scout, 2, "throughput")

        def judge_valley(mid_key: str, mid_score: Optional[float], now: float) -> None:
            item = valley.pop(mid_key, None)
            if item is None:
                return
            config, score, basin_id = item
            b = basins.get(basin_id)
            if b is None or b.retired:
                return
            if mid_score is not None and mid_score >= min(score, b.best) - 0.005:   # no dip: the same hill
                b.sources.append(config)
                top_cfg, top = (config, score) if score >= mid_score else (candidates.get(mid_key, (None, 0.0)))
                if top_cfg is not None and b.adopt(top_cfg, top, now):          # the higher point on the hill is its centre
                    was_centre.add(top_cfg.key)
                    log.info("engine: basin %d -> %.4f, the %s on its hill (midpoint %.4f) t+%.1fm", b.id, top,
                             "candidate" if top_cfg is config else "midpoint", mid_score, elapsed() / 60)
                else:
                    log.info("engine: candidate %.4f sits on basin %d's hill (midpoint %.4f) t+%.1fm", score, b.id,
                             mid_score, elapsed() / 60)
                return
            if len(live()) < k_max or score > max(x.best for x in live()) + LEAD_MARGIN:
                log.info("engine: a valley between %.4f and basin %d (midpoint %s) t+%.1fm", score, b.id,
                         f"{mid_score:.4f}" if mid_score is not None else "failed", elapsed() / 60)
                found(config, score, now)

        def found(config: Config, score: float, now: float) -> None:
            lb = live()
            if len(lb) >= k_max:                                  # no slot: the weakest basin that had its chance
                settled = [b for b in lb if now - b.born > REPLACE_AGE_S or now - b.improved_at > STALL_S]
                if not settled:
                    return
                weakest = min(settled, key=lambda b: b.best)
                last_gain = (weakest.history[-1][1] - weakest.history[-2][1]) if len(weakest.history) >= 2 else 0.0
                if score - weakest.best < max(LEAD_MARGIN, last_gain):
                    return
                weakest.retired = True
                retired_centres.append((weakest.centre, weakest.best))
                log.info("engine: basin %d (%.4f) replaced by a better region at t+%.1fm", weakest.id, weakest.best,
                         elapsed() / 60)
            birth(config, score, now)

        def births_and_retirements() -> None:
            now = clock()
            t = left()
            lb = live()
            # merge basins that converged (fewer than MIN_FAR live axes apart: one hill)
            for a, b in [(a, b) for a in lb for b in lb if a.id < b.id]:
                if not a.retired and not b.retired and same_region(a.centre, b.centre):
                    worse = a if a.best < b.best else b
                    worse.retired = True
                    retired_centres.append((worse.centre, worse.best))
                    log.info("engine: basin %d merged into its neighbour at t+%.1fm", worse.id, elapsed() / 60)
            lb = live()
            ld = leader()
            for b in lb:
                if b.exhausted() or (ld is not None and b is not ld and now - b.improved_at > RETIRE_PATIENCE_S
                                     and b.best < ld.best - RETIRE_BEHIND):
                    b.retired = True
                    retired_centres.append((b.centre, b.best))
                    log.info("engine: basin %d retired at %.4f (%s) t+%.1fm", b.id, b.best,
                             "exhausted" if b.exhausted() else "behind", elapsed() / 60)
            if t is not None and t < BIRTH_UNTIL_S:
                return
            lb = live()
            ld = leader()
            opening = elapsed() < OPENING_S
            ranked = sorted(candidates.values(), key=lambda c: -c[1])
            tests = 0
            for config, score in ranked[:50]:
                lb = live()
                if config.key in was_centre:
                    continue
                here = [b for b in lb if same_region(b.centre, config)]
                if here:                                         # on a known hill: higher than its centre, it IS the centre
                    b = max(here, key=lambda x: x.best)
                    if b.adopt(config, score, now):
                        was_centre.add(config.key)
                        log.info("engine: basin %d -> %.4f, a scout on its hill, at t+%.1fm", b.id, score, elapsed() / 60)
                    continue
                if config.key in valley_done:                    # tested already: on some hill, or waiting for its midpoint
                    continue
                if ld is not None and score < max(ld.best, best_known()) - WORTH_A_HILL:
                    break                                        # the rest is too low to be worth a run each
                if any(same_region(c, config) and score <= best + LEAD_MARGIN for c, best in retired_centres):
                    continue                                    # an explored region, and no better than it ended
                if not lb:                                       # the first basin needs no valley
                    birth(config, score, now)
                    continue
                if not opening and ld is not None and score <= ld.best + LEAD_MARGIN:
                    break                                        # after the opening only a clear leader founds a basin
                if len(lb) >= k_max and not opening:
                    break
                nearest = min(lb, key=lambda b: b.far(config, live_names()))
                if test_valley(config, score, nearest, now):     # the hill question is answered by one run
                    tests += 1
                if tests >= 4:
                    break

        def birth(config: Config, score: float, now: float) -> None:
            b = Basin(next_basin[0], space, book, config, score, now, now)
            next_basin[0] += 1
            basins[b.id] = b
            was_centre.add(config.key)
            b.refill(rng, [o.centre for o in live() if o is not b], now, ranker)
            log.info("engine: basin %d born at %.4f (%d moves) t+%.1fm", b.id, score, len(b.pool), elapsed() / 60)

        def maybe_fit() -> None:
            if ranker.stale() and (fitting[0] is None or fitting[0].done()):
                fitting[0] = asyncio.create_task(asyncio.to_thread(ranker.fit))

        # ---- the loop ----------------------------------------------------------------------------------------
        opening()
        fill()
        while True:
            if pending:
                done, _ = await asyncio.wait(list(pending), return_when=asyncio.FIRST_COMPLETED)
                for fut in done:
                    unit, is_head = pending.pop(fut)
                    try:
                        note_result(unit, is_head, fut.result())
                    except Exception:  # noqa: BLE001 - one bad result must not stop the round
                        log.exception("engine: result handling failed")
                births_and_retirements()
                maybe_fit()
            n = fill()
            if clock() - last_log[0] > 120:
                last_log[0] = clock()
                lb = live()
                log.info("engine: t+%.1fm basins %s | scouts in flight %d | candidates %d | pending %d",
                         elapsed() / 60, " ".join(f"{b.id}:{b.best:.3f}/{b.inflight}" for b in lb),
                         inflight_scouts[0], len(candidates), len(pending))
                client.notes["engine"] = {                      # what the dashboard shows of the search
                    "basins": [{"id": b.id, "best": round(b.best, 4), "generation": b.generation, "inflight": b.inflight,
                                "born_min": round((b.born - t0) / 60, 1)} for b in lb],
                    "scouts_in_flight": inflight_scouts[0], "candidates": len(candidates), "pending": len(pending),
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
