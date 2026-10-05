"""The search engine (v4): a global genetic opening, ParamILS climbs on up to six hills, merges by recombination. Owner: search.

  opening     wave 0 = seeds + a Latin hypercube over every setting (4 threshold riders each), sent while the answer
              set is still being built; then GENERATIONS genetic generations with niching (genetic.py), each launched
              once GEN_BACK of the previous one is back or workers sit idle (no barrier). Stops early when a
              generation gains less than GEN_MIN_GAIN
  founding    the hill registry (hills.py) groups the good configs into hills; one domain per hill, up to N_DOMAINS
  climb       every domain runs iterated local search (climber.py); winning settings are shared between domains
              every SHARE_EVERY_S; the scout keeps SCOUT_SHARE of the workers on the whole space for the whole round
  merges      one operation everywhere: pairing (every combination of the settings where two centres differ, when
              few) plus a genetic generation of the two (genetic.merge_children); the best of parents and children is
              the merged domain's centre, which inherits both parents' tested configs and steps. Merges happen
                - on evidence, at any time: two domains on one hill, or a frozen domain on a live one's hill
                - on schedule: closest pairs first from REDUCE_TO_3[0] so at most 3 domains remain at REDUCE_TO_3[1];
                  all remaining domains at REDUCE_TO_1[0] (a three-way merge) so one remains for the final
  final       the last FINAL_S on one point: halving single steps and all pairs, one last pairing with the runner-up,
              and the plateau-centre pick: among configs within PICK_BAND of the best, the one whose neighbours score
              highest (client.pick_key; main submits it)

One controller pass after every batch of results (and at least every TICK_S) decides everything from facts: what is
in flight, each domain's progress rate (own gain per minute over the last 10 minutes), the registry. Schedule times
are fractions of the search's time, set for a 72-minute round.
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
from search.families import genotyping_riders
from search.ranker import Ranker
from search.space import GENOTYPING_PARAMS, Space, bam_caps
from search import genetic
from search.climber import ACCEPT_MIN, Climber, Move
from search.hills import MIN_FAR, Hills

log = logging.getLogger(__name__)

ROUND_S = 72 * 60                    # the schedule below is written for this much search time; it scales with it
GENERATIONS = 3
GEN_BACK = 0.6                       # a generation launches once this share of the previous wave is back
GEN_IDLE = 0.25                      # ... or this share of the workers is idle with half of it back
GEN_MIN_GAIN = 0.05                  # a generation gaining less than this ends the opening early
WAVE_BACK = 0.8
OPENING_MAX_S = 14 * 60
N_DOMAINS = 6
WORTH_A_HILL = 0.5                   # no domain is founded this far below the best
REDUCE_TO_3 = (34 * 60, 40 * 60)     # scheduled merges start / at most 3 domains
REDUCE_TO_1 = (56 * 60, 62 * 60)     # the three-way merge starts / one domain
FINAL_S = 10 * 60
MERGE_BACK = 0.8
MERGE_CLOSE, MIDPOINT_RETEST_S = 0.02, 10 * 60
SHARE_EVERY_S, SHARE_AXES = 5 * 60, 3
SCOUT_SHARE, SCOUT_KAPPA = 0.15, 1.0
RIDERS_WAVE0, RIDERS_GLOBAL = 4, 2
PICK_BAND, PICK_NEAR, PICK_MIN, PICK_S = 0.05, 2, 3, 90.0
TICK_S = 2.0
SMALL_FLEET = 64

# Known strong configs (a seed is a start, never a claim about the round).
OFFICIAL_SEED = {
    "min_base_quality_score": 16, "min_mapping_quality_score": 46, "base_quality_score_threshold": 22,
    "standard_min_confidence_threshold_for_calling": 80.0, "pcr_indel_model": "NONE", "min_pruning": 2,
    "max_alternate_alleles": 6, "min_dangling_branch_length": 20, "recover_all_dangling_branches": False,
    "max_num_haplotypes_in_population": 128, "active_probability_threshold": 0.002, "min_assembly_region_size": 50,
    "max_assembly_region_size": 270, "assembly_region_padding": 100, "pair_hmm_gap_continuation_penalty": 7,
    "phred_scaled_global_read_mismapping_rate": 45, "heterozygosity": 0.002154, "indel_heterozygosity": 0.000125,
    "contamination_fraction_to_filter": 0.04, "max_reads_per_alignment_start": 116,
}
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
    role: str                        # wave | scout | move | merge | midpoint
    head: Config
    domain: Optional[int] = None
    move: Optional[Move] = None
    merge: Optional[int] = None
    wave: int = 0


@dataclass
class MergeJob:
    id: int
    parents: List[int]
    why: str
    sent: int = 0
    back: int = 0
    best: Optional[Tuple[Config, float]] = None
    started: float = 0.0


async def settled_capacity(client: EvalClient, workers: int, wait_s: float = 120.0, settle_s: float = 20.0) -> int:
    """The workers the round has, once a fleet's nodes have joined (local and virtual executors answer at once)."""
    cap = client.capacity()
    nodes = getattr(client.executor, "nodes", None)
    if nodes is None:
        return cap or workers or 8
    t0 = stable = time.time()
    while time.time() - t0 < wait_s:
        await asyncio.sleep(5.0)
        c = client.capacity() or 0
        if c != cap:
            cap, stable = c, time.time()
            continue
        joining = [n for n in nodes() if not n["ready"] and not n["rejected"]]
        if cap and not joining and time.time() - stable >= settle_s:
            break
    log.info("engine: fleet capacity settled at %d workers after %.0f s", cap or 0, time.time() - t0)
    return cap or workers or 8


def searcher(ranges: ParamRanges, workers: int, *, tool: str = "gatk", seed: int = 0, tag: str = "engine") -> Strategy:
    rng = random.Random(seed)

    async def run(client: EvalClient) -> None:
        r = client.round
        caps = bam_caps(r.bam, r.chrom, r.start, r.end) if r.bam else {}
        space = Space(ranges, tool, caps=caps)
        book = AxisBook(space)
        ranker = Ranker(space, seed=seed)
        clock = client.clock
        t0 = clock()
        cap = await settled_capacity(client, workers)
        tl = client.time_left()
        scale_t = (tl if tl is not None and tl > 0 else ROUND_S) / ROUND_S     # the schedule follows the round's length
        T = lambda s: s * scale_t
        n_domains = max(1, min(N_DOMAINS, cap // 8))
        virtual = hasattr(client.executor, "now")
        genotyping = [a for a in GENOTYPING_PARAMS if a in space.names]

        pending: Dict[asyncio.Future, Tuple[Unit, bool]] = {}
        proposed: Set[str] = set()
        domains: Dict[int, Climber] = {}
        status: Dict[int, str] = {}                  # live | merging | retired
        merges: Dict[int, MergeJob] = {}
        ids = [0, 0]                                 # next domain id, next merge id
        candidates: Dict[str, Tuple[Config, float]] = {}
        was_centre: Set[str] = set()
        wave_heads: List[int] = []                   # heads sent per wave (0 = hypercube, 1.. = generations)
        wave_back: List[int] = []
        wave_best: List[float] = []                  # best score seen when each wave launched
        st = {"founded": False, "gens_done": False, "r3": False, "r3_done": False, "r1": False, "final": False,
              "last_share": 0.0, "last_log": t0, "scouts": 0, "runner_up": None, "picked": False}
        midpoint: Dict[str, Tuple[int, int]] = {}
        find_test: Dict[str, Tuple[Config, float, int]] = {}     # midpoint key -> (scout find, its score, nearest domain)
        find_done: Set[str] = set()                              # finds already tested
        tested_pair: Dict[Tuple[int, int], float] = {}
        scout_queue: List[Config] = []
        fitting: List[Optional[asyncio.Task]] = [None]
        wake = asyncio.Event()
        client.on_result.append(lambda _o: wake.set())

        def elapsed() -> float:
            return clock() - t0

        def left() -> Optional[float]:
            return client.time_left()

        def fresh(c: Optional[Config]) -> bool:
            return c is not None and c.key not in proposed and (c.key, FULL) not in client.latest

        def names() -> List[str]:
            return book.live(clock()) + genotyping

        hills = Hills(space, names)

        def live() -> List[Climber]:
            return [d for i, d in domains.items() if status[i] == "live" and d.state != "exhausted"]

        def standing() -> List[Climber]:
            return [d for i, d in domains.items() if status[i] != "retired"]

        def centres(which: Sequence[Climber]) -> Dict[int, Config]:
            return {d.id: d.centre for d in which}

        def owner_of(c: Config) -> Optional[int]:
            return hills.owner(c, centres([d for i, d in domains.items() if status[i] == "live"]))

        def top() -> float:
            return max((s for _, s in candidates.values()), default=0.0)

        # ---- proposing ----------------------------------------------------------------------------------------------
        def propose(unit: Unit, kids: Sequence[Config], priority: int) -> int:
            batch = [unit.head] + [c for c in kids if fresh(c)]
            if not fresh(unit.head) and unit.role != "move":
                batch = batch[1:]
            if not batch:
                return 0
            kind = unit.move.kind if unit.move is not None else unit.role
            # a result the search waits on (a move, a merge child, a generation child, a midpoint or valley test) is a
            # small batch: one config per worker, in parallel, carrying only its cheap siblings. Wave 0 and the scout
            # pack freely (big): nothing waits on any single one of them.
            small = unit.role in ("move", "merge", "midpoint", "findtest") or (unit.role == "wave" and unit.wave > 0)
            batches = [("small" if small else "big")] + ["big"] * (len(batch) - 1)
            if batch[0] is not unit.head:
                batches = ["big"] * len(batch)
            futs = client.propose_many(batch, priority=priority, tag=f"{tag}:{kind}", batches=batches)
            for c, f in zip(batch, futs):
                proposed.add(c.key)
                pending[f] = (unit, c is unit.head)
            return len(futs)

        def send_global(c: Config, role: str, riders: int, priority: int, wave: int = 0) -> bool:
            return propose(Unit(role, c, wave=wave), genotyping_riders(space, c, riders, fresh), priority) > 0

        # ---- opening --------------------------------------------------------------------------------------------
        def seeds() -> List[Config]:
            out: List[Config] = []
            for p in ({}, OFFICIAL_SEED, FAMILY_SEED):
                c = space._assemble(dict(p))
                if c is not None and all(c.key != x.key for x in out):
                    out.append(c)
            for raw in client.store.portfolio(limit=max(4, cap // 20), tool=tool):
                c = space._assemble(dict(raw.params))
                if c is not None and all(c.key != x.key for x in out):
                    out.append(c)
            return out

        def launch_wave(configs: Sequence[Config], riders: int, label: str) -> None:
            w = len(wave_heads)
            wave_heads.append(0); wave_back.append(0); wave_best.append(top())
            for c in configs:
                if send_global(c, "wave", riders, 5 if w == 0 else 4, wave=w):
                    wave_heads[w] += 1
            log.info("engine: wave %d (%s) at t+%.1fm: %d points", w, label, elapsed() / 60, wave_heads[w])

        def opening_step() -> None:
            if st["founded"]:
                return
            w = len(wave_heads) - 1
            back = wave_back[w] / max(1, wave_heads[w])
            idle = client.idle_slots()
            if elapsed() >= T(OPENING_MAX_S):                  # the opening's limit is hard: found now
                st["gens_done"] = True
                found()
                return
            if not st["gens_done"]:
                can = back >= GEN_BACK or (back >= 0.5 and idle is not None and idle >= GEN_IDLE * cap)
                if can and len(candidates) >= 4:
                    if w >= 2 and top() - wave_best[w] < GEN_MIN_GAIN:        # the last generation gained too little
                        st["gens_done"] = True
                        log.info("engine: generations stop after %d (last gained %.3f) t+%.1fm", w, top() - wave_best[w],
                                 elapsed() / 60)
                    else:
                        children = [c for c in genetic.generation(space, list(candidates.values()), rng, cap, hills, names())
                                    if fresh(c)]
                        launch_wave(children, RIDERS_GLOBAL, f"generation {w + 1}")
                        st["gens_done"] = len(wave_heads) - 1 >= GENERATIONS
                return
            if back >= WAVE_BACK or not pending:
                found()

        def found() -> None:
            st["founded"] = True
            ranked = sorted(candidates.values(), key=lambda t: -t[1])
            for c, s in hills.cluster(ranked, n_domains, WORTH_A_HILL):
                birth(c, s, "the opening")
            log.info("engine: founding at t+%.1fm: %d domains from %d scored configs (waves %s)", elapsed() / 60,
                     len(live()), len(candidates), " ".join(f"{b}/{h}" for b, h in zip(wave_back, wave_heads)))

        def birth(c: Config, s: float, why: str, parents: Sequence[Climber] = ()) -> Climber:
            d = Climber(ids[0], space, book, c, s, clock(), owner_of, extra_axes=genotyping)
            ids[0] += 1
            for p in parents:
                d.inherit(p.tried, p.steps, p.axis_gain)
                if p.mode == "fine":
                    d.mode = "fine"
            domains[d.id], status[d.id] = d, "live"
            was_centre.add(c.key)
            log.info("engine: domain %d born at %.4f (%s) t+%.1fm", d.id, s, why, elapsed() / 60)
            return d

        # ---- merges -----------------------------------------------------------------------------------------------
        def start_merge(parents: List[Climber], why: str) -> None:
            parents = sorted(parents, key=lambda d: -d.best)
            for p in parents:
                status[p.id] = "merging"
                p.pool = []
            mj = MergeJob(ids[1], [p.id for p in parents], why, started=clock(),
                          best=(parents[0].centre, parents[0].best))
            ids[1] += 1
            merges[mj.id] = mj
            if len(parents) == 2:
                kids = genetic.merge_children(space, parents[0].centre, parents[1].centre, rng, names())
            else:
                kids = genetic.merge3_children(space, [p.centre for p in parents], rng, names())
            tried = set().union(*(p.tried for p in parents))
            kids = [c for c in kids if c.key not in tried and fresh(c)]
            limit = max(24, cap // 4)                       # one wave's worth on a big fleet, a few on a small one
            if len(kids) > limit:
                kids = rng.sample(kids, limit)
            for c in kids:
                if c.key in tried or not fresh(c):
                    continue
                if propose(Unit("merge", c, merge=mj.id), [], 6):
                    mj.sent += 1
            log.info("engine: merge %d (%s): domains %s at %s, %d children, t+%.1fm", mj.id, why,
                     ",".join(str(p.id) for p in parents), " ".join(f"{p.best:.4f}" for p in parents), mj.sent,
                     elapsed() / 60)
            if mj.sent == 0:
                resolve(mj)

        def resolve(mj: MergeJob) -> Climber:
            merges.pop(mj.id, None)
            parents = [domains[i] for i in mj.parents]
            pool = [(p.centre, p.best, f"domain {p.id}") for p in parents]
            if mj.best is not None:
                pool.append((mj.best[0], mj.best[1], "a child"))
            c, s, src = max(pool, key=lambda t: t[1])
            for p in parents:
                status[p.id] = "retired"
            rest = sorted((p for p in parents if p.centre.key != c.key), key=lambda p: -p.best)
            if rest:
                st["runner_up"] = rest[0].centre
            d = birth(c, s, f"merge {mj.id}, {mj.why}: best from {src}, {mj.back}/{mj.sent} children back", parents)
            return d

        def same_hill_pairs(which: Sequence[Climber]) -> List[Tuple[int, Climber, Climber]]:
            out = []
            for i, a in enumerate(which):
                for b in which[i + 1:]:
                    out.append((hills.distance(a.centre, b.centre), a, b))
            return sorted(out, key=lambda t: t[0])

        def test_midpoint(a: Climber, b: Climber) -> None:
            p = dict(a.centre.params)
            for n in names():
                s = space.specs[n]
                if s.allowed or s.type == "bool":
                    continue
                u = (space._to_unit(s, a.centre.params.get(n, s.default)) + space._to_unit(s, b.centre.params.get(n, s.default))) / 2
                p[n] = space._from_unit(s, u)
            mid = space._assemble(p)
            tested_pair[(a.id, b.id)] = clock()
            if mid is not None and fresh(mid) and propose(Unit("midpoint", mid), [], 6):
                midpoint[mid.key] = (a.id, b.id)

        def test_find(c: Config, s: float) -> None:
            """One run halfway between a scout find and the nearest domain: a dip means a hill of its own."""
            find_done.add(c.key)
            sd = [d for d in standing()]
            if not sd:
                birth(c, s, "a scout find")
                return
            d = min(sd, key=lambda x: hills.distance(x.centre, c))
            p = dict(c.params)
            for n in names():
                spec = space.specs[n]
                if spec.allowed or spec.type == "bool":
                    continue
                u = (space._to_unit(spec, c.params.get(n, spec.default)) + space._to_unit(spec, d.centre.params.get(n, spec.default))) / 2
                p[n] = space._from_unit(spec, u)
            mid = space._assemble(p)
            if mid is not None and fresh(mid) and propose(Unit("findtest", mid), [], 5):
                find_test[mid.key] = (c, s, d.id)

        def judge_find(key: str, mid_score: Optional[float]) -> None:
            item = find_test.pop(key, None)
            if item is None:
                return
            c, s, did = item
            d = domains.get(did)
            ref = d.best if d is not None else s
            if mid_score is None or mid_score < min(s, ref) - ACCEPT_MIN:          # a valley: a hill of its own
                if elapsed() < T(REDUCE_TO_3[0]) and len(live()) + len(merges) < n_domains:
                    birth(c, s, f"a scout find, valley at {mid_score if mid_score is not None else 'a failed run'}")
            elif d is not None and status[d.id] == "live" and d.offer(c, s, clock()):
                log.info("engine: domain %d -> %.4f, a scout find on its hill (midpoint %.4f) t+%.1fm", d.id, s,
                         mid_score, elapsed() / 60)

        def judge_midpoint(key: str, score: Optional[float]) -> None:
            pair = midpoint.pop(key, None)
            if pair is None:
                return
            a, b = domains.get(pair[0]), domains.get(pair[1])
            if a is None or b is None or status[a.id] != "live" or status[b.id] != "live":
                return
            if score is not None and score >= min(a.best, b.best) - ACCEPT_MIN:
                start_merge([a, b], f"one hill by midpoint {score:.4f}")

        # ---- the controller ------------------------------------------------------------------------------------------
        def control() -> None:
            now = clock()
            # merges that are back
            for mj in list(merges.values()):
                if mj.back >= MERGE_BACK * mj.sent:
                    resolve(mj)
            # frozen domains: merge into a live domain on the same hill, else retire (their best stays a candidate)
            for d in [x for i, x in domains.items() if status[i] == "live" and x.state == "exhausted"]:
                others = [x for x in live() if x is not d]
                near = min(others, key=lambda x: hills.distance(x.centre, d.centre), default=None)
                if near is not None and hills.same(near.centre, d.centre):
                    start_merge([near, d], "a frozen domain on a live one's hill")
                else:
                    status[d.id] = "retired"
                    log.info("engine: domain %d retired at %.4f (its climb is exhausted) t+%.1fm", d.id, d.best, elapsed() / 60)
            # evidence: two domains on one hill
            lv = [d for d in live() if d.state != "final"]
            busy: Set[int] = set()
            for dist, a, b in same_hill_pairs(lv):
                if a.id in busy or b.id in busy:
                    continue
                if dist < MIN_FAR:
                    start_merge([a, b], f"one hill, {dist} settings apart")
                    busy |= {a.id, b.id}
                elif abs(a.best - b.best) <= MERGE_CLOSE and dist < MIN_FAR + 2 and \
                        now - tested_pair.get((a.id, b.id), -1e9) > MIDPOINT_RETEST_S:
                    test_midpoint(a, b)
            # schedule: at most 3 from REDUCE_TO_3[1]
            if elapsed() >= T(REDUCE_TO_3[0]) and not st["r3"]:
                st["r3"] = True
                lv = live()
                count = len(lv) + len(merges)
                taken: Set[int] = set()
                for dist, a, b in same_hill_pairs(lv):
                    if count <= 3:
                        break
                    if a.id in taken or b.id in taken:
                        continue
                    start_merge([a, b], f"scheduled 6 to 3, {dist} settings apart")
                    taken |= {a.id, b.id}
                    count -= 1
            if elapsed() >= T(REDUCE_TO_3[1]) and not st["r3_done"]:
                st["r3_done"] = True
                for mj in list(merges.values()):
                    resolve(mj)
                lv = live()
                if len(lv) > 3:
                    holder = max(lv, key=lambda d: d.best)
                    rest = sorted((d for d in lv if d is not holder), key=lambda d: (-d.progress(now), -d.best))
                    for d in rest[2:]:
                        status[d.id] = "retired"
                        log.info("engine: domain %d (%.4f, progress %.4f/min) dropped to keep 3 t+%.1fm", d.id, d.best,
                                 d.progress(now), elapsed() / 60)
            # schedule: one domain for the final
            if elapsed() >= T(REDUCE_TO_1[0]) and not st["r1"]:
                st["r1"] = True
                for mj in list(merges.values()):
                    resolve(mj)
                lv = live()
                if len(lv) > 1:
                    start_merge(lv, "scheduled 3 to 1")
            # free slots, until the scheduled merges begin: the scout's finds on hills nobody climbs
            if st["founded"] and elapsed() < T(REDUCE_TO_3[0]) and len(live()) + len(merges) < n_domains:
                ranked = sorted(candidates.values(), key=lambda t: -t[1])
                cs = centres(standing())
                weakest = min((d.best for d in live()), default=float("-inf"))
                for c, s in ranked:
                    if s < max(top() - WORTH_A_HILL, weakest) or len(find_test) >= 4:
                        break                                  # a find below every domain is not worth a domain
                    if c.key in was_centre or c.key in find_done or hills.owner(c, cs) is not None:
                        continue
                    test_find(c, s)                            # a new hill only if a valley separates it
            # sharing between domains: their strongest settings, offered as moves
            if now - st["last_share"] >= SHARE_EVERY_S and len(live()) > 1:
                st["last_share"] = now
                for d in live():
                    for a in sorted(d.axis_gain, key=lambda x: -d.axis_gain[x])[:SHARE_AXES]:
                        v = d.centre.params.get(a)
                        for e in live():
                            if e is not d and not space._same(space.specs[a], e.centre.params.get(a), v):
                                e.queue("share", {a: v})
            # the final
            t = left()
            if not st["final"] and st["founded"] and (t is not None and t < T(FINAL_S) or elapsed() >= T(REDUCE_TO_1[1])):
                final_round()
            if st["final"] and not st["picked"] and t is not None and t < PICK_S:
                pick()

        def final_round() -> None:
            st["final"] = True
            for mj in list(merges.values()):
                resolve(mj)
            sd = standing()
            if not sd:
                return
            win = max(sd, key=lambda d: d.best)
            others = sorted((d for d in sd if d is not win), key=lambda d: -d.best)
            for d in others:
                status[d.id] = "retired"
            if st["runner_up"] is None and others:
                st["runner_up"] = others[0].centre
            status[win.id] = "live"
            win.state = "climb"
            win.start_final()
            if st["runner_up"] is not None:
                ru = st["runner_up"]
                diff = genetic.differing(space, names(), win.centre, ru)
                if 0 < len(diff) <= genetic.EXHAUSTIVE_MAX:
                    for c in genetic.merge_children(space, win.centre, ru, rng, names(), n_mutated=0):
                        win.queue("pairing", {n: c.params.get(n) for n in diff})
            log.info("engine: final at t+%.1fm on domain %d (%.4f)%s", elapsed() / 60, win.id, win.best,
                     ", with a pairing against the runner-up" if st["runner_up"] is not None else "")

        def pick() -> None:
            """The plateau-centre pick: among configs within PICK_BAND of the best, the one whose neighbours (at most
            PICK_NEAR settings apart) score highest - the middle of the plateau, least exposed to the answer set's blur."""
            st["picked"] = True
            ranked = sorted(candidates.values(), key=lambda t: -t[1])
            if not ranked:
                return
            best = ranked[0][1]
            near_top = [(c, s) for c, s in ranked if s >= best - PICK_BAND][:60]
            pool = ranked[:600]
            chosen, chosen_v = ranked[0][0], None
            for c, s in near_top:
                nb = [s2 for c2, s2 in pool if c2.key != c.key and far_axes(space, space.names, c, c2) <= PICK_NEAR]
                if len(nb) < PICK_MIN:
                    continue
                v = sum(nb) / len(nb)
                if chosen_v is None or v > chosen_v:
                    chosen, chosen_v = c, v
            client.pick_key = chosen.key
            log.info("engine: plateau-centre pick %s (%.4f; best %.4f; neighbour mean %s) t+%.1fm", chosen.key,
                     candidates.get(chosen.key, (None, 0.0))[1], best,
                     f"{chosen_v:.4f}" if chosen_v is not None else "n/a", elapsed() / 60)

        # ---- the scout ------------------------------------------------------------------------------------------------
        def scout_pick(n: int) -> List[Config]:
            cs = centres(standing())
            pool = [c for c in latin_hypercube(space, rng, 2 * n) if fresh(c) and hills.owner(c, cs) is None]
            cents = list(cs.values())
            for _ in range(n):                                 # crossovers between the bests of different hills
                if len(cents) < 2:
                    break
                a, b = rng.sample(cents, 2)
                c = genetic.uniform_child(space, a, b, rng, names())
                if c is not None and rng.random() < genetic.MUTATE_FRAC:
                    c = genetic.mutate(space, c, rng, names()) or c
                if fresh(c):
                    pool.append(c)
            if ranker.ready and pool:
                mean, spread = ranker.predict(pool)
                order = sorted(range(len(pool)), key=lambda i: -(mean[i] + SCOUT_KAPPA * spread[i]))
                pool = [pool[i] for i in order]
            return pool[:n]

        def send_scout(idle: Optional[int]) -> bool:
            if not scout_queue:
                scout_queue.extend(scout_pick(max(1, min(64, idle or 64))))
            while scout_queue:
                c = scout_queue.pop(0)
                if fresh(c) and send_global(c, "scout", RIDERS_GLOBAL, 1):
                    st["scouts"] += 1
                    return True
            return False

        # ---- allocation ------------------------------------------------------------------------------------------------
        def fill() -> int:
            n = 0
            now = clock()
            while not client.exhausted:
                idle = client.idle_slots()
                if idle is not None and idle <= 0:
                    break
                if idle is None and len(pending) >= cap:
                    break
                if not st["founded"]:
                    break                                      # the opening's waves are queued by opening_step
                if st["scouts"] < max(1, int(cap * SCOUT_SHARE)) and send_scout(idle):
                    n += 1
                    continue
                lv = live()
                if lv:
                    rates = [max(d.progress(now), 1e-4) for d in lv]
                    floor = 0.25 * sum(rates) / len(rates)
                    w = [max(x, floor) for x in rates]
                    share = cap * (1 - SCOUT_SHARE)
                    targets = {d.id: share * x / sum(w) for d, x in zip(lv, w)}
                    for d in lv:
                        if not d.pool and d.inflight < targets[d.id]:
                            d.refill(names(), rng, int(targets[d.id] - d.inflight) + 1)
                    ready = [d for d in lv if d.pool]
                    if ready:
                        d = max(ready, key=lambda x: targets[x.id] - x.inflight)
                        move = d.take()
                        kids = d.siblings(move, rng)
                        if propose(Unit("move", move.config, d.id, move), kids, 3):
                            n += 1
                        else:
                            d.dropped(move)
                        continue
                if send_scout(idle):                          # idle capacity never waits: the scout takes it
                    n += 1
                    continue
                break
            return n

        # ---- results -------------------------------------------------------------------------------------------------
        def note(unit: Unit, is_head: bool, o: Outcome) -> None:
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
                if len(candidates) > 1200:
                    for key in sorted(candidates, key=lambda x: -candidates[x][1])[600:]:
                        del candidates[key]
            if unit.role in ("wave", "scout"):
                if is_head:
                    if unit.role == "wave":
                        wave_back[unit.wave] += 1
                    else:
                        st["scouts"] = max(0, st["scouts"] - 1)
                if st["founded"] and not failed:
                    i = owner_of(o.config)
                    if i is not None and domains[i].offer(o.config, score, now):
                        log.info("engine: domain %d -> %.4f, a global point on its hill t+%.1fm", i, score, elapsed() / 60)
                return
            if unit.role == "midpoint":
                judge_midpoint(unit.head.key, None if failed else score)
                return
            if unit.role == "findtest":
                judge_find(unit.head.key, None if failed else score)
                return
            if unit.role == "merge":
                mj = merges.get(unit.merge)
                if mj is not None:
                    mj.back += 1
                    if not failed and (mj.best is None or score > mj.best[1]):
                        mj.best = (o.config, score)
                return
            d = domains.get(unit.domain) if unit.domain is not None else None
            if d is None or unit.move is None or status[d.id] == "retired":
                return
            if d.consider(unit.move, is_head, o.config, score, failed, now):
                log.info("engine: domain %d -> %.4f via %s %s t+%.1fm", d.id, d.best, unit.move.kind,
                         ",".join(unit.move.axes), elapsed() / 60)

        def maybe_fit() -> None:
            if ranker.stale() and (fitting[0] is None or fitting[0].done()):
                fitting[0] = asyncio.create_task(asyncio.to_thread(ranker.fit))

        def status_line() -> None:
            st["last_log"] = clock()
            lv = [d for i, d in domains.items() if status[i] != "retired"]
            log.info("engine: t+%.1fm domains %s | merges %s | scouts %d | candidates %d | pending %d", elapsed() / 60,
                     " ".join(f"{d.id}:{d.best:.3f}/{d.state[:5]}/{d.mode[:1]}/{d.inflight}{'' if status[d.id] == 'live' else '/m'}"
                              for d in lv), " ".join(f"{m.id}:{m.back}/{m.sent}" for m in merges.values()) or "-",
                     st["scouts"], len(candidates), len(pending))
            client.notes["engine"] = {
                "version": "v4", "domains": [{"id": d.id, "best": round(d.best, 4), "state": d.state, "mode": d.mode,
                                              "inflight": d.inflight, "status": status[d.id],
                                              "born_min": round((d.born - t0) / 60, 1)} for d in lv],
                "merges": len(merges), "scouts": st["scouts"], "candidates": len(candidates), "pending": len(pending),
                "caps": caps, "dead_axes": sorted(set(book.names) - set(book.live(clock())))}

        # ---- the loop ------------------------------------------------------------------------------------------------
        sd = seeds()
        cells = [c for c in latin_hypercube(space, rng, max(0, cap - len(sd))) if all(c.key != x.key for x in sd)]
        log.info("engine: capacity %d -> up to %d domains; schedule x%.2f (a %d-minute round)", cap, n_domains, scale_t,
                 round(ROUND_S * scale_t / 60))
        launch_wave(sd + cells, RIDERS_WAVE0, "seeds + hypercube")
        while True:
            if pending:
                done, _ = await asyncio.wait(list(pending), timeout=TICK_S, return_when=asyncio.FIRST_COMPLETED)
                for fut in done:
                    unit, is_head = pending.pop(fut)
                    try:
                        note(unit, is_head, fut.result())
                    except Exception:  # noqa: BLE001
                        log.exception("engine: result handling failed")
                maybe_fit()
            try:
                opening_step()
                if st["founded"]:
                    control()
            except Exception:  # noqa: BLE001 - a controller bug must not end the round
                log.exception("engine: controller step failed")
            n = fill()
            if clock() - st["last_log"] > 120:
                status_line()
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
