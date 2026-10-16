"""The search engine (v4.22): twenty-four centres, then twelve, then one.
Owner: search.

v4.22 keeps v4.21. The opening founds 24 centres. Halfway to the final, the 12 that have moved the most keep
the workers. The final still starts at minute 56. The strongest of those twelve keeps the seat, and each of the
others hands it every setting value that was in place when that basin's score rose.

The search engine (v4.21): the final starts at minute 56 and takes the steps that raised the other basins.
Owner: search.

v4.21 keeps v4.20. The final starts at minute 56. The strongest of the six basins keeps the seat. Each of the
other five hands it every setting value that was in place when that basin's score rose, one setting at a time,
including a value it later replaced.

The search engine (v4.20): a move keeps core and the false-positive count together.
Owner: search.

v4.20 keeps v4.19. The opening still breeds high core and high false positives apart, then crosses them.
Every seat, including the core seat, keeps a move in one of two ways: core rises and the false-positive
count does not, or the false-positive count falls and core holds. A higher blended score that adds false
positives is refused. A cleaner score that gives up core is refused. A refused config still hands on the
single setting that raised core or lowered the false-positive count.

The search engine (v4.19): a config that loses still keeps the settings that helped.
Owner: search.

v4.19 keeps v4.18. A scored config the centre does not take can still hold one setting that raised core.
That setting is tried, alone, on the basins that stay. At the cut, each retired basin hands over the same
kind of setting: the values that raised core, and the centre's value on every axis whose moves gained.
The settings that did not gain stay behind with the config.

The search engine (v4.18): a scan that failed is sent again.
Owner: search.

v4.18 keeps v4.17. A full-range scan whose every result failed (a full disk, a crash) is not a flat peak.
The basin sends that scan again, up to two extra times, and only then holds the peak. A scan that scored
and gained nothing still holds, as before.

The search engine (v4.17): the opening generations actually breed core and false positives apart.
Owner: search.

v4.17 keeps v4.16. With islands on, which is the default, breed() returned the island children before the
core and false-positive split, so generations 1–3 never used it. Those three generations now take the split
first. A breed with no generation number still uses islands.

The search engine (v4.16): the opening breeds core and false positives apart, then mixes them.
Owner: search.

v4.16 keeps v4.15. Generations 1 and 2 each split in half. One half is children of the highest measured core,
the other half children of the highest false-positive rate. Generation 3 crosses one parent from each group.
The portfolio keeps the same three kinds of config from a reviewed run: near-best score, high core, and high
false positives, and wave 0 draws from all three. After the cut, a stall still retires a basin that has
stopped gaining, but never the live basin with the highest score: its workers stay until the final.

The search engine (v4.15): one founding seat chases the variants the clump missed.
Owner: search.

v4.15 keeps v4.14. Eleven seats are still the highest centres far enough apart. The twelfth is the scored
opening config at least domains.found_min_far settings from those eleven with the highest measured core, and
its blended score may sit below the clump. That seat keeps a move when the score rises, including when false
positives rise. The other seats still refuse a higher score that adds false positives. The midpoint cut is
unchanged, so a seat whose score never moves loses its workers.

The search engine (v4.14): the first scan's winners are tried together.
Owner: search.

v4.14 keeps v4.13. When a basin's first scan comes back, the centre still moves to the best point that raises
the score without adding false positives. The scan's other winners that passed that test are applied together
as one config, and that config goes out before the neighbourhood walk. It is kept only when its own score rises
and its false positives do not.

The search engine (v4.13): a stall waits for the cut.
Owner: search.

v4.13 keeps v4.12. All twelve centres climb until the cut. A basin that has stopped gaining keeps its seat until
then, so the ranking still sees it. The basins the cut keeps start their stall clock over, so that ranking is not
undone in the next minutes. After the cut, a stall retires the basin and does not birth another domain. A higher
score moves a centre only when it does not add false positives: the climb keeps a cleaner call of the variants it
already finds, and a gain bought with more false positives is not a move.

The search engine (v4.12): twelve basins climb, then the moving six keep the workers.
Owner: search.

v4.12 keeps v4.11's one scan per basin. The opening births every centre it found, up to twelve, and one of them
is the scored config farthest from the others, down to 5 points below the opening best. Halfway from founding to
the final, once every first scan is back, the six with the highest score-plus-gain keep climbing. The rest retire
and their workers join those six. A coarse sweep that finishes goes on to the fine neighbourhood and does not
launch the wide all-pairs round.

The search engine (v4.11): one full-range scan per basin.
Owner: search.

v4.11 keeps v4.10. A basin scans every setting once, at founding. A later gain, including one of 0.05, moves
the centre and the neighbourhood climb continues. It does not start the scan again.

The search engine (v4.10): three opening generations, and the final breeds every domain.
Owner: search.

v4.10 keeps v4.9's climb. The opening runs opening.generations (3) generations after wave 0, and a small gain
does not end it before those three have been launched. At the final, one capped generation takes each setting
from any domain's centre, including domains retired earlier in the round. The keeper is still the strongest
basin still climbing. The other live basins still hand over their differing settings as moves.

The search engine (v4.9): the changes made while v4.8 was running.
Owner: search.

v4.9 keeps v4.8's six climbing basins, six held candidates, and the ten-minute stall. Two changes landed while
that replay was still going, and both are this version. The answer set separates neighbours at climb.accept_min
(0.005), so a scan or a chain that gains that much becomes the centre and the neighbourhood climb continues. A
gain of climb.decisive (0.05) is the one that restarts the full-range scan. For the five minutes before the
final (schedule.no_new_before_final_min), a stalled basin retires and its workers join the basins already
climbing. A new centre is not founded: its scan would not finish before the merge.

The search engine (v4.8): a basin that stops moving for ten minutes gives its workers to a new centre.
Owner: search.

v4.8 keeps v4.7's six climbing basins and six held candidates. After a basin's first scan, ten minutes
(domains.stall_min) without a new best retires it. The next held candidate takes the seat, or the next unused
centre at least six settings away when the held list is empty. Its best config stays a candidate for the pick.

The search engine (v4.7): six basins climb, six candidates wait, and a failure is swapped one at a time.
Owner: search.

v4.7 identifies up to twelve opening centres and starts six of them. The other six are held, not given workers.
A basin whose first scan gains under climb.decisive is retired. The next held candidate takes that seat, and the
one after it waits until that scan is back. The even split, the single merge at the final, and the raw-best pick
are the v4.6 behaviour this sits on.

The search engine (v4.6): the basins stay separate until one merge at the final.
Owner: search.

v4.6 keeps v4.5's twelve opening basins and starts all of them. There is no merge at minute 34 and none at minute 50.
Workers stay split evenly across the live
basins until the final window. When schedule.final_min remains, the strongest basin keeps its seat, the others
retire, and it takes their differing settings as moves plus one capped set of recombination children. Those
children go out before the final's halved steps and before any chain. The final is then iterated local search on
that one centre. The pick is still the highest local score.

The search engine (v4.5): twelve opening basins, and one replacement of those that do not move.
Owner: search.

v4.5 founds up to domains.max basins (12 on a fleet; one per eight workers on a small machine). Seats inside
domains.worth_a_hill are taken first; the rest may sit up to 5 points below the best, still domains.found_min_far
apart. While any basin's first scan is still out, workers are split evenly. A gain under climb.decisive does not
count as progress: that basin is retired once every first scan is back, or at schedule.reduce_to_3_start_min, and
the next unused candidate takes its seat. That happens once. The merge down to one domain at minute 50 and the
raw-best pick are unchanged from v4.4. The ledger from v4.4 is still written.

The search engine (v4.4): the v4.3 search, plus a ledger of every basin's first scan.
Owner: search.

v4.4 does not change founding, the scan, the 0.05 rule, progress-weighted workers, or the 50-minute merge.
It records, for each domain, the opening score, the best point on its first scan, and the score it holds when
the scheduled merge to one domain begins. Those config keys are copied onto the run manifest so the review
scores them on the released truth. The sections below are the v4.3 search this ledger sits on.

v4.3, on top of the v4.2 switches (common/tuning.py): the submission is the highest local score. A domain's first
sweep is the scan, and it is not interrupted by a gain under climb.decisive; a gain that large moves the centre and
the scan runs again. When the scan comes back under that bar the peak is finished: no coarser climb, and no new
domain is founded for the rest of the round. An off-hill find is recombined with the nearest domain, one find at a
time, and those children wait behind that domain's scan (priority of a climb move, not ahead of it). A scheduled
merge still retires the weaker domain; the survivor keeps the transplant moves and one capped set of recombination
children, also behind its scan. Queued children of a retired domain are cancelled. The scout sampler is off
(scout.share and scout.late_share are 0). Before any of that, three settings are reduced to the values this round's
BAM can tell apart (search/bamgroups.py). The sections below describe the v4.1 behaviour the older switches return to.

  opening     wave 0 = the defaults, the official seed, up to opening.portfolio_share of the workers from the portfolio
              (search/portfolio.py: the legacy pipeline's 749 configs plus our reviewed runs, the round's read-length
              profile first) and a Latin hypercube for the rest (4 threshold riders each), sent while the answer set
              is still being built; then opening.generations genetic generations with niching (genetic.py), each
              launched once opening.launch_back of the previous one is back or workers sit idle (no barrier). Stops
              early when a generation gains less than opening.min_gain; hard limit opening.max_min
  founding    the hill registry (hills.py) groups the good configs into hills; one domain per hill, up to domains.max.
              Once the ranker is fitted, the hill rule counts only the settings it finds important
              (domains.importance_share), so settings that barely change the calls do not make separate hills
  climb       every domain runs iterated local search (climber.py): a full-range scan of every setting first. A gain
              under climb.decisive waits until climb.back_frac of that scan is back, and a flat scan holds the peak.
              A gain of climb.decisive moves the centre and arms the scan again. After a flat scan the coarser sweeps
              and the chains do not run. Winning settings are shared between domains every sharing.every_min
  finds       an off-hill config is recombined with the nearest live domain (pairing plus mutated children, at most
              merges.find_children), one such find per domain at a time. The children are sent by fill() behind that
              domain's scan, at the priority of a climb move. Nothing founded after the opening is a new domain. A
              find tested once is never tested again; a point on a merging domain's hill is handed to the merged domain
  merges      one operation everywhere: pairing (every combination of the settings where two centres differ, when
              few) plus a genetic generation of the two (genetic.merge_children); the best of parents and children is
              the merged domain's centre, which inherits both parents' tested configs and steps. Merges happen
                - on evidence, at any time: two domains on one hill
                - on schedule: closest pairs first from schedule.reduce_to_3_start_min so at most 3 domains remain at
                  schedule.reduce_to_3_by_min; all remaining domains at schedule.reduce_to_1_start_min (a three-way
                  merge) so one remains for the final
  final       the last schedule.final_min on one point: halving single steps and all pairs, one last pairing with the
              runner-up, chains with small shakes; then the pick is the highest local score (client.pick_key; main
              submits it)

One controller pass after every batch of results (and at least every TICK_S) decides everything from facts: what is
in flight, each domain's progress rate (own gain per minute over the last 10 minutes), the registry. Schedule times
are fractions of the search's time from its start, set for a 72-minute round. Every tunable number is a setting in
common/tuning.py (read with tune(), changeable on main while a round runs); the constants below are not tunable.
"""
from __future__ import annotations

import asyncio
import collections
import json
import logging
import math
import os
import random
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Set, Tuple

from common.params import ParamRanges
from common.scoring import GERMLINE_FP_SCALE, V2_CORE_WEIGHT, V2_GERMLINE_WEIGHT
from common.tuning import tune
from common.types import FULL, Config
from search.v422.axes import AxisBook, latin_hypercube
from search.v422.base import Strategy
from search.v422.client import EvalClient, Outcome
from search.v422.families import genotyping_riders
from search.v422.ranker import Ranker
from search.v422.reductions import Reductions
from search.v422.bamgroups import bam_bounds
from search.v422.space import GENOTYPING_PARAMS, Space, bam_caps
from search.v422 import genetic, portfolio
from search.v422.climber import Climber, Move
from search.v422.hills import MIN_FAR, Hills

log = logging.getLogger(__name__)

def core_from_score(score: float, fp_per_target: float) -> float:
    """The core implied by a v2 score and an FP rate: the inverse of score_from_parts."""
    germline = math.exp(-float(fp_per_target) / GERMLINE_FP_SCALE)
    return (float(score) / 100.0 - V2_GERMLINE_WEIGHT * germline) / V2_CORE_WEIGHT


def core_seat(ranked: Sequence[Tuple[Config, float]], anchors: Sequence[Tuple[Config, float]],
              fp_of: Dict[str, float], hills: Hills, min_far: int):
    """The scored config at least `min_far` settings from every anchor with the highest measured core.

    Its blended score may sit below the anchors: a call that found more variants and paid for them with false
    positives is the seat the score ranking skips. None when no such config has a measured FP rate."""
    anchor_keys = {c.key for c, _ in anchors}
    best = None
    for c, s in ranked:
        if c.key in anchor_keys:
            continue
        fp = fp_of.get(c.key)
        if fp is None:
            continue
        if any(hills.distance(c, a) < min_far for a, _ in anchors):
            continue
        core = core_from_score(s, fp)
        if best is None or core > best[0] or (core == best[0] and s > best[1]):
            best = (core, s, c, fp)
    return best


ROUND_S = 72 * 60                    # the schedule below is written for this much search time; it scales with it
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


@dataclass
class Unit:
    role: str                        # wave | scout | move | merge | findmerge | midpoint
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
    guest: Optional[Tuple[Config, float]] = None      # a scout find merged into its parent (not a domain)


@dataclass
class FindMerge:
    """A scout find recombined with the nearest domain's centre while that domain keeps climbing (merges.on_evidence 0
    era): a child better than the domain moves its centre; nothing waits for the job."""
    id: int
    domain: int
    find: Tuple[Config, float]
    sent: int = 0
    back: int = 0


async def settled_capacity(client: EvalClient, workers: int, wait_s: float = 60.0, settle_s: float = 20.0) -> int:
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
        bounds = bam_bounds(r.bam, ranges, tool) if r.bam else None
        caps = {} if bounds is not None else (bam_caps(r.bam, r.chrom, r.start, r.end) if r.bam else {})
        profile = portfolio.bam_profile(r.bam, r.chrom, r.start, r.end) if r.bam else None
        space = Space(ranges, tool, caps=caps, bounds=bounds)
        reductions = Reductions()
        space.reductions = reductions
        if bounds is not None:
            log.info("bam groups: %s", bounds.report)
        book = AxisBook(space)
        ranker = Ranker(space, seed=seed)
        clock = client.clock
        t0 = clock()
        tl = client.time_left()                      # the schedule spans the whole search, the node wait included
        cap = await settled_capacity(client, workers)
        scale_t = (tl if tl is not None and tl > 0 else ROUND_S) / ROUND_S     # the schedule follows the round's length
        T = lambda s: s * scale_t
        n_domains = max(1, min(tune("domains.max"), cap // 8))
        at = lambda step: T(60 * tune(f"schedule.{step}_min"))       # a schedule time, scaled to this round
        virtual = hasattr(client.executor, "now")
        genotyping = [a for a in GENOTYPING_PARAMS if a in space.names]
        # MINOS_FOCUS_AXES=a,b and MINOS_FOCUS_CENTRE=params.json: one hour on those axes, from that config.
        # Unset, this is the normal search. The centre is the finished run's best; every other setting stays put.
        focus_axes = [a.strip() for a in os.environ.get("MINOS_FOCUS_AXES", "").split(",") if a.strip()]
        focus_axes = [a for a in focus_axes if a in space.names]
        focus_centre: Optional[dict] = None
        centre_path = os.environ.get("MINOS_FOCUS_CENTRE", "")
        if focus_axes and centre_path:
            with open(centre_path) as fh:
                loaded = json.load(fh)
            focus_centre = loaded.get("params", loaded) if isinstance(loaded, dict) else None

        pending: Dict[asyncio.Future, Tuple[Unit, bool]] = {}
        proposed: Set[str] = set()
        domains: Dict[int, Climber] = {}
        status: Dict[int, str] = {}                  # live | merging | retired
        merges: Dict[int, MergeJob] = {}
        find_merges: Dict[int, FindMerge] = {}
        # (role, config, find-merge id or 0) waiting until that domain's scan is back. findmerge | recomb
        behind: Dict[int, List[Tuple[str, Config, int]]] = {}
        ids = [0, 0]                                 # next domain id, next merge id
        candidates: Dict[str, Tuple[Config, float]] = {}
        fp_of: Dict[str, float] = {}                     # config key -> false positives per target, when measured
        core_of: Dict[str, float] = {}                   # config key -> core, when measured
        part_keep: Dict[str, List[Tuple[float, Config]]] = {"core": [], "fp": []}   # top measured core and FP rate
        island_of: Dict[str, int] = {}               # opening.islands: config key -> island
        was_centre: Set[str] = set()
        wave_heads: List[int] = []                   # heads sent per wave (0 = hypercube, 1.. = generations)
        wave_back: List[int] = []
        wave_best: List[float] = []                  # best score seen when each wave launched
        wave_started: List[float] = []               # seconds into the search when each wave was launched
        reserve: List[Tuple[Config, float]] = []
        st = {"founded": False, "founded_at": 0.0, "cut": False, "gens_done": False, "r3": False, "r3_done": False,
              "replaced": False, "r1": False,
              "owed": 0,
              "final": False,
              "last_share": 0.0, "last_log": t0, "last_publish": t0 - TICK_S, "scouts": 0, "runner_up": None, "picked": False}
        midpoint: Dict[str, Tuple[int, int]] = {}
        find_test: Dict[str, Tuple[Config, float, int, str]] = {}   # midpoint key -> (find, score, nearest domain, mode)
        deferred: List[Tuple[Config, float]] = []                # global points on a merging domain's hill
        find_done: Set[str] = set()                              # finds already tested
        tested_pair: Dict[Tuple[int, int], float] = {}
        scout_queue: List[Config] = []
        origin: Dict[int, str] = {}                              # domain -> why it was born
        lineage: Dict[int, List[int]] = {}                       # domain -> the domains it was merged from
        ended: Dict[int, float] = {}                             # domain -> seconds into the search it retired
        alloc: Dict[int, float] = {}                             # domain -> workers it is allotted
        events: collections.deque = collections.deque(maxlen=200)
        picked: Dict[str, object] = {}
        fitting: List[Optional[asyncio.Task]] = [None]
        wake = asyncio.Event()
        client.on_result.append(lambda _o: wake.set())

        def say(msg: str, *args) -> None:
            """A decision: logged, and kept for the dashboard's event feed."""
            log.info("engine: " + msg, *args)
            try:
                events.append([round(time.time(), 1), msg % args if args else msg])
            except (TypeError, ValueError):
                events.append([round(time.time(), 1), msg])

        def elapsed() -> float:
            return clock() - t0

        def left() -> Optional[float]:
            return client.time_left()

        def fresh(c: Optional[Config]) -> bool:
            return c is not None and c.key not in proposed and (c.key, FULL) not in client.latest

        def names() -> List[str]:
            live = book.live(clock()) + genotyping
            if not focus_axes:
                return live
            return [a for a in focus_axes if a in live or a in space.names]

        important: Dict[str, object] = {"at": -1, "names": set()}

        def hill_names() -> List[str]:
            """The settings the hill rule counts: once the ranker is fitted, those whose importance is at least
            domains.importance_share of the mean (settings that barely change the calls no longer make two configs
            separate hills); all live settings before that, or when too few would be left."""
            base = names()
            share = tune("domains.importance_share")
            if share <= 0 or not ranker.ready:
                return base
            if important["at"] != ranker.fitted_on:
                imp = ranker.importance() or {}
                mean = sum(imp.get(n, 0.0) for n in base) / max(1, len(base))
                important.update(at=ranker.fitted_on, names={n for n in base if imp.get(n, 0.0) >= share * mean})
            keep = [n for n in base if n in important["names"]]
            return keep if len(keep) >= MIN_FAR + 2 else base

        hills = Hills(space, hill_names)

        def window() -> bool:
            """Scout finds may found domains freely until the scheduled merges begin."""
            return elapsed() < at("reduce_to_3_start")

        def minimum() -> int:
            """The fewest domains this stage of the schedule keeps (refilled from finds)."""
            if st["final"]:
                return 1
            if window():
                return tune("domains.min_early")
            return tune("domains.min_late")

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
            small = unit.role in ("move", "merge", "midpoint", "findtest", "findmerge", "recomb") \
                or (unit.role == "wave" and unit.wave > 0)
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

        def cancel_behind(domain_id: int, roles: Sequence[str] = ("findmerge", "recomb")) -> None:
            """Drop children not yet sent, and cancel the ones still in the executor's queue. A run that has
            already started is left to finish: its score can still be the one that is submitted."""
            q = behind.get(domain_id)
            if q is not None:
                behind[domain_id] = [item for item in q if item[0] not in roles]
                if not behind[domain_id]:
                    behind.pop(domain_id, None)
            targets = []
            for fut, (unit, _) in pending.items():
                if unit.domain != domain_id or unit.role not in roles or fut.done():
                    continue
                for jid, (f, _, _) in client._pending.items():
                    if f is fut:
                        targets.append((fut, unit, jid))
                        break
            if targets:
                client.executor.cancel(jid for _, _, jid in targets)
                queued = set()
                for item in getattr(client.executor, "_queue", ()) or ():
                    job = item[2] if isinstance(item, tuple) and len(item) >= 3 else item
                    jid = getattr(job, "job_id", None)
                    if jid:
                        queued.add(jid)
                for fut, unit, jid in targets:
                    if jid in queued and not fut.done():
                        fut.set_result(Outcome(unit.head, None, "cancelled", ""))
            if "findmerge" in roles:
                for fm in list(find_merges.values()):
                    if fm.domain == domain_id:
                        find_merges.pop(fm.id, None)

        def retire(d: Climber) -> None:
            status[d.id] = "retired"
            cancel_behind(d.id)

        # ---- opening --------------------------------------------------------------------------------------------
        def seeds() -> List[Config]:
            """The defaults, the official seed and up to opening.portfolio_share of wave 0 from the portfolio
            (search/portfolio.py: the round's read-length profile first). A focused run sends only its centre."""
            if focus_axes and focus_centre:
                c = space._assemble(dict(focus_centre))
                if c is not None:
                    say("focus centre %s on %s", c.key, ",".join(focus_axes))
                    return [c]
                say("focus centre did not assemble; opening as usual")
            out: List[Config] = []
            for p in ({}, OFFICIAL_SEED):
                c = space._assemble(dict(p))
                if c is not None and all(c.key != x.key for x in out):
                    out.append(c)
            want, n = int(cap * tune("opening.portfolio_share")), 0
            for c in portfolio.pick(space, profile, want):
                if all(c.key != x.key for x in out):
                    out.append(c)
                    n += 1
            say("seeds: %d from the portfolio (profile %s)", n, profile or "unknown")
            return out

        def launch_wave(configs: Sequence[Config], riders: int, label: str) -> None:
            w = len(wave_heads)
            wave_heads.append(0); wave_back.append(0); wave_best.append(top()); wave_started.append(elapsed())
            for c in configs:
                if send_global(c, "wave", riders, 5 if w == 0 else 4, wave=w):
                    wave_heads[w] += 1
            say("wave %d (%s) at t+%.1fm: %d points", w, label, elapsed() / 60, wave_heads[w])

        def opening_step() -> None:
            if st["founded"]:
                return
            if focus_axes and wave_heads and wave_back[0] >= wave_heads[0] and candidates:
                st["gens_done"] = True
                found()
                return
            w = len(wave_heads) - 1
            back = wave_back[w] / max(1, wave_heads[w])
            idle = client.idle_slots()
            if elapsed() >= T(60 * tune("opening.max_min")):                  # the opening's limit is hard: found now
                st["gens_done"] = True
                found()
                return
            if not st["gens_done"] and tune("opening.steady"):
                breed_steady(idle)
                return
            if not st["gens_done"]:
                can = back >= tune("opening.launch_back") or (back >= 0.5 and idle is not None and idle >= tune("opening.launch_idle") * cap)
                if can and len(candidates) >= 4:
                    if w >= tune("opening.generations") and top() - wave_best[w] < tune("opening.min_gain"):
                        st["gens_done"] = True
                        say("generations stop after %d (last gained %.3f) t+%.1fm", w, top() - wave_best[w],
                                 elapsed() / 60)
                    else:
                        children = [c for c in breed(cap, gen=w + 1) if fresh(c)]
                        launch_wave(children, tune("riders.global"), f"generation {w + 1}")
                        st["gens_done"] = len(wave_heads) - 1 >= tune("opening.generations")
                return
            if back >= tune("opening.found_back") or not pending:
                found()

        def remember(kind: str, metric: float, config: Config) -> None:
            board = part_keep[kind]
            for i, (old, c) in enumerate(board):
                if c.key == config.key:
                    if metric > old:
                        board[i] = (metric, config)
                    return
            board.append((metric, config))
            board.sort(key=lambda t: -t[0])
            del board[40:]

        def breed_parts(n: int, gen: int) -> List[Config]:
            """Half the children from the highest core, half from the highest false-positive rate."""
            core = [(c, m) for m, c in part_keep["core"]]
            fp = [(c, m) for m, c in part_keep["fp"]]
            if len(core) < 2 or len(fp) < 2:
                say("generation %d: blended, only %d core and %d false-positive parents", gen, len(core), len(fp))
                return genetic.generation(space, list(candidates.values()), rng, n, hills, names())
            half = n // 2
            a = genetic.generation(space, core, rng, half, hills, names())
            b = genetic.generation(space, fp, rng, n - half, hills, names())
            say("generation %d: %d children of the high core, %d of the high false positives", gen, len(a), len(b))
            return a + b

        def breed_cross(n: int) -> List[Config]:
            """One child takes settings from a high-core parent and a high-FP parent."""
            cores = [c for _, c in part_keep["core"]]
            fps = [c for _, c in part_keep["fp"]]
            kids = genetic.cross_generation(space, cores, fps, rng, n, names())
            if len(kids) < 4:
                say("generation 3: blended, the cross made %d children", len(kids))
                return genetic.generation(space, list(candidates.values()), rng, n, hills, names())
            say("generation 3: %d children of a core parent and an fp parent", len(kids))
            return kids

        def breed(n: int, gen: Optional[int] = None) -> List[Config]:
            """`n` children of the scored configs: on islands (opening.islands) or with soft niching.
            Generations 1 and 2 split core from false positives; generation 3 crosses the two.
            That split wins over islands, which otherwise breed every generation themselves."""
            if gen in (1, 2):
                return breed_parts(n, gen)
            if gen == 3:
                return breed_cross(n)
            if tune("opening.islands"):
                kids = genetic.island_generation(space, list(candidates.values()), island_of, n_domains, rng, n, hills,
                                                 names(), tune("domains.found_min_far"))
                for c, i in kids:
                    island_of[c.key] = i
                return [c for c, _ in kids]
            return genetic.generation(space, list(candidates.values()), rng, n, hills, names())

        def breed_steady(idle: Optional[int]) -> None:
            """opening.steady: children go to idle workers as they free up, opening.generations x the workers in all;
            stops early when a workers' worth of children gains less than opening.min_gain."""
            if len(wave_heads) < 2:
                if wave_back[0] < 0.5 * max(1, wave_heads[0]) or len(candidates) < 4:
                    return
                launch_wave([], 0, "bred into idle workers")
                st["bred"], st["mark"], st["marks"] = 0, top(), 0
            if st["bred"] >= tune("opening.generations") * cap:
                st["gens_done"] = True
                say("breeding done: %d children t+%.1fm", st["bred"], elapsed() / 60)
                return
            if st["bred"] >= cap * (st["marks"] + 1):
                st["marks"] += 1
                if st["marks"] >= tune("opening.generations") and top() - st["mark"] < tune("opening.min_gain"):
                    st["gens_done"] = True
                    say("breeding stops after %d children (the last %d gained %.3f) t+%.1fm", st["bred"], cap,
                        top() - st["mark"], elapsed() / 60)
                    return
                st["mark"] = top()
            if idle is None or idle < max(8, int(0.02 * cap)):
                return
            for c in breed(idle):
                if fresh(c) and send_global(c, "wave", tune("riders.global"), 4, wave=1):
                    wave_heads[1] += 1
                    st["bred"] += 1

        def found() -> None:
            st["founded"] = True
            ranked = sorted(candidates.values(), key=lambda t: -t[1])
            band, min_far = tune("domains.worth_a_hill"), tune("domains.found_min_far")
            if tune("opening.islands") and island_of:
                firsts: Dict[int, Tuple[Config, float]] = {}       # each island's best, then everything else
                for c, s in ranked:
                    i = island_of.get(c.key)
                    if i is not None and i not in firsts:
                        firsts[i] = (c, s)
                order = sorted(firsts.values(), key=lambda t: -t[1]) + ranked
            else:
                order = ranked
            best = ranked[0][1] if ranked else 0.0
            near: List[Tuple[Config, float]] = []
            for c, s in order:
                if len(near) >= n_domains or s < best - band:
                    continue
                if all(c.key != o.key and hills.distance(c, o) >= min_far for o, _ in near):
                    near.append((c, s))
            picked = list(near)
            if len(picked) < n_domains:
                picked += hills.take(ranked, n_domains - len(picked), min_far, best - 5.0, [c for c, _ in picked])
            near_keys = {c.key for c, _ in near}
            far_key = None
            core_key = None
            if n_domains >= 2 and picked:
                anchors = picked[:n_domains - 1] if len(picked) >= n_domains else list(picked)
                seat = core_seat(ranked, anchors, fp_of, hills, min_far) if anchors else None
                if seat is not None:
                    core, s, c, fp = seat
                    picked = list(anchors) + [(c, s)]
                    core_key = c.key
                    say("core seat %.4f core, fp %.3f, score %.4f, %d settings from the clump",
                        core, fp, s, min(hills.distance(c, a) for a, _ in anchors))
                elif len(picked) >= 2:                    # no measured core outside: the farthest config within 5 points
                    anchors = picked[:-1]
                    anchor_keys = {c.key for c, _ in anchors}
                    current = picked[-1]
                    cur_dist = min(hills.distance(current[0], a) for a, _ in anchors)
                    best_far = None
                    best_dist = cur_dist
                    for c, s in ranked:
                        if c.key in anchor_keys or s < best - 5.0:
                            continue
                        dist = min(hills.distance(c, a) for a, _ in anchors)
                        if dist <= cur_dist:
                            continue
                        if best_far is None or dist > best_dist or (dist == best_dist and s > best_far[1]):
                            best_dist = dist
                            best_far = (c, s)
                    if best_far is not None and best_far[0].key != current[0].key:
                        picked = anchors + [best_far]
                        far_key = best_far[0].key
            for c, s in picked:
                if c.key == core_key:
                    why = "the opening, core outside the clump"
                elif c.key == far_key:
                    why = "the opening, far from the clump"
                elif c.key in near_keys:
                    why = "the opening"
                else:
                    why = "the opening, past the band"
                birth(c, s, why, allow_fp=c.key == core_key)
            st["founded_at"] = elapsed()
            say("founding at t+%.1fm: %d domains, %d held (%d within %.2f of the best) from %d scored configs (waves %s)",
                elapsed() / 60, len(live()), len(reserve), len(near), band, len(candidates),
                " ".join(f"{b}/{h}" for b, h in zip(wave_back, wave_heads)))

        def birth(c: Config, s: float, why: str, parents: Sequence[Climber] = (), allow_fp: bool = False) -> Climber:
            d = Climber(ids[0], space, book, c, s, clock(), owner_of, extra_axes=genotyping, scan=not parents,
                        allow_fp=allow_fp)
            ids[0] += 1
            for p in parents:
                d.inherit(p.tried, p.steps, p.axis_gain)
                if p.mode == "fine":
                    d.mode = "fine"
            domains[d.id], status[d.id] = d, "live"
            origin[d.id], lineage[d.id] = why, [p.id for p in parents]
            if st["final"]:
                d.start_final()                              # born from a merge during the final: it stays final
            was_centre.add(c.key)
            if c.key in fp_of:
                d.best_fp = fp_of[c.key]
            if c.key in core_of:
                d.best_core = core_of[c.key]
            say("domain %d born at %.4f (%s) t+%.1fm", d.id, s, why, elapsed() / 60)
            return d

        # ---- merges -----------------------------------------------------------------------------------------------
        def start_merge(parents: List[Climber], why: str, guest: Optional[Tuple[Config, float]] = None) -> None:
            parents = sorted(parents, key=lambda d: -d.best)
            if tune("merges.retire"):
                retire_merge(parents, why, guest)
                return
            for p in parents:
                status[p.id] = "merging"
                p.pool = []
            first = (parents[0].centre, parents[0].best)
            if guest is not None and guest[1] > first[1]:
                first = guest
            mj = MergeJob(ids[1], [p.id for p in parents], why, started=clock(), best=first, guest=guest)
            ids[1] += 1
            merges[mj.id] = mj
            if guest is not None:
                kids = genetic.merge_children(space, parents[0].centre, guest[0], rng, names())
            elif len(parents) == 2:
                kids = genetic.merge_children(space, parents[0].centre, parents[1].centre, rng, names())
            else:
                kids = genetic.merge3_children(space, [p.centre for p in parents], rng, names())
            tried = set().union(*(p.tried for p in parents))
            kids = [c for c in kids if c.key not in tried and fresh(c)]
            limit = max(tune("merges.children_min"), int(cap * tune("merges.children_share")))
            if len(kids) > limit:
                kids = rng.sample(kids, limit)
            for c in kids:
                if c.key in tried or not fresh(c):
                    continue
                if propose(Unit("merge", c, merge=mj.id), [], 6):
                    mj.sent += 1
            say("merge %d (%s): domains %s at %s, %d children, t+%.1fm", mj.id, why,
                     ",".join(str(p.id) for p in parents), " ".join(f"{p.best:.4f}" for p in parents), mj.sent,
                     elapsed() / 60)
            if mj.sent == 0:
                resolve(mj)

        def retire_merge(parents: List[Climber], why: str, guest: Optional[Tuple[Config, float]]) -> None:
            """v3.2's merge: the strongest domain stays, the others retire and hand it their differing settings as
            moves (a find merged in: the domain moves to it when it is better, and the other point's settings
            become moves). No children are run."""
            keep, rest = parents[0], parents[1:]
            n = 0
            if guest is not None:
                old = keep.centre
                if keep.offer(guest[0], guest[1], clock(), fp_of.get(guest[0].key), core_of.get(guest[0].key)):
                    n += keep.transplant(old, names())
                else:
                    n += keep.transplant(guest[0], names())
            for p in rest:
                retire(p)
                keep.inherit(p.tried, p.steps, p.axis_gain)
                n += keep.transplant(p.centre, names())
                n += keep.hand(p.raised, [keep])          # every value that raised this basin, including earlier steps
            if rest:
                st["runner_up"] = rest[0].centre
            pool: List[Tuple[Config, float]] = []
            seen: Set[str] = set()
            for d in sorted(domains.values(), key=lambda d: -d.best):
                if d.centre.key in seen:
                    continue
                seen.add(d.centre.key)
                pool.append((d.centre, d.best))
            n_kids = 0
            if len(pool) >= 2:
                tried = keep.tried
                limit = max(tune("merges.children_min"), int(cap * tune("merges.children_share")))
                kids = [k for k in genetic.across_children(space, [c for c, _ in pool], rng, names(), limit)
                        if k.key not in tried and fresh(k)]
                behind.setdefault(keep.id, []).extend(("recomb", k, 0) for k in kids)
                n_kids = len(kids)
            n_retired = sum(1 for s in status.values() if s == "retired")
            say("merge (%s): domain %d keeps %.4f, %s; %d transplant moves, %d children from %d centres "
                "(%d retired) behind its scan t+%.1fm",
                why, keep.id, keep.best,
                ", ".join(f"domain {p.id} ({p.best:.4f}) retired" for p in rest) or "a find handed over",
                n, n_kids, len(pool), n_retired, elapsed() / 60)

        def resolve(mj: MergeJob) -> Climber:
            merges.pop(mj.id, None)
            parents = [domains[i] for i in mj.parents]
            pool = [(p.centre, p.best, f"domain {p.id}") for p in parents]
            if mj.guest is not None:
                pool.append((mj.guest[0], mj.guest[1], "the find"))
            if mj.best is not None and all(mj.best[0].key != x[0].key for x in pool):
                pool.append((mj.best[0], mj.best[1], "a child"))
            c, s, src = max(pool, key=lambda t: t[1])
            for p in parents:
                retire(p)
            rest = sorted((p for p in parents if p.centre.key != c.key), key=lambda p: -p.best)
            if rest:
                st["runner_up"] = rest[0].centre
            d = birth(c, s, f"merge {mj.id}, {mj.why}: best from {src}, {mj.back}/{mj.sent} children back", parents)
            for g in list(deferred):                         # points that arrived on the merging hill
                if hills.same(g[0], d.centre):
                    deferred.remove(g)
                    if d.offer(g[0], g[1], clock(), fp_of.get(g[0].key), core_of.get(g[0].key)):
                        say("domain %d -> %.4f, a point that arrived during merge %d t+%.1fm", d.id, g[1], mj.id,
                            elapsed() / 60)
            return d

        def merge_find(c: Config, s: float) -> bool:
            """A strong off-hill find with no slot of its own: merged into the nearest live domain (pairing plus
            mutated children; the best of the find, the centre and the children becomes the centre). Only a find
            within merges.close of that domain's best: a weaker one would mostly spend children on its own hill."""
            lv = live()
            if not lv:
                return False
            d = min(lv, key=lambda x: hills.distance(x.centre, c))
            find_done.add(c.key)
            if s < d.best - tune("merges.close"):
                return False
            if tune("merges.on_evidence"):
                start_merge([d], f"a scout find at {s:.4f} merged in", guest=(c, s))
                return True
            find_merge(d, c, s)
            return True

        def find_merge(d: Climber, c: Config, s: float) -> None:
            """One off-hill find per domain, recombined with its centre. The children wait behind that domain's scan
            and are sent by fill() at climb priority. A better find replaces the one already queued."""
            existing = [fm for fm in find_merges.values() if fm.domain == d.id]
            if existing and s <= existing[0].find[1]:
                return
            if existing:
                cancel_behind(d.id, ("findmerge",))
            if d.offer(c, s, clock(), fp_of.get(c.key), core_of.get(c.key)):
                say("domain %d -> %.4f, an off-hill find merged in t+%.1fm", d.id, s, elapsed() / 60)
            fm = FindMerge(ids[1], d.id, (c, s))
            ids[1] += 1
            kids = [k for k in genetic.merge_children(space, d.centre, c, rng, names())
                    if k.key not in d.tried and fresh(k)]
            if len(kids) > tune("merges.find_children"):
                kids = rng.sample(kids, tune("merges.find_children"))
            if kids:
                find_merges[fm.id] = fm
                behind.setdefault(d.id, []).extend(("findmerge", k, fm.id) for k in kids)
            say("find merge %d: an off-hill find at %.4f into domain %d (%.4f), %d children behind its scan "
                "t+%.1fm", fm.id, s, d.id, d.best, len(kids), elapsed() / 60)

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

        def test_find(c: Config, s: float, mode: str = "found") -> None:
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
                find_test[mid.key] = (c, s, d.id, mode)

        def judge_find(key: str, mid_score: Optional[float]) -> None:
            item = find_test.pop(key, None)
            if item is None:
                return
            c, s, did, mode = item
            d = domains.get(did)
            ref = d.best if d is not None else s
            if mid_score is None or mid_score < min(s, ref) - tune("climb.accept_min"):          # a valley: a hill of its own
                valley = f"valley at {mid_score:.4f}" if mid_score is not None else "the midpoint failed"
                count = len(live()) + len(merges)
                if count < max(n_domains if window() else 0, minimum()):
                    birth(c, s, f"a scout find, {valley}")
                    return
                lv = [x for x in live() if x.state != "final"]
                weak = min(lv, key=lambda x: (x.progress(clock()), x.best), default=None)
                if mode == "replace" and weak is not None and s >= weak.best + tune("domains.replace_margin"):
                    retire(weak)
                    say("domain %d (%.4f) retired: a find at %.4f replaces it t+%.1fm", weak.id, weak.best, s,
                        elapsed() / 60)
                    birth(c, s, f"a scout find, {valley}, replacing domain {weak.id}")
                    return
                merge_find(c, s)                               # no slot: into the nearest domain
            elif d is not None and status[d.id] == "live" and d.offer(c, s, clock(), fp_of.get(c.key), core_of.get(c.key)):
                say("domain %d -> %.4f, a scout find on its hill (midpoint %.4f) t+%.1fm", d.id, s,
                         mid_score, elapsed() / 60)

        def judge_midpoint(key: str, score: Optional[float]) -> None:
            pair = midpoint.pop(key, None)
            if pair is None:
                return
            a, b = domains.get(pair[0]), domains.get(pair[1])
            if a is None or b is None or status[a.id] != "live" or status[b.id] != "live":
                return
            if score is not None and score >= min(a.best, b.best) - tune("climb.accept_min"):
                start_merge([a, b], f"one hill by midpoint {score:.4f}")

        # ---- the controller ------------------------------------------------------------------------------------------
        def control() -> None:
            now = clock()
            # merges that are back
            for mj in list(merges.values()):
                if mj.back >= tune("merges.resolve_back") * mj.sent:
                    resolve(mj)
            # frozen domains: merge into a live domain on the same hill, else retire (their best stays a candidate)
            for d in [x for i, x in domains.items() if status[i] == "live" and x.state == "exhausted"]:
                others = [x for x in live() if x is not d]
                near = min(others, key=lambda x: hills.distance(x.centre, d.centre), default=None)
                if tune("merges.on_evidence") and near is not None and hills.same(near.centre, d.centre):
                    start_merge([near, d], "a frozen domain on a live one's hill")
                else:
                    retire(d)
                    say("domain %d retired at %.4f (its climb is exhausted) t+%.1fm", d.id, d.best, elapsed() / 60)
            # evidence: two domains on one hill
            lv = [d for d in live() if d.state != "final"]
            busy: Set[int] = set()
            for dist, a, b in same_hill_pairs(lv) if tune("merges.on_evidence") else ():
                if a.id in busy or b.id in busy:
                    continue
                if dist < MIN_FAR:
                    start_merge([a, b], f"one hill, {dist} settings apart")
                    busy |= {a.id, b.id}
                elif abs(a.best - b.best) <= tune("merges.close") and dist < MIN_FAR + 2 and \
                        now - tested_pair.get((a.id, b.id), -1e9) > 60 * tune("merges.midpoint_retest_min"):
                    test_midpoint(a, b)
            # The climb ends at the final. In the last minutes a new centre cannot finish a scan, so a stall
            # frees its workers for the basins already climbing and no replacement is born.
            t_now = left()
            until_final = (t_now - at("final")) if t_now is not None else 1e9
            until_by = at("reduce_to_1_by") - elapsed()
            no_new = min(until_final, until_by) <= T(60 * tune("schedule.no_new_before_final_min"))
            # Halfway from founding to the final, once every first scan is back, the twelve that are highest
            # after adding their gain since that scan keep the workers. The others retire.
            if not st["final"] and st["founded"] and not st["cut"]:
                climbing = [d for d in live() if d.state != "final"]
                final_at = (tl if tl is not None and tl > 0 else ROUND_S) - at("final")
                half = st["founded_at"] + 0.5 * max(0.0, final_at - st["founded_at"])
                if len(climbing) > 12 and elapsed() >= half and all(d.first_scan_score is not None for d in climbing):
                    ranked_d = sorted(climbing, key=lambda d: (2 * d.best - d.first_scan_score, d.best), reverse=True)
                    keep, drop = ranked_d[:12], ranked_d[12:]
                    st["cut"] = True
                    for d in keep:
                        d.history.append((clock(), d.best, "cut"))   # the stall clock starts over for the basins that stay
                    for d in drop:
                        gain = d.best - d.first_scan_score
                        handed = d.hand(d.gained(), keep)
                        say("domain %d retired at %.4f (cut, gain %.4f); %d gained settings go to the basins that stay t+%.1fm",
                            d.id, d.best, gain, handed, elapsed() / 60)
                        retire(d)
                    say("cut at t+%.1fm: keeping %s", elapsed() / 60,
                        " ".join("%d:%.3f" % (d.id, d.best) for d in keep))
            # A failed first scan loses its seat. The next held candidate takes it, one at a time: the one after
            # waits until this scan is back. Sooner than minute 34 when every live scan is already back.
            if not st["final"] and st["founded"] and reserve and not st["cut"]:
                lv = [d for d in live() if d.state != "final"]
                swap_out = any(origin.get(d.id) == "the next candidate" and d.first_scan_score is None for d in lv)
                scans_back = bool(lv) and all(d.first_scan_score is not None for d in lv)
                if not swap_out and (scans_back or elapsed() >= at("reduce_to_3_start")):
                    decisive = tune("climb.decisive")
                    flat = [d for d in lv if d.first_scan_score is not None
                            and d.first_scan_score < d.birth_score + decisive]
                    for d in flat:
                        say("domain %d retired at %.4f (first scan %.4f, opened %.4f) t+%.1fm", d.id, d.best,
                            d.first_scan_score, d.birth_score, elapsed() / 60)
                        retire(d)
                        st["owed"] += 1
                    if no_new or st["cut"]:
                        st["owed"] = 0
                    elif st["owed"] and reserve and not any(
                            origin.get(d.id) == "the next candidate" and d.first_scan_score is None for d in live()):
                        c, s = reserve.pop(0)
                        birth(c, s, "the next candidate")
                        st["owed"] -= 1
                        say("next candidate at %.4f, %d held, %d swaps still waiting t+%.1fm", s, len(reserve),
                            st["owed"], elapsed() / 60)
            # A basin that has stopped moving gives its workers to the next candidate, or, in the last minutes of
            # the climb, to the basins already climbing. The clock starts at the last
            # time its best changed, and only after the first scan has come back. The highest live basin
            # keeps its workers through a stall; the final still includes it.
            if not st["final"] and st["founded"]:
                limit = T(60 * tune("domains.stall_min"))
                lv = live()
                lead = max((d.best for d in lv), default=None)
                stalled = [d for d in lv if d.state != "final" and d.first_scan_score is not None
                           and not d.scan_open() and now - d.history[-1][0] >= limit
                           and lead is not None and d.best < lead]
                for d in stalled:
                    if not st["cut"]:
                        continue                               # the cut ranks every basin; a stall waits for it
                    if (no_new or st["cut"]) and not any(x is not d for x in live()):
                        continue                               # the last basin keeps its workers
                    say("domain %d retired at %.4f (no gain for %.0f min) t+%.1fm", d.id, d.best,
                        tune("domains.stall_min"), elapsed() / 60)
                    retire(d)
                    if no_new or st["cut"]:
                        say("its workers join the basins still climbing t+%.1fm", elapsed() / 60)
                        continue
                    nxt = reserve.pop(0) if reserve else None
                    if nxt is None:
                        ranked = [(c, s) for c, s in sorted(candidates.values(), key=lambda t: -t[1])
                                  if c.key not in was_centre]
                        avoid = [x.centre for x in domains.values()]
                        picked = hills.take(ranked, 1, tune("domains.found_min_far"), top() - 5.0, avoid)
                        nxt = picked[0] if picked else None
                    if nxt is None:
                        say("domain %d's seat stays empty (no candidate far enough) t+%.1fm", d.id, elapsed() / 60)
                        continue
                    birth(nxt[0], nxt[1], "the next candidate")
                    say("next candidate at %.4f takes a stalled seat, %d held t+%.1fm", nxt[1], len(reserve),
                        elapsed() / 60)
            # finds: off-hill configs, best first, recombined into the nearest domain. A new basin is founded
            # after the opening only when a flat one is replaced.
            if st["founded"] and live():
                ranked = sorted(candidates.values(), key=lambda t: -t[1])
                cs = centres(standing())
                weakest = min((d.best for d in live()), default=float("-inf"))
                for c, s in ranked:
                    if s < max(top() - tune("domains.worth_a_hill"), weakest):
                        break
                    if c.key in was_centre or c.key in find_done or hills.owner(c, cs) is not None:
                        continue
                    if merge_find(c, s):
                        cs = centres(standing())
                        if not live():
                            break
            # sharing between domains: their strongest settings, offered as moves
            if now - st["last_share"] >= 60 * tune("sharing.every_min") and len(live()) > 1:
                st["last_share"] = now
                for d in live():
                    pairs = [(a, v) for a, v in d.pros.items()]
                    for a in sorted(d.axis_gain, key=lambda x: -d.axis_gain[x])[:tune("sharing.settings")]:
                        if d.axis_gain[a] <= 0 or a in d.pros:
                            continue
                        v = d.centre.params.get(a)
                        if v is not None:
                            pairs.append((a, v))
                    for a, v in pairs:
                        for e in live():
                            if e is not d and not space._same(space.specs[a], e.centre.params.get(a), v):
                                e.queue("share", {a: v})
            # the final
            t = left()
            if not st["final"] and st["founded"] and (t is not None and t < at("final") or elapsed() >= at("reduce_to_1_by")):
                final_round()
            if st["final"] and not st["picked"] and t is not None and t < tune("schedule.pick_s"):
                pick()

        def final_round() -> None:
            """One merge, at the start of the final window: the strongest basin keeps the seat and takes the
            others' differing settings, plus one capped set of children. The halved-step search starts after."""
            st["final"] = True
            for mj in list(merges.values()):
                resolve(mj)
            lv = [d for d in live() if d.state != "final"]
            if len(lv) > 1:
                say_basins("final merge")
                ordered = sorted(lv, key=lambda d: -d.best)
                retire_merge(ordered, "final", None)
                win = ordered[0]
            else:
                sd = standing()
                if not sd:
                    return
                win = max(sd, key=lambda d: d.best)
                others = sorted((d for d in sd if d is not win), key=lambda d: -d.best)
                for d in others:
                    retire(d)
                if st["runner_up"] is None and others:
                    st["runner_up"] = others[0].centre
            status[win.id] = "live"
            win.start_final()
            say("final at t+%.1fm on domain %d (%.4f)", elapsed() / 60, win.id, win.best)

        def pick() -> None:
            """The highest local score. Main submits client.pick_key."""
            st["picked"] = True
            ranked = sorted(candidates.values(), key=lambda t: -t[1])
            if not ranked:
                return
            chosen, best = ranked[0]
            client.pick_key = chosen.key
            picked.update(key=chosen.key, score=best, best=best)
            say("pick %s (%.4f) t+%.1fm", chosen.key, best, elapsed() / 60)
            say_basins("pick")

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
                order = sorted(range(len(pool)), key=lambda i: -(mean[i] + tune("scout.kappa") * spread[i]))
                pool = [pool[i] for i in order]
            return pool[:n]

        def explore(idle: Optional[int]) -> bool:
            """The scout sampler. Off unless scout.share is set above 0."""
            return send_scout(idle)

        def send_scout(idle: Optional[int]) -> bool:
            if not scout_queue:
                scout_queue.extend(scout_pick(max(1, min(64, idle or 64))))
            while scout_queue:
                c = scout_queue.pop(0)
                if fresh(c) and send_global(c, "scout", tune("riders.global"), 1):
                    st["scouts"] += 1
                    return True
            return False

        # ---- allocation ------------------------------------------------------------------------------------------------
        def fill() -> int:
            n = 0
            while not client.exhausted:
                idle = client.idle_slots()
                if idle is not None and idle <= 0:
                    break
                if idle is None and len(pending) >= cap:
                    break
                if not st["founded"]:
                    break                                      # the opening's waves are queued by opening_step
                share_now = 0.0 if focus_axes else (tune("scout.share") if window() else tune("scout.late_share"))
                if share_now > 0 and not focus_axes and st["scouts"] < int(cap * share_now) and explore(idle):
                    n += 1
                    continue
                lv = live()
                if lv:
                    share = cap * (1 - share_now)
                    targets = {d.id: share / len(lv) for d in lv}   # even split until the final merge leaves one
                    alloc.clear(); alloc.update(targets)
                    waiting = [d for d in lv if behind.get(d.id) and not d.scan_open()]
                    for d in lv:
                        d.eager_wide = len(lv) <= tune("climb.eager_wide_domains")
                        if d in waiting:
                            continue                       # the merge's children go out before the next sweep
                        if not d.main_pending:                 # the main line refills whenever it is ready; the
                            d.refill(names(), rng, max(0, int(targets[d.id] - d.inflight) + 1))   # share caps chains
                    if waiting:
                        d = max(waiting, key=lambda x: targets[x.id] - x.inflight)
                        q = behind.get(d.id)
                        if q:
                            role, cfg, mid = q.pop(0)
                            if not q:
                                behind.pop(d.id, None)
                            if fresh(cfg) and propose(Unit(role, cfg, domain=d.id, merge=mid or None), [], 3):
                                fm = find_merges.get(mid)
                                if fm is not None:
                                    fm.sent += 1
                                n += 1
                        continue
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
                    for d in sorted(lv, key=lambda x: targets[x.id] - x.inflight, reverse=True):
                        if d.scan_open():
                            continue                       # the scan goes out before any recombination child
                        q = behind.get(d.id)
                        if not q:
                            continue
                        role, cfg, mid = q.pop(0)
                        if not q:
                            behind.pop(d.id, None)
                        if not fresh(cfg):
                            continue
                        if propose(Unit(role, cfg, domain=d.id, merge=mid or None), [], 3):
                            fm = find_merges.get(mid)
                            if fm is not None:
                                fm.sent += 1
                            n += 1
                            break
                    else:
                        if share_now > 0 and not focus_axes and explore(idle):
                            n += 1
                            continue
                        break
                    continue
                if share_now > 0 and not focus_axes and explore(idle):
                    n += 1
                    continue
                break
            return n

        # ---- results -------------------------------------------------------------------------------------------------
        def note(unit: Unit, is_head: bool, o: Outcome) -> None:
            now = clock()
            score = o.score if o.ok else None
            failed = (not o.ok) or score is None or score <= 0
            fid = o.result.fidelity if o.result is not None else None
            cs = client.callsets.get((o.config.key, fid)) if fid else None
            learned = reductions.note(
                o.config, None if failed else score,
                None if cs is None else cs.haplotype_cap)
            if learned:
                say("reductions: %s t+%.1fm", learned, elapsed() / 60)
            rr = o.result
            fp = (rr.region_fp / rr.n_targets) if (rr and not failed and rr.n_targets) else None
            core = rr.core if (rr and not failed and rr.core is not None) else None
            if core is not None:
                core_of[o.config.key] = core
            ranker.add(o.config, core, fp, failed)
            if fp is not None:
                fp_of[o.config.key] = fp
            if not failed and rr is not None and rr.core is not None and fp is not None:
                remember("core", rr.core, o.config)
                remember("fp", fp, o.config)
            if not failed:
                k = o.config.key
                if k not in candidates or score > candidates[k][1]:
                    candidates[k] = (o.config, score)
                if len(candidates) > 1200:
                    for key in sorted(candidates, key=lambda x: -candidates[x][1])[600:]:
                        del candidates[key]
            if unit.role in ("wave", "scout", "support"):
                if is_head:
                    if unit.role == "wave":
                        wave_back[unit.wave] += 1
                    else:
                        st["scouts"] = max(0, st["scouts"] - 1)
                if st["founded"] and not failed:
                    i = owner_of(o.config)
                    if i is not None and domains[i].offer(o.config, score, now, fp, core):
                        say("domain %d -> %.4f, a global point on its hill t+%.1fm", i, score, elapsed() / 60)
                    elif i is None:
                        j = hills.owner(o.config, centres(standing()))
                        if j is not None and status[j] == "merging":
                            deferred.append((o.config, score))   # handed to the merged domain when it is born
                            del deferred[:-200]
                return
            if unit.role == "midpoint":
                judge_midpoint(unit.head.key, None if failed else score)
                return
            if unit.role == "findtest":
                judge_find(unit.head.key, None if failed else score)
                return
            if unit.role in ("findmerge", "recomb"):
                fm = find_merges.get(unit.merge) if unit.merge else None
                if fm is not None:
                    fm.back += 1
                    queued = any(mid == fm.id for _, _, mid in behind.get(fm.domain, ()))
                    if fm.back >= fm.sent and not queued:
                        find_merges.pop(fm.id, None)
                d = domains.get(unit.domain) if unit.domain is not None else None
                if not failed and d is not None and status.get(d.id) == "live" and d.offer(o.config, score, now, fp, core):
                    say("domain %d -> %.4f, a recombination child t+%.1fm", d.id, score, elapsed() / 60)
                    cancel_behind(d.id)                    # the centre moved: the rest were built for the old one
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
            if d.consider(unit.move, is_head, o.config, score, failed, now, fp, core):
                say("domain %d -> %.4f via %s %s t+%.1fm", d.id, d.best, unit.move.kind,
                         ",".join(unit.move.axes), elapsed() / 60)
                cancel_behind(d.id)                        # a new centre: children of the old one are dropped
            fresh = d.take_pros()
            if fresh:
                handed = d.hand(fresh, [e for e in live() if e is not d])
                if handed:
                    say("domain %d keeps %s from a config it did not take; %d basins try it t+%.1fm",
                        d.id, ",".join(a for a, _ in fresh), handed, elapsed() / 60)

        def maybe_fit() -> None:
            if ranker.stale() and (fitting[0] is None or fitting[0].done()):
                fitting[0] = asyncio.create_task(asyncio.to_thread(ranker.fit))

        def status_line() -> None:
            st["last_log"] = clock()
            lv = [d for i, d in domains.items() if status[i] != "retired"]
            log.info("engine: t+%.1fm domains %s | merges %s | find merges %s | scouts %d | candidates %d | pending %d",
                     elapsed() / 60,
                     " ".join(f"{d.id}:{d.best:.3f}/{d.state[:5]}/{d.mode[:1]}/{d.inflight}{'' if status[d.id] == 'live' else '/m'}"
                              for d in lv), " ".join(f"{m.id}:{m.back}/{m.sent}" for m in merges.values()) or "-",
                     " ".join(f"{m.id}>{m.domain}:{m.back}/{m.sent}" for m in find_merges.values()) or "-",
                     st["scouts"], len(candidates), len(pending))

        def basin_rows() -> List[dict]:
            rows = []
            for i, d in sorted(domains.items()):
                rows.append({
                    "id": i, "opened": round(d.birth_score, 4),
                    "first": None if d.first_scan_score is None else round(d.first_scan_score, 4),
                    "first_key": d.first_scan_key, "held": d.first_scan_held,
                    "best": round(d.best, 4), "best_key": d.centre.key, "status": status.get(i),
                })
            return rows

        def say_basins(when: str) -> None:
            parts = []
            for r in basin_rows():
                first = "pending" if r["first"] is None else f"{r['first']:.4f}"
                parts.append("d%d opened %.4f first %s best %.4f %s" % (r["id"], r["opened"], first, r["best"], r["best_key"]))
            say("basins at %s: %s", when, " | ".join(parts) or "none")

        def phase() -> str:
            if not st["founded"]:
                return "opening"
            if st["picked"]:
                return "picked"
            if st["final"]:
                return "final"
            e = elapsed()
            if any(len(m.parents) > 2 for m in merges.values()):
                return "merging to 1"
            if not st["replaced"] and e >= at("reduce_to_3_start"):
                return "replacing"
            return "climbing"

        def publish() -> None:
            """What the dashboard shows (client.notes["engine"]), every TICK_S at most."""
            now = clock()
            st["last_publish"] = now
            e = elapsed()
            for i, s_ in status.items():
                if s_ == "retired" and i not in ended:
                    ended[i] = round(e, 1)
            length = tl if tl is not None and tl > 0 else ROUND_S
            client.notes["engine"] = {
                "version": "v4.22", "phase": phase(), "started_ts": round(time.time() - e, 1), "elapsed_s": round(e, 1),
                "basins": basin_rows(),
                "schedule": {"length_s": round(length, 1), "scale": round(scale_t, 4),
                             "opening_max_s": round(T(60 * tune("opening.max_min")), 1),
                             "reduce_to_3_start_s": round(at("reduce_to_3_start"), 1),
                             "reduce_to_3_by_s": round(at("reduce_to_3_by"), 1),
                             "reduce_to_1_start_s": round(at("reduce_to_1_start"), 1),
                             "reduce_to_1_by_s": round(at("reduce_to_1_by"), 1),
                             "final_by_s": round(min(at("reduce_to_1_by"), length - at("final")), 1),
                             "pick_s": round(length - tune("schedule.pick_s"), 1)},
                "capacity": cap, "max_domains": n_domains,
                "waves": [{"wave": w, "sent": h, "back": b, "best_at_launch": round(x, 4),
                           "started_s": round(wave_started[w], 1)}
                          for w, (h, b, x) in enumerate(zip(wave_heads, wave_back, wave_best))],
                "domains": [{"id": d.id, "status": status[d.id], "state": d.state, "mode": d.mode,
                             "best": round(d.best, 4), "centre": d.centre.key, "inflight": d.inflight,
                             "allotted": round(alloc.get(d.id, 0.0), 1), "progress_per_min": round(d.progress(now), 4),
                             "chains": len(d.chains), "perturb_rounds": d.rounds,
                             "tried": len(d.tried), "born_s": round(d.born - t0, 1),
                             "best_s": round(next((t for t, _s, why in reversed(d.history) if why != "cut"), d.born) - t0, 1),
                             "ended_s": ended.get(d.id),
                             "origin": origin.get(d.id, ""), "parents": lineage.get(d.id, [])}
                            for d in domains.values()],
                "merges": [{"id": m.id, "parents": list(m.parents), "why": m.why, "sent": m.sent, "back": m.back,
                            "started_s": round(m.started - t0, 1),
                            "best": round(m.best[1], 4) if m.best is not None else None} for m in merges.values()],
                "scout": {"in_flight": st["scouts"], "finds_tested": len(find_done), "finds_testing": len(find_test),
                          "share": tune("scout.share")},
                "best": round(top(), 4), "candidates": len(candidates), "pending": len(pending),
                "pick": dict(picked) or None, "events": list(events),
                "caps": caps, "dead_axes": sorted(set(book.names) - set(book.live(clock())))}

        # ---- the loop ------------------------------------------------------------------------------------------------
        sd = seeds()
        cells = [] if focus_axes else [c for c in latin_hypercube(space, rng, max(0, cap - len(sd))) if all(c.key != x.key for x in sd)]
        say("capacity %d -> up to %d domains; schedule x%.2f (a %d-minute round)", cap, n_domains, scale_t,
                 round(ROUND_S * scale_t / 60))
        launch_wave(sd + cells, 0 if focus_axes else tune("riders.wave0"),
                    "focus centre" if focus_axes else "seeds + hypercube")
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
            for d in domains.values():                        # the climbers' notes join the event feed
                while d.events:
                    say("domain %d: %s", d.id, d.events.pop(0))
            if clock() - st["last_log"] > 120:
                status_line()
            if clock() - st["last_publish"] >= TICK_S:
                try:
                    publish()
                except Exception:  # noqa: BLE001 - the dashboard must never end the round
                    log.exception("engine: publishing the status failed")
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
