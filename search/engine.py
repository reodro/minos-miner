"""Search engine v9.28 (search/DESIGN.md): A candidates -> B climb hills -> C one fine-tune.

    Every tag, 2026-10-16: one independent search per islands.workers workers (search/islands.py; 810 workers is
    two islands of 405, each a whole A -> B -> C run on its own share, and the pick is the best island's). Each
    climber also sends centre jumps with its grids: one read-pass or assembly setting moved 2 or more rungs, chosen
    by a bandit over setting x distance that learns from the jumps' results (search/leaps.py); a peak is named only
    after one jump per distance class and direction, then two-setting jumps.

    v9.28: while a grid is walking, a new high in the pile does not replace the fingerprint or the centre.
    The grid writes the setting it found.
    v9.27: an elite probe is sent even when it is within four settings of another probe already in flight.
    Those probes are the search. Dropping them leaves the fleet idle and the same rejection runs again.
    After that window, a worker the planned grid did not fill walks further up the same line.
    v9.26: a candidate does not take a climber onto a centre a peaked climber still holds.
    While the GA owns the fleet, every child copies the best config and changes one to three settings.
    That window does not restart. After it, the grid keeps the fleet and the GA does not breed.
    v9.25: while a GA child or a candidate is still moving a climber, the climbers send nothing.
    The fleet breeds, so those updates reach their limit sooner. Once they have been quiet, the
    grid is one rung past the compass and stays there. It does not widen to five rungs.
    v9.24: a config at least 0.10 from every fingerprint, and above the lowest climbing pile, takes
    that climber as soon as it scores. The pile it leaves becomes a candidate. A candidate does
    not take a climber onto a centre another climber already holds. The fine-tune restarts from
    a higher score when one arrives.
    v9.23: climbers keep their seats until the fine-tune. A peaked seat goes to the best candidate.
    The fine-tune climbs the best config and does not breed. A candidate whose score is higher
    than the lowest climbing pile takes that climber.
    The first half of Step B keeps a one-rung line. The second half spreads five rungs along each side.
    v9.22: six piles. Four climb. Two candidates wait, and the best one takes a seat when a climber reaches
    the top. The highest score within 0.10 of a pile becomes that climber's centre. The fingerprint moves
    only when that centre would stay at least 0.10 from every other pile. A read grid whose settings are
    already dead is skipped, and the next grid is assembly.
    v9.21: the highest score within 0.10 of a pile becomes that pile's reference, and its climber develops from it.
    v9.20: a run this climber sent that comes back above the centre takes the centre. A fingerprint
    whose climber has already reached its peak is not climbed again. Every fingerprint a pile
    stands on is logged, and a fingerprint within 0.10 of that history is not explored.
    A pile breeds from up to 32 parents.
    v9.19 seats one climber on each callset pile. A pile is the ball of radius 0.05 around its best fingerprint.
    Pile fingerprints stay at least 0.10 apart. Half of leftover GA breeds inside the piles. Half mutates in
    search of a fingerprint that is at least 0.10 from every pile.


    Step A  search/candidates.py  wave 0 (seeds, portfolio, hypercube) and one GA generation. The 60 parents
                                  are the best scores that are pairwise at least the seat distance apart. The both /
                                  core / fp groups are built inside those 60. A seat must
                                  be within 0.10 of the best score, then far apart by callset distance. Safety limit
                                  candidates.max_min
    Step B  search/climb.py       the first workers // candidates.workers_per_candidate of that pick climb at once.
                                  Main points go out first; genotyping rides inside each point's run. PairHMM,
                                  cheap pairs, and the rest of the cheap ladders run on the centre only, after
                                  those points are out. The GA is last. A mutation moves one rung. The share of
                                  Step B children that mutate starts at 0.3 and rises linearly to 0.7 by the end
                                  of Step B. When the best
                                  score is shared, the climber on that score probes one locked setting by one rung,
                                  once, and the other free workers keep breeding. A child has to be at least
                                  the seat distance from every climber. A child that differs in fewer than
                                  four settings from another climber, or from another GA child still in flight, is not sent. Step B opens with the whole fleet on
                                  the GA. Grids start once a free GA score beats the lowest climber, or after
                                  climb.ga_open_min. After that a free GA that is on a near climber's hill
                                  goes to the lowest near seat that it beats; a far seat takes only a GA that
                                  is still in the edge family (a higher AGGRESSIVE child that differs in
                                  four genotyping settings is this shelf, not a second hill). A GA child already sitting
                                  on a climber is skipped. It does not wait for a grid to end.
                                  A value that has already tied a config's exact score is not sent again.
                                  An untried value of that setting still is.
                                  Two flat grids on a seat
                                  0.10 or more behind the best seat hand that seat to the next candidate.
                                  Compass then a cube on the axes that beat the centre (v9.17 / B28),
                                  lines, step growth
                                  (long ladders 1, 2, 3, 5), 20 random jumps at a stall. When a hill is found,
                                  that seat climbs the next Step A result. When none of those is far from
                                  the climbers still seated, the best scored config that is far takes the seat.
                                  One of the opening seats is the edge family (validator-legal values the
                                  published API does not propose), even when it is up to 20 behind. Family A
                                  is a measured local opt on this shelf. Scouts snap onto the edge box. Leftover
                                  GA children are half elite-anchored (v9.17 / B26): the current best is one
                                  parent, and the child must differ in at least four searched settings.
                                  Leftover cores take GA children on the same pass as the compass (v9.18 / B29);
                                  the 60 far parents rebuild only when the best score moves.
                                  A far seat only climbs a GA child that is far from every near centre (B30);
                                  a near-hill GA goes to a near seat, so the far hill is not overwritten.
                                  Ends when the fine-tune starts, eight scaled minutes before the pick.
                                  A peaked seat goes to the best candidate through that whole climb.
    Step C  search/far.py refine  one climb of the best config until the pick. No GA generation, and the
                                  fine-tune does not breed. far.pick_s before the end names the highest
                                  local score in client.pick_key

The schedule is written for a 72-minute round and scales with the real one. Notes for main and the review:
client.notes["engine"] (version, phase, basins: one row per climber).
"""
from __future__ import annotations

import asyncio
import json
import logging
import math
import os
import random
import time
from typing import Dict, List, Optional, Sequence, Set, Tuple

from common.params import ParamRanges
from common.tuning import tune
from common.types import FULL, Config
from search import far, islands, portfolio
from search.bamgroups import bam_bounds
from search.base import Strategy
from search.candidates import StepA
from search.client import EvalClient, Outcome
from search.crashguard import CrashGuard
from search.duplicates import PARKED, DuplicateGuard
from search.climb import Climber, occupy
from search.distance import Distance
from search.engine_tag import VERSION
from search.repeats import TiedCallsets
from search.space import Space, bam_caps
from search.steps import Ladders

log = logging.getLogger(__name__)
ROUND_S = 72 * 60                    # the schedule is written for this much search time; it scales with it
TICK_S = 2.0


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


SWEEP_MAX = 60                  # riders per trade-off sweep (genotyping ~1% of a run each, PairHMM ~12%)
SWEEP_EVERY_S = 60.0
_NOT_CRASH = ("cancel", "timeout", "timed out", "lost", "killed", "worker", "shutdown")


def _not_a_crash(error: Optional[str]) -> bool:
    """A failure that says nothing about the config (the run was cancelled, timed out, or its worker went away)."""
    e = (error or "").lower()
    return any(w in e for w in _NOT_CRASH)


class Context:
    """What the steps share: the client, the round's space and ladders, distances, the clock, proposing and waiting."""

    def __init__(self, client: EvalClient, space: Space, rng: random.Random, profile: Optional[str]):
        self.client, self.space, self.rng, self.profile = client, space, rng, profile
        self.ladders = Ladders(space)
        self.dist = Distance(client, space)
        self.clock = client.clock
        self.t0 = self.clock()
        tl = client.time_left()
        self.length = tl if tl is not None and tl > 0 else ROUND_S
        self.scale = self.length / ROUND_S
        self.cap = 8
        self.phase = "starting"
        self.proposed: Set[str] = set()
        self.climbers: List[Climber] = []
        self.events: List[list] = []
        self.step_a: Optional[StepA] = None
        self.c_log: List[dict] = []                  # Step C: one row per pass on a peak
        self.ga_futs: list = []                      # idle-core GA children; lowest priority, not a scheduled wave
        self.ga_keys: Set[str] = set()               # configs the GA proposed; a best among these moves the lowest climber
        self.marks: Dict[str, float] = {}            # seconds into the search each step ended
        self.picked: Optional[dict] = None
        self.started_ts = time.time()
        self.last_publish = -1e18
        self.tied = TiedCallsets(space.names)
        self.repeat_skips = 0
        self.failures: Dict[str, int] = {}             # failed runs by source (wave 0, GA, scout, climber kind)
        self.fail_reasons: Dict[str, int] = {}
        self.guard = CrashGuard(space, self.ladders)  # learns this round's GATK crashes; fresh() asks it
        self.dups = DuplicateGuard(space.names, random.Random(17))   # parks predicted repeats of a callset
        self._top: Optional[float] = None            # best score so far (the duplicate guard's "near the best")
        self._released: list = []                    # parked configs to send now (a broken rule, a near twin)
        self._sweep_next: Optional[Config] = None    # the newest best still to sweep (climb.tradeoff_sweep)
        self._swept: Set[str] = set()
        self._sweep_at = -1e18
        client.on_result.append(self._learn_score)

    # ---- time ------------------------------------------------------------------------------------------------------
    def T(self, seconds: float) -> float:
        return seconds * self.scale

    def elapsed(self) -> float:
        return self.clock() - self.t0

    def pick_at(self) -> float:
        return self.t0 + self.length - float(tune("far.pick_s"))

    def over(self) -> bool:
        t = self.client.time_left()
        return t is not None and t <= 0

    # ---- configs ---------------------------------------------------------------------------------------------------
    def note_failure(self, source: str, reason: str) -> None:
        """One failed run, by source; the reason is kept up to its first ' | '."""
        self.failures[source] = self.failures.get(source, 0) + 1
        why = (reason or "unknown").split(" | ")[0][:90]
        self.fail_reasons[why] = self.fail_reasons.get(why, 0) + 1

    def _learn_score(self, outcome) -> None:
        r = outcome.result
        if r is None and getattr(outcome, "error", None):
            tag = getattr(outcome, "tag", "") or ""
            if not tag.startswith("B"):                    # climbers count their own, by kind (Climber.drop_failed)
                self.note_failure({"A": "step A", "G": "GA", "F": "scout"}.get(tag, tag or "other"), outcome.error)
            if outcome.config is not None and not _not_a_crash(outcome.error):
                self.guard.note(outcome.config.params, True, outcome.config.key)
        elif r is not None and outcome.config is not None:
            self.guard.note(outcome.config.params, False, outcome.config.key)
        score = getattr(r, "score", None)
        if score is None or getattr(r, "gate_passed", True) is False or outcome.config is None:
            return
        self.tied.note(outcome.config.key, outcome.config.params, score)
        if self._top is None or score > self._top:
            self._top = score
            self._sweep_next = outcome.config
        self._released.extend(self.dups.note_score(outcome.config, score, self._top))

    def fresh(self, c: Optional[Config]) -> bool:
        if c is None or c.key in self.proposed or (c.key, FULL) in self.client.latest:
            return False
        if self.guard.blocks(c):
            n = len(self.guard.skipped)
            if n == 1 or (n % 200 == 0 and n != getattr(self, "_guard_said", 0)):
                self._guard_said = n
                self.say("crash guard: %d configs not sent (predicted to crash)", n)
            return False
        return True

    def score(self, c: Config) -> Optional[float]:
        r = self.client.latest.get((c.key, FULL))
        if r is None or r.score is None or r.gate_passed is False:
            return None
        return r.score

    def scored(self):
        """(config, score) of every scored config."""
        for (key, fid), r in self.client.latest.items():
            if fid == FULL and r.score is not None and r.gate_passed is not False:
                c = self.client.configs.get(key)
                if c is not None:
                    yield c, r.score

    def best(self) -> Optional[Tuple[Config, float]]:
        top = None
        for (key, fid), r in self.client.latest.items():
            if fid != FULL or r.score is None or r.gate_passed is False:
                continue
            if top is None or r.score > top[1]:
                c = self.client.configs.get(key)
                if c is not None:
                    top = (c, r.score)
        return top

    def propose(self, batch: Sequence[Config], tag: str, batches: Sequence[str], priority: int = 0) -> List[asyncio.Future]:
        """One GATK run's worth: a head and its riders. A higher priority is taken first.

        A config the duplicate guard predicts to repeat a callset is parked (search/duplicates.py): its future is
        done at once with a PARKED outcome, and the guard sends it later if the prediction breaks or could matter.
        Nothing is parked while the pick is certified."""
        batch, batches = list(batch), list(batches)
        park = not str(self.phase).startswith("C: certify")
        out: List[Optional[asyncio.Future]] = [None] * len(batch)
        keep = []
        for i, c in enumerate(batch):
            self.proposed.add(c.key)
            if park and self.dups.decide(c, tag, self._top):
                fut = asyncio.get_event_loop().create_future()
                fut.set_result(Outcome(c, None, PARKED, tag))
                out[i] = fut
            else:
                keep.append(i)
        if keep:
            labels = [batches[i] for i in keep]
            if keep[0] != 0 and batches:
                labels[0] = batches[0]                    # the first config sent leads the run the head would have
            fs = self.client.propose_many([batch[i] for i in keep], tag=tag, batches=labels, priority=priority)
            for i, f in zip(keep, fs):
                out[i] = f
                self.dups.note_sent(batch[i])
        return out  # type: ignore[return-value]

    def sweep_riders(self, centre: Config) -> List[Config]:
        """The trade-off sweep of one config: every genotyping setting within 4 rungs and at every 4th rung of its
        ladder, every PairHMM setting within 2 rungs; nothing already scored, running or parked."""
        from search.space import GENOTYPING_PARAMS, PAIRHMM_PARAMS
        out: List[Config] = []
        for n in GENOTYPING_PARAMS + PAIRHMM_PARAMS:
            if n not in self.space.names or not self.ladders.movable(n):
                continue
            rungs = self.ladders.rungs[n]
            i = int(round(self.ladders.index(n, centre.params.get(n))))
            near = 4 if n in GENOTYPING_PARAMS else 2
            idx = set(range(i - near, i + near + 1))
            if n in GENOTYPING_PARAMS:
                idx |= set(range(0, len(rungs), 4))
            for j in sorted(idx, key=lambda j: (abs(j - i), j)):
                if 0 <= j < len(rungs) and j != i:
                    c = self.space._assemble({**centre.params, n: rungs[j]})
                    if c is not None and c.key != centre.key and self.fresh(c) and all(c.key != x.key for x in out):
                        out.append(c)
        return out[:SWEEP_MAX]

    def send_sweep(self) -> None:
        """climb.tradeoff_sweep: the round's newest best gets one batch of cheap riders over the settings that trade
        core for false positives (2026-10-17: around every recent pick the calling threshold had been tried at 1 of 74
        values; the 7th place's last step was a two-setting genotyping jump). At most one sweep a minute, Step B and C
        only; a higher result moves a climber the way a GA child does."""
        c = self._sweep_next
        if c is None or not int(tune("climb.tradeoff_sweep")) or self.over():
            return
        if not (str(self.phase).startswith("B") or str(self.phase).startswith("C")) or c.key in self._swept:
            return
        if self.clock() - self._sweep_at < SWEEP_EVERY_S:
            return
        self._sweep_next = None
        self._swept.add(c.key)
        self._sweep_at = self.clock()
        riders = self.sweep_riders(c)
        if not riders:
            return
        fs = self.propose(riders, "S", ["big"] * len(riders), priority=1)
        for r in riders:
            self.ga_keys.add(r.key)
        self.ga_futs.extend(fs)
        self.say("trade-off sweep: %d riders around the best %s", len(riders), c.key[:8])

    def send_released(self) -> None:
        """Parked configs whose prediction broke, or whose twin scored near the best: sent now. A higher one moves a
        climber the way a GA child does (adopt_ga_best)."""
        if not self._released:
            return
        todo, self._released = self._released, []
        send = [c for c, _ in todo if (c.key, FULL) not in self.client.latest]
        if not send:
            return
        fs = self.client.propose_many(send, tag="D", batches=["big"] * len(send), priority=1)
        for c in send:
            self.ga_keys.add(c.key)
            self.dups.note_sent(c)
        self.ga_futs.extend(fs)

    async def wait_share(self, futs: Sequence[asyncio.Future], share: float, until: float) -> None:
        """Until `share` of `futs` are done, `until` passes, or the search's time is over."""
        if not futs:
            return
        need = math.ceil(share * len(futs))
        pending = {f for f in futs if not f.done()}
        while len(futs) - len(pending) < need and pending:
            if self.clock() >= until or self.over():
                break
            _, pending = await asyncio.wait(pending, timeout=TICK_S, return_when=asyncio.FIRST_COMPLETED)
            self.publish()

    def owner(self, c: Config, me: Climber) -> Optional[Climber]:
        """Another live climber already centred on `c`."""
        for o in self.climbers:
            if o is not me and o.state in ("climbing", "peak") and o.centre.key == c.key:
                return o
        return None

    # ---- telling -----------------------------------------------------------------------------------------------------
    def say(self, msg: str, *args) -> None:
        island = getattr(self, "island", None)
        if island is not None:
            msg = f"[i{island.index}] " + msg
        log.info("engine: " + msg, *args)
        try:
            text = msg % args if args else msg
        except (TypeError, ValueError):
            text = msg
        self.events.append([round(time.time(), 1), text])
        del self.events[:-200]

    def name_pick(self) -> None:
        """The highest local score goes in client.pick_key (main submits it), from far.pick_s before the end on."""
        top = self.best()
        if top is not None and top[0].key != self.client.pick_key:
            self.client.pick_key = top[0].key
            self.picked = {"key": top[0].key, "score": round(top[1], 4), "at_s": round(self.elapsed(), 1)}
            self.say("pick %s (%.4f) at t+%.1fm", top[0].key, top[1], self.elapsed() / 60)

    def publish(self, force: bool = False) -> None:
        """client.notes["engine"]: what main and the review read (the basins ledger)."""
        self.send_released()
        self.send_sweep()
        now = self.clock()
        if now >= self.pick_at():                    # a climb still waiting for its wave must not miss the pick
            self.name_pick()
        if not force and now - self.last_publish < TICK_S:
            return
        self.last_publish = now
        try:
            top = self.best()
            self.client.notes["engine"] = {
                "version": VERSION, "phase": self.phase, "elapsed_s": round(self.elapsed(), 1),
                "started_ts": round(self.started_ts, 1), "capacity": self.cap,
                "best": None if top is None else round(top[1], 4),
                "schedule": {"length_s": round(self.length, 1), "scale": round(self.scale, 4),
                             "a_max_s": round(self.T(60 * float(tune("candidates.max_min"))), 1),
                             "gen_s": round(self.T(60 * float(tune("candidates.generation_min"))), 1),
                             "b_end_s": round(max(0.0, self.length - float(tune("far.pick_s")) - self.T(8 * 60)), 1),
                             "ga_open_s": round(self.T(60 * float(tune("climb.ga_open_min"))), 1),
                             "c_reserve_s": round(self.T(8 * 60), 1),
                             "pick_s": round(self.length - float(tune("far.pick_s")), 1),
                             **{f"{k}_done_s": round(v, 1) for k, v in self.marks.items()}},
                "ga": {"sent": len(self.ga_futs), "back": sum(1 for f in self.ga_futs if f.done()),
                       "open_reason": getattr(self, "ga_open_reason", None)},
                "waves": self.step_a.wave_rows() if self.step_a is not None else [],
                "candidates": list(self.step_a.picked) if self.step_a is not None else [],
                "basins": [c.row() for c in self.climbers],
                "fingerprints": _fingerprint_history(self),
                "piles": _pile_rows(self),
                "step_c": list(self.c_log),
                "pick": self.picked,
                "failures": dict(self.failures),
                "head_reserve_ws": round(getattr(self, "reserve_ws", 0.0), 1),
                "leaps": _leaps_summary(self),
                "duplicates": {"parked": self.dups.n_parked, "audits": self.dups.n_audits, "near": self.dups.n_near,
                               "rules_broken": self.dups.n_broken, "released": self.dups.n_released,
                               "still_parked": self.dups.still_parked()},
                "crash_guard": {"crashes": len(self.guard.crashed), "completed": len(self.guard.ok),
                                "not_sent": len(self.guard.skipped), "fitted": self.guard.model is not None,
                                "audits": self.guard.n_audits, "audits_crashed": self.guard.n_audit_crashed,
                                "rules_off": self.guard.n_off},
                "events": list(self.events[-200:]),
            }
        except Exception:  # noqa: BLE001 - the notes must never end the round
            log.exception("engine: publishing the notes failed")


def _leaps_summary(ctx) -> dict:
    from search import leaps
    return leaps.summary(ctx)


def _fingerprint_history(ctx) -> list:
    from search.piles import history_rows
    return history_rows(ctx)


def _pile_rows(ctx) -> list:
    from search.piles import pile_rows
    return pile_rows(ctx)


def searcher(ranges: ParamRanges, workers: int, *, tool: str = "gatk", seed: int = 0) -> Strategy:
    rng = random.Random(seed)

    async def run(client: EvalClient) -> None:
        r = client.round
        bounds = bam_bounds(r.bam, ranges, tool) if r.bam else None
        caps = {} if bounds is not None else (bam_caps(r.bam, r.chrom, r.start, r.end) if r.bam else {})
        profile = portfolio.bam_profile(r.bam, r.chrom, r.start, r.end) if r.bam else None
        # One Space for every island: building a second one widens the shared ranges again (Space._widen).
        space = Space(ranges, tool, caps=caps, bounds=bounds)
        if bounds is not None:
            log.info("engine: bam groups: %s", bounds.report)
        cap = await settled_capacity(client, workers)
        polish = os.environ.get("MINOS_POLISH")
        if polish:                                     # a test-only focused search from given centres (search/polish.py)
            from search import polish as fleet_polish
            await fleet_polish.run(client, space, polish)
            return
        probe = os.environ.get("MINOS_PROBE")
        if probe:                                      # a test-only fleet probe of the jump policies; no search
            from search import probe as fleet_probe
            await fleet_probe.run(client, space, probe)
            return
        k = islands.count(cap)
        if k <= 1:
            await search_one(client, space, profile, rng, cap)
            return
        log.info("engine: %d islands on %d workers (%s); each is a whole search on its own share", k, cap,
                 " + ".join(str(q) for q in islands.quotas(cap, k)))
        await islands.run(client, cap, k, lambda view, island: search_one(
            view, space, profile, random.Random(seed * 1000 + 1 + island.index), view.quota, island))

    return run


def arms_path() -> str:
    from common.settings import SETTINGS
    return os.path.join(SETTINGS.data_dir, "search", "jump_arms.json")


def carry_arms(client, arms) -> None:
    """Keep this round's centre-jump arms for the next round's priors (search/leaps.py). A simulated round
    (an executor with a virtual clock) keeps nothing."""
    if arms is None or not carry_enabled(client):
        return
    from search.leaps import save_arms
    try:
        save_arms(arms, arms_path())
    except Exception:  # noqa: BLE001 - the round's result does not depend on this file
        log.exception("engine: saving the jump arms failed")


def round_arms(client):
    """The arms a round starts from: earlier rounds' counts when a real round can read them, plain priors otherwise."""
    from search.leaps import Arms, load_arms
    if not carry_enabled(client):
        return Arms()
    return load_arms(arms_path())


def carry_enabled(client) -> bool:
    """Arms are carried between real rounds. Not in a simulation, and not with MINOS_ARMS_CARRY=0 (replays of one
    round: a replay must not learn from the previous replay of the same round)."""
    if hasattr(getattr(client, "executor", None), "now"):
        return False
    return os.environ.get("MINOS_ARMS_CARRY", "1") != "0"


async def search_one(client, space: Space, profile: Optional[str], rng: random.Random, cap: int,
                     island: Optional["islands.Island"] = None) -> None:
    """One whole search: Step A, Step B, Step C and the pick. With islands, `client` is the island's view."""
    ctx = Context(client, space, rng, profile)
    ctx.island = island
    if island is not None:
        ctx.climber_base = 100 * island.index
        if getattr(island, "arms", None) is not None:
            ctx._jump_arms = island.arms                   # the islands learn their centre jumps together
    else:
        ctx._jump_arms = round_arms(client)
    ctx.cap = cap
    ctx.say("%s: %d workers, %d settings searched, a %d-minute round (schedule x%.2f)", VERSION, ctx.cap,
            len(space.names), round(ctx.length / 60), ctx.scale)

    # ---- Step A ----------------------------------------------------------------------------------------------------
    # SEARCH_START is a one-replay override: climb that frozen config and skip the GA. The normal path is unchanged.
    ctx.phase = "A: candidates"
    ctx.step_a = StepA(ctx)
    forced = os.environ.get("SEARCH_START")
    known_score: Dict[str, float] = {}
    if forced:
        raw = json.loads(open(forced).read())
        seat = space._assemble(raw["params"])
        starts = [seat] if seat is not None else []
        if seat is not None and seat.key == raw.get("key"):
            known_score[seat.key] = float(raw["score"])
        ctx.say("step A skipped: climbing frozen seat %s opened at %.4f (%s)",
                seat.key[:8] if seat is not None else "?", float(raw.get("score") or 0), raw.get("group") or "")
    else:
        try:
            starts = await ctx.step_a.run()
        except Exception:  # noqa: BLE001 - a bug in a step must not end the round
            log.exception("engine: step A failed")
            top = ctx.best()
            starts = [top[0]] if top else []
        else:
            from search.separate import enabled as separate_enabled, gate
            if separate_enabled() and starts:
                try:
                    starts = await gate(ctx, list(starts))
                except Exception:  # noqa: BLE001
                    log.exception("engine: the valley gate failed")
    if island is not None:
        starts = island.apart(ctx, starts)
    ctx.marks["a"] = ctx.elapsed()
    ctx.publish(force=True)
    if os.environ.get("MINOS_ONLY_WAVE0"):
        ctx.say("generation 1 is inside the peak interval; step B and step C are not run")
        ctx.phase = "picked"
        ctx.name_pick()
        ctx.publish(force=True)
        while not ctx.over() and client.in_flight:
            await asyncio.sleep(0 if hasattr(client.executor, "now") else TICK_S)
            ctx.name_pick()
            ctx.publish()
        return

    # ---- Step B ----------------------------------------------------------------------------------------------------
    # The first seats climb at once. When one finds its hill, it climbs the next Step A result.
    # They keep the seats until the fine-tune. A peaked seat goes to the best candidate.
    ctx.phase = "B: climbing"
    end_b = far.fine_tune_at(ctx)
    if end_b < ctx.clock():
        end_b = ctx.clock()
    ctx.b_span = (ctx.clock(), end_b)
    ready: List[Tuple[Config, float, str]] = []
    for c in starts:
        s = ctx.score(c)
        if s is None:
            s = known_score.get(c.key)
        if s is not None:
            if c.key in known_score:
                origin = "frozen seat"
            else:
                row = next((r for r in (getattr(ctx.step_a, "picked", None) or []) if r.get("key") == c.key), None)
                group = (row or {}).get("group") or ""
                if group.startswith("pile"):
                    origin = group
                else:
                    origin = "far" if group == "far" else "step A"
            ready.append((c, s, origin))
    try:
        await occupy(ctx, ready, end_b)
    except Exception:  # noqa: BLE001 - Step B failing must not skip Step C
        log.exception("engine: step B failed")
    ctx.marks["b"] = ctx.elapsed()
    ctx.say("step B done at t+%.1fm: %s", ctx.elapsed() / 60,
            " | ".join(f"c{d.id} {d.opened:.4f} -> {d.score:.4f} ({d.state}, {d.iteration} it)"
                       for d in ctx.climbers) or "no climber")
    ctx.publish(force=True)

    # ---- Step C: one fine-tune of the best config. The climbers have already used the time before it. ------------
    try:
        await far.refine(ctx)
    except Exception:  # noqa: BLE001
        log.exception("engine: step C failed")
    try:
        from search.local_opt import certify_round
        await certify_round(ctx)
    except Exception:  # noqa: BLE001
        log.exception("engine: 1-rung certify failed")
    ctx.marks["c"] = ctx.elapsed()
    ctx.publish(force=True)

    if island is None:
        carry_arms(client, getattr(ctx, "_jump_arms", None))
    if ctx.failures:
        ctx.say("failed runs by source: %s | reasons: %s", ", ".join(f"{k} {v}" for k, v in sorted(ctx.failures.items())),
                "; ".join(f"{v} x {k}" for k, v in sorted(ctx.fail_reasons.items(), key=lambda kv: -kv[1])[:3]))
    js = _leaps_summary(ctx)
    t2 = js["two_setting"]
    if js["tries"] or t2["tries"]:
        ctx.say("centre jumps: %d sent, %d won (%.0f%%), gain %+.4f in all, %+.4f per jump sent; two-setting: %d sent, "
                "%d won, gain %+.4f", js["tries"], js["wins"], 100.0 * js["wins"] / max(1, js["tries"]), js["gain"],
                js["gain"] / max(1, js["tries"]), t2["tries"], t2["wins"], t2["gain"])
    if ctx.guard.summary():
        ctx.say("%s", ctx.guard.summary())
    if ctx.dups.summary():
        ctx.say("%s", ctx.dups.summary())
    if getattr(ctx, "reserve_ws", 0.0) > 0:
        span = max(1.0, ctx.marks.get("b", ctx.elapsed()) - ctx.marks.get("a", 0.0))
        ctx.say("head reserve kept %.0f worker-minutes free (%.1f workers on average over Step B)",
                ctx.reserve_ws / 60, ctx.reserve_ws / span)

    # ---- the pick --------------------------------------------------------------------------------------------------
    ctx.phase = "picked"
    ctx.name_pick()
    ctx.publish(force=True)
    while not ctx.over() and client.in_flight:      # what is still running may still be submitted by main
        await asyncio.sleep(0 if hasattr(client.executor, "now") else TICK_S)
        ctx.name_pick()
        ctx.publish()
