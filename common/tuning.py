"""Tunable settings of the search, the evaluator and main: defaults here, overrides in a JSON file on main.

The file (`TUNING_FILE`: tuning.json in the repo folder, git-ignored; MINOS_TUNING overrides the path) holds only
overrides, as flat names: {"climb.random_jumps": 30, "climb.end_min": 55}. The dashboard edits it, and a
person may too. Only main turns the file on (`TUNING.use_file()`); tests, simulations and scripts run on defaults.

Reading: `tune("climb.random_jumps")` returns the value in force. The file is re-read when it changes (checked at most once a
second), so a "live" setting takes effect at its next use. A "round" setting is read once at a round's start: a change
made during a round waits in `pending` until the next round begins. A value that fails its check is rejected with a
reason and the previous value stays in force, so nothing written to the file can stop mining.
"""
from __future__ import annotations

import collections
import json
import logging
import math
import os
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

log = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
TUNING_FILE = Path(os.environ.get("MINOS_TUNING") or REPO / "tuning.json")
CHECK_EVERY_S = 1.0
LIVE, ROUND = "live", "round"


@dataclass(frozen=True)
class Knob:
    name: str                       # "<group>.<setting>"
    default: Any                    # int, float or a list of ints
    lo: float
    hi: float
    when: str                       # LIVE: at its next use; ROUND: from the next round's start
    unit: str                       # min | s | share | score | runs | count | x | ""
    help: str
    max_len: int = 8                # list settings: at most this many values

    @property
    def group(self) -> str:
        return self.name.split(".", 1)[0]

    @property
    def kind(self) -> str:
        if isinstance(self.default, list):
            return "list"
        return "int" if isinstance(self.default, int) and not isinstance(self.default, bool) else "float"

    def check(self, value: Any) -> Tuple[bool, Any, str]:
        """(ok, the value in its proper type, why not)."""
        def one(v, as_int: bool):
            if isinstance(v, bool) or not isinstance(v, (int, float)):
                return None, f"{v!r} is not a number"
            if isinstance(v, float) and not math.isfinite(v):
                return None, f"{v!r} is not a finite number"
            if as_int:
                if float(v) != int(v):
                    return None, f"{v!r} is not a whole number"
                v = int(v)
            else:
                v = float(v)
            if not self.lo <= v <= self.hi:
                return None, f"{v} is outside {fmt(self.lo)}..{fmt(self.hi)}"
            return v, ""
        if self.kind == "list":
            if not isinstance(value, list) or not 1 <= len(value) <= self.max_len:
                return False, None, f"needs a list of 1 to {self.max_len} whole numbers"
            out = []
            for v in value:
                x, why = one(v, True)
                if x is None:
                    return False, None, why
                out.append(x)
            return True, out, ""
        x, why = one(value, self.kind == "int")
        return (x is not None), x, why


def fmt(x: float) -> str:
    return str(int(x)) if float(x) == int(x) else str(x)


GROUPS = {
    "candidates": "Step A: candidates (wave 0, the genetic generations, the pick)",
    "climb": "Step B: each candidate's grid climb",
    "islands": "Islands: independent searches sharing one fleet",
    "dups": "Duplicate callsets: predicted repeats are parked, audited and released",
    "far": "Step C: what Step B's grid did not look at, and the pick",
    "evaluator": "Evaluator: batches and routing on the fleet",
    "main": "Main: submissions",
    "answer_set": "Answer set: the targets the local score divides by",
}

KNOBS_LIST: List[Knob] = [
    Knob("candidates.generalist_share", 0.0, 0.0, 1.0, ROUND, "share",
         "Share of wave 0 seeded with the cross-round seed portfolio ($DATA_DIR/portfolio/generalists.json, from "
         "test/bench/portfolio_eval.py), before the portfolio. 0 (default since 2026-10-17): none - on six unseen rounds "
         "it was only +0.04 per config better than our existing seeds, and on 20261016T130600Z its best was -1.15."),
    Knob("candidates.portfolio_share", 0.4, 0.0, 1.0, ROUND, "share",
         "Share of wave 0 taken by portfolio seeds (the round's read-length profile first, never the round being "
         "replayed); the rest is a hypercube."),
    Knob("candidates.riders", 4, 0, 6, ROUND, "count", "Calling-threshold riders on each wave 0 point."),
    Knob("candidates.generations", 1, 0, 6, LIVE, "count",
         "Genetic generations after wave 0. One: on 20261008T152900Z and 20261008T175300Z generation 1 raised the "
         "best by 0.15-0.18, generation 2 raised it by nothing, and generation 3 by 0.037. Later generations are "
         "not worth a fleet wave. Step C still breeds the cross (core x fp) itself."),
    Knob("candidates.generation_min", 12.0, 2.0, 30.0, LIVE, "min",
         "Minutes generation 1 is given after it is launched (scaled). The wave is waited until "
         "candidates.launch_back of it is back. candidates.max_min does not abandon it."),
    Knob("candidates.wave0_settle", 1, 0, 1, ROUND, "count",
         "1: wave 0 ends once half of it is back and the best has stood still for 90 s (at the latest at "
         "candidates.launch_back). 0: wave 0 waits for candidates.launch_back of it."),
    Knob("candidates.launch_back", 0.8, 0.1, 1.0, LIVE, "share",
         "A generation launches once this share of the previous wave is back."),
    Knob("candidates.mutate_frac", 0.3, 0.0, 1.0, LIVE, "share",
         "Share of the children that get candidates.mutate_axes settings moved by one rung."),
    Knob("candidates.mutate_axes", 2, 1, 10, LIVE, "count", "Settings a mutated child moves by one rung."),
    Knob("candidates.elite_frac", 0.5, 0.0, 1.0, LIVE, "share",
         "Step B leftover GA: this share of children use the current best as one parent (v9.17). "
         "The rest stay within-group pairs. 0 turns elite-anchored breeding off."),
    Knob("candidates.group_size", 60, 2, 200, LIVE, "count", "Parents per group (both, core, fp), each far from the others."),
    Knob("candidates.core_band", 0.01, 0.0, 0.2, LIVE, "score",
         "The core and fp groups take configs whose core is within this of the best core."),
    Knob("candidates.max_min", 20.0, 2.0, 50.0, LIVE, "min",
         "Safety limit of wave 0 only (scaled). On a heavy region that wave took about 10 minutes to come back. "
         "Generation 1 is separate: candidates.generation_min, and it is waited out."),
    Knob("candidates.workers_per_candidate", 243, 8, 5000, ROUND, "count",
         "One climbing candidate per this many workers (at least one). 243 seats four on a 972-worker fleet."),
    Knob("candidates.min_distance", 0.05, 0.0, 1.0, LIVE, "x",
         "Two candidates are at least this far apart by callset distance (class-weighted 1 - Jaccard)."),
    Knob("candidates.min_setting_distance", 0.05, 0.0, 1.0, LIVE, "x",
         "The same when a callset has no sketch: mean normalised setting difference."),
    Knob("climb.group_size", 5, 1, 7, LIVE, "count",
         "Settings in each grid (every combination at -d, 0, +d: 3^n points)."),
    Knob("climb.accept_min", 0.005, 0.0, 0.1, LIVE, "score",
         "An iteration gaining less than this has failed; both grid kinds failing in a row is a stall."),
    Knob("climb.random_jumps", 20, 0, 500, LIVE, "count",
         "At a stall: random local jumps (2-4 settings, 2-5 rungs each); one that gains climb.accept_min resumes the grid."),
    Knob("climb.leap_share", 0.25, 0.0, 1.0, LIVE, "share",
         "Starting share of each climber's workers that send centre jumps with every grid (one read-pass or assembly "
         "setting moved 2+ rungs, chosen by a bandit). It doubles after three waves with a jump win and halves after "
         "two without (0.1-0.5). Before a peak: one jump per distance class and direction, then two-setting jumps. "
         "0: none."),
    Knob("climb.own_hill_adopt", 1, 0, 1, LIVE, "count",
         "1: a higher GA config on a climber's own hill becomes that climber's centre. 0: it is left there (before "
         "2026-10-16). For a fleet A/B: it made the fastest early climb, and maybe a lower final peak."),
    Knob("dups.enabled", 1, 0, 1, LIVE, "count",
         "1: a config predicted to repeat a callset already sent is parked, audited, and sent the moment the "
         "prediction breaks or its twin scores near the best (search/duplicates.py). 0: every config is sent."),
    Knob("dups.min_ties", 5, 1, 100, LIVE, "count",
         "Ties (and no different score) before a one-setting change is predicted to repeat a callset. 3 lost two "
         "better configs in the v4 replay; 5 lost none."),
    Knob("dups.audit", 0.1, 0.0, 1.0, LIVE, "share",
         "Share of predicted repeats sent anyway, so a wrong rule breaks and releases what it parked."),
    Knob("dups.near", 0.2, 0.0, 5.0, LIVE, "score",
         "A predicted repeat whose twin scores within this of the island's best is sent anyway: only those could "
         "move the round's best."),
    Knob("climb.tradeoff_sweep", 0, 0, 1, LIVE, "count",
         "1: the round's newest best (Step B and C) gets one batch of cheap riders over the genotyping and PairHMM "
         "settings (search/engine.py Context.send_sweep), at most once a minute. 0: off."),
    Knob("climb.model_filter", 1, 1, 10, LIVE, "count",
         "k > 1: the GA breeds k times as many children and sends the 1/k a per-class score model of the round ranks "
         "highest (search/model_filter.py). 1: every child is sent (before 2026-10-17)."),
    Knob("climb.stall_waves", 2, 0, 20, LIVE, "count",
         "Finished grids in a row, each below climb.accept_min, after which a climber certifies its peak (every "
         "untried one-rung neighbour and its centre jumps in one wave). 0: only when every group is flat (before "
         "2026-10-16)."),
    Knob("climb.all_steps", 0, 0, 1, LIVE, "count",
         "1: every grid wave sends each setting at every step of its plan at once (+-1 and the big steps of "
         "step_plan, e.g. 1,3,5,7,9 rungs), not one step size per wave; the cube keeps the scheduled step. 0: one step."),
    Knob("islands.by_climber", 0, 0, 1, LIVE, "count",
         "1: the fleet is shared equally by the climbing climbers, whichever island they are on: an island's quota is "
         "its climbers' share (one climber's worth with none). 0: equal islands."),
    Knob("climb.heads_first", 0, 0, 1, LIVE, "count",
         "1: a grid wave's heads go out even past climb.grid_cores (they decide when the wave ends); extra grid points "
         "stay within the quota."),
    Knob("climb.flat_shift", 0.0, 0.0, 1.0, LIVE, "share",
         "With the core quotas: a climber whose last grid wave gained nothing moves this share of its grid runs to "
         "centre jumps until it moves again."),
    Knob("climb.grid_cores", 0.0, 0.0, 1.0, LIVE, "share",
         "Share of an island's workers the grid may keep busy in Step B (heads, line points, packs, more rungs). With "
         "climb.jump_cores both > 0, they replace climb.grid_fill; workers the two leave free go to RF-UCB and the GA."),
    Knob("climb.jump_cores", 0.0, 0.0, 1.0, LIVE, "share",
         "Share of an island's workers centre jumps may keep busy in Step B (see climb.grid_cores)."),
    Knob("climb.grid_fill", 0.0, 0.0, 1.0, LIVE, "share",
         "Share of the free workers left after the climbers' heads, packs and centre jumps that go to more grid points "
         "(more rungs of the wave's settings, Climber.extend) before RF-UCB and the GA. 0: the GA takes them all."),
    Knob("climb.rf_ucb_share", 0.0, 0.0, 1.0, LIVE, "share",
         "Share of the free workers left after the grids that go to RF-UCB proposals (search/rf_ucb.py step_b); the GA "
         "breeds on the rest. 0: no RF-UCB in Step B."),
    Knob("climb.rf_ucb_kappa", 0.0, 0.0, 5.0, LIVE, "share",
         "RF-UCB ranks by predicted score + kappa x the spread across trees. 2026-10-19 offline check on four "
         "full-fleet rounds (test/bench/rf_ucb_offline.py): improvers among the top 100 summed 146 at 0, 90 at 0.25, "
         "49 at 0.5, 10 at 1."),
    Knob("climb.head_reserve", 24, 0, 200, LIVE, "count",
         "Workers the GA leaves free for a climber's next heads once its wave is half back (at most this many per "
         "climber). 0: the GA fills every free worker."),
    Knob("islands.workers", 324, 50, 100000, ROUND, "count",
         "One independent search per this many workers (810: two islands of 405). Below twice this, one search."),
    Knob("islands.min_settings", 8, 0, 30, LIVE, "count",
         "A Step A seat that differs in fewer settings than this from another island's seat is passed over "
         "while the island has another."),
    Knob("climb.wait_riders", 0, 0, 1, LIVE, "count",
         "1: a wave decides once its riders are back too. 0: once its grid points are back; riders that come back "
         "later count at the next decision."),
    Knob("climb.quorum", 0.8, 0.5, 1.0, LIVE, "share",
         "A wave counts as back at this share: the slowest fifth no longer holds the next wave (their scores still "
         "count when they return). 1.0 before 2026-10-16."),
    Knob("climb.move_at_once", 1, 0, 1, LIVE, "count",
         "1: a run that beats the centre by climb.accept_min moves it at once (the wave closes early; no cube on the "
         "old centre once a compass point clears the bar). 0: the climber moves only when its whole wave is back."),
    Knob("climb.end_min", 55.0, 10.0, 72.0, LIVE, "min",
         "No longer the cut. Step B runs until the fine-tune, eight scaled minutes before the pick. "
         "Kept so an older override still loads."),
    Knob("climb.ga_quiet_min", 4.0, 1.0, 20.0, LIVE, "min",
         "After a GA child or a candidate moves a climber, the climbers send nothing for this many "
         "minutes (scaled) and the fleet breeds. Then one rung past the compass walks to the peak."),
    Knob("climb.after_peak", 1, 0, 1, ROUND, "count",
         "1: a climber that named a 1-rung peak keeps sending centre jumps (no grid) instead of stopping; after three "
         "jump waves without a win (and 10+ minutes left) it kicks 4-6 settings by 3-9 rungs and climbs again. 0: it stops."),
    Knob("climb.after_peak_jumps", 24, 4, 200, ROUND, "count", "Centre jumps per wave of a peaked climber."),
    Knob("climb.kick_configs", 12, 2, 100, ROUND, "count", "Configs scored for one kick; the climber moves to the best."),
    Knob("climb.big_steps", 1, 0, 1, ROUND, "count",
         "1: each grid setting starts at a step that fits its ladder (<=40 rungs: 1; 41-150: 3; >150: 3), widens "
         "after an iteration that did not move it (to 5 / 9 / climb.big_step_large_max) and then narrows back to 1. "
         "0: every setting steps one rung, as before 2026-10-18."),
    Knob("climb.big_step_large_max", 27, 9, 100, ROUND, "rungs",
         "The widest grid step of a setting with more than 150 rungs (padding, assembly region sizes)."),
    Knob("climb.ga_open_min", 0.0, 0.0, 20.0, ROUND, "min",
         "At the start of Step B the whole fleet breeds before any grid. 0 (the default since 2026-10-18): no such "
         "step, the grids start at once. Over 33 openings of 2026-10-17/18 it gained +0.062 in 2.4 min on average "
         "against +0.073 for the grids in the same time after it, and nothing at all in 16 of them."),
    Knob("far.coarse", 25, 3, 500, LIVE, "count",
         "Step C: rungs tried beyond the nearest tried value on each side of a peak, at least (more when the fleet has room)."),
    Knob("far.second_margin", 0.1, 0.0, 2.0, LIVE, "score",
         "Step C: a second high point on a setting's profile within this of the peak is refined and climbed."),
    Knob("far.landings", 3, 0, 20, LIVE, "count", "Step C: second high points per peak that get a short climb."),
    Knob("far.landing_iterations", 2, 1, 10, LIVE, "count", "Step C: grids in a short climb."),
    Knob("far.pick_s", 90.0, 10.0, 600.0, LIVE, "s",
         "The highest-scoring config is named this many seconds before the search ends (not scaled)."),
    Knob("evaluator.small_batch", 1.3, 1.0, 3.0, LIVE, "runs",
         "Most work in a small batch, in solo runs: a waited-for result and its cheap riders."),
    Knob("evaluator.big_batch", 2.3, 1.0, 10.0, LIVE, "runs", "Most work in a big batch, in solo runs."),
    Knob("evaluator.overdue", 1.5, 1.1, 5.0, LIVE, "x",
         "A small batch running this many times its expected time is duplicated on a free fast node."),
    Knob("evaluator.slow_factor", 1.5, 1.0, 5.0, LIVE, "x",
         "A node slower than this times the fastest is slow: small batches prefer the others."),
    Knob("evaluator.sick_factor", 3.0, 1.5, 20.0, LIVE, "x",
         "A node slower than this times the fleet's median gets no small batches at all, only big ones."),
    Knob("evaluator.keep_rounds", 3, 1, 50, LIVE, "count",
         "Rounds the fleet keeps on main and the nodes. Older ones are pruned once none of their jobs is queued or "
         "running. Each kept round holds its BAM and callsets on every node: raise it for multi-round tuning."),
    Knob("evaluator.max_copy_share", 0.02, 0.0, 1.0, LIVE, "share",
         "Most of the ready workers that slow-batch copies may hold at once (0 turns copies off)."),
    Knob("evaluator.split_k", 1, 1, 16, ROUND, "count",
         "Region pieces per full-region GATK run (evaluator/split.py): each config runs as this many padded pieces at "
         "once on one node and is merged, so it comes back several times sooner and each node offers workers / k "
         "slots. 1: off. Read when the node sessions start (a miner restart)."),
    Knob("evaluator.split_pad", 2000, 0, 20000, ROUND, "bp",
         "Padding of each region piece beyond its own part, on both sides (evaluator.split_k > 1)."),
    Knob("evaluator.lost_factor", 6.0, 2.0, 100.0, LIVE, "x",
         "A GATK run outstanding this many times its expected time (and at least 40 minutes) is lost: its jobs are "
         "queued again and the node is told to drop it."),
    Knob("main.resubmit_interval_s", 300.0, 0.0, 3600.0, ROUND, "s", "Least time between two submissions."),
    Knob("main.resubmit_margin", 0.0, 0.0, 5.0, ROUND, "score", "A new best must gain this much to be resubmitted."),
    Knob("main.stop_before_deadline_s", 60.0, 10.0, 900.0, ROUND, "s", "The search stops this long before the deadline."),
    Knob("answer_set.n", 548, 50, 2000, ROUND, "count",
         "How many targets the next round's answer set writes. The local score divides false calls by this. "
         "A change made during a round waits for the next round. Used when main.n_recent_rounds gives no N."),
    Knob("main.n_recent_rounds", 10, 0, 100, ROUND, "count",
         "N for a round the platform gives no N for: the median target count of this many recent rounds whose truth "
         "was out when the round started (main/n_estimate.py). 0: the fixed answer_set.n (before 2026-10-18)."),
]
KNOBS: Dict[str, Knob] = {k.name: k for k in KNOBS_LIST}
assert len(KNOBS) == len(KNOBS_LIST) and all(k.group in GROUPS for k in KNOBS_LIST)
DEFAULTS: Dict[str, Any] = {k.name: (list(k.default) if k.kind == "list" else k.default) for k in KNOBS_LIST}

# Checks across settings: (the settings involved, test on the would-be values, why it fails).
CROSS: List[Tuple[Sequence[str], Callable[[Dict[str, Any]], bool], str]] = [
    (("evaluator.small_batch", "evaluator.big_batch"),
     lambda v: v["evaluator.small_batch"] <= v["evaluator.big_batch"],
     "a small batch cannot be bigger than a big one"),
]


def validate(raw: Dict[str, Any], base: Optional[Dict[str, Any]] = None) -> Tuple[Dict[str, Any], List[dict]]:
    """Check overrides. Returns (the accepted overrides in their proper types, the rejected ones with reasons).
    `base` is what a rejected setting falls back to for the checks across settings (the values in force)."""
    base = dict(base or DEFAULTS)
    ok: Dict[str, Any] = {}
    rejected: List[dict] = []
    for name, value in raw.items():
        k = KNOBS.get(name)
        if k is None:
            rejected.append({"name": name, "value": value, "reason": "unknown setting"})
            continue
        good, v, why = k.check(value)
        if good:
            ok[name] = v
        else:
            rejected.append({"name": name, "value": value, "reason": why})
    would = {n: ok.get(n, DEFAULTS[n]) for n in KNOBS}
    for r in rejected:
        if r["name"] in KNOBS:
            would[r["name"]] = base[r["name"]]
    for names, test, why in CROSS:
        if not test(would):
            for n in names:
                if n in ok:
                    rejected.append({"name": n, "value": ok.pop(n), "reason": why})
                    would[n] = base[n]
    return ok, rejected


def read_overrides(path: Path) -> Tuple[Dict[str, Any], Optional[str]]:
    """The file's overrides as written ({} when there is no file), and a parse error if any."""
    try:
        text = path.read_text()
    except FileNotFoundError:
        return {}, None
    except OSError as e:
        return {}, f"cannot read {path}: {e}"
    if not text.strip():
        return {}, None
    try:
        data = json.loads(text)
    except ValueError as e:
        return {}, f"{path.name} is not valid JSON: {e}"
    if not isinstance(data, dict):
        return {}, f"{path.name} must hold one JSON object of setting names to values"
    return data, None


def write_overrides(path: Path, overrides: Dict[str, Any]) -> None:
    """Replace the file atomically (a reader never sees half a file)."""
    path = Path(path)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(dict(sorted(overrides.items())), indent=1) + "\n")
    os.replace(tmp, path)


def describe() -> List[dict]:
    """Every setting's metadata, for the dashboard."""
    return [{"name": k.name, "group": k.group, "default": k.default, "lo": k.lo, "hi": k.hi, "when": k.when,
             "unit": k.unit, "help": k.help, "kind": k.kind} for k in KNOBS_LIST]


class Tuning:
    """The values in force. `values` is replaced as a whole on every change, so readers in other threads always see
    one consistent set without a lock."""

    def __init__(self) -> None:
        self.values: Dict[str, Any] = dict(DEFAULTS)
        self.path: Optional[Path] = None
        self.overrides: Dict[str, Any] = {}             # the file's accepted overrides
        self.pending: Dict[str, Any] = {}               # round settings waiting for the next round
        self.rejected: List[dict] = []
        self.error: Optional[str] = None
        self.loaded_at: Optional[float] = None
        self.changes: collections.deque = collections.deque(maxlen=200)
        self.in_round = False
        self._stamp: Optional[Tuple[int, int]] = None
        self._checked = 0.0
        self._lock = threading.RLock()

    # ---- reading ------------------------------------------------------------------------------------------------------
    def get(self, name: str) -> Any:
        if self.path is not None:
            now = time.monotonic()
            if now - self._checked >= CHECK_EVERY_S:
                self._checked = now
                self._maybe_reload()
        return self.values[name]

    # ---- the file -----------------------------------------------------------------------------------------------------
    def use_file(self, path: Optional[Path] = None) -> None:
        """Main only: read overrides from the file from now on."""
        with self._lock:
            self.path = Path(path) if path is not None else TUNING_FILE
            self._stamp = None
            self._maybe_reload()

    def _maybe_reload(self) -> None:
        with self._lock:
            if self.path is None:
                return
            try:
                st = self.path.stat()
                stamp = (st.st_mtime_ns, st.st_size)
            except FileNotFoundError:
                stamp = (0, -1)
            except OSError:
                return
            if stamp == self._stamp:
                return
            self._stamp = stamp
            self._reload()

    def _reload(self) -> None:
        raw, error = read_overrides(self.path)
        self.loaded_at = time.time()
        if error:
            if error != self.error:
                log.warning("tuning: %s - the values in force are kept", error)
            self.error = error
            return
        self.error = None
        ok, rejected = validate(raw, self.values)
        for r in rejected:
            if r not in self.rejected:
                log.warning("tuning: %s = %r rejected: %s", r["name"], r["value"], r["reason"])
        self.overrides, self.rejected = ok, rejected
        self._settle()

    def _settle(self) -> None:
        """Apply the file's state: live settings now; round settings now between rounds, else at the next round."""
        rejected = {r["name"] for r in self.rejected}
        values = dict(self.values)
        pending: Dict[str, Any] = {}
        for name, k in KNOBS.items():
            want = self.overrides.get(name, values[name] if name in rejected else DEFAULTS[name])
            if want == values[name]:
                continue
            if k.when == ROUND and self.in_round:
                pending[name] = want
                continue
            self.changes.append({"t": time.time(), "name": name, "old": values[name], "new": want})
            log.info("tuning: %s %r -> %r", name, values[name], want)
            values[name] = want
        self.pending = pending
        self.values = values

    # ---- rounds -------------------------------------------------------------------------------------------------------
    def round_started(self) -> None:
        with self._lock:
            self._maybe_reload()
            self.in_round = False
            self._settle()                                     # pending round settings take effect now
            self.in_round = True

    def round_ended(self) -> None:
        with self._lock:
            self.in_round = False
            self._settle()

    # ---- status -------------------------------------------------------------------------------------------------------
    def snapshot(self) -> dict:
        with self._lock:
            return {"file": str(self.path) if self.path else None, "loaded_at": self.loaded_at,
                    "values": dict(self.values), "overrides": dict(self.overrides), "pending": dict(self.pending),
                    "rejected": list(self.rejected), "error": self.error, "changes": list(self.changes)[-50:]}


TUNING = Tuning()


def tune(name: str) -> Any:
    """The value in force of one setting."""
    return TUNING.get(name)
