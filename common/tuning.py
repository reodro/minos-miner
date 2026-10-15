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
    "far": "Step C: what Step B's grid did not look at, and the pick",
    "excess": "Step B: lanes that climb on core above the run's core-vs-FP curve, then on the score",
    "evaluator": "Evaluator: batches and routing on the fleet",
    "main": "Main: submissions",
    "answer_set": "Answer set: the targets the local score divides by",
}

KNOBS_LIST: List[Knob] = [
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
    Knob("climb.plateau_interval", 1, 0, 1, LIVE, "flag",
         "1: on an exact score tie, the plateau probe centres on the tied config nearest the per-setting intervals "
         "and sends the rung toward each interval first (search/plateau.py). It never accepts a move."),
    Knob("climb.tie_walk", 4, 0, 10, LIVE, "rungs",
         "A grid that gains nothing walks each one-setting line that tied the centre exactly, every untried rung out "
         "to this many rungs, in one more wave (lines toward search/plateau.py's interval first). 0: off (v9.18.2)."),
    Knob("climb.random_jumps", 20, 0, 500, LIVE, "count",
         "At a stall: random local jumps (2-4 settings, 2-5 rungs each); one that gains climb.accept_min resumes the grid."),
    Knob("climb.wait_riders", 0, 0, 1, LIVE, "count",
         "1: a wave decides once its riders are back too. 0: once its grid points are back; riders that come back "
         "later count at the next decision."),
    Knob("climb.quorum", 1.0, 0.5, 1.0, LIVE, "share", "A wave counts as back at this share."),
    Knob("climb.end_min", 55.0, 10.0, 72.0, LIVE, "min",
         "No longer the cut. Step B runs until the fine-tune, eight scaled minutes before the pick. "
         "Kept so an older override still loads."),
    Knob("climb.ga_quiet_min", 4.0, 1.0, 20.0, LIVE, "min",
         "After a GA child or a candidate moves a climber, the climbers send nothing for this many "
         "minutes (scaled) and the fleet breeds. Then one rung past the compass walks to the peak."),
    Knob("climb.ga_open_min", 8.0, 1.0, 20.0, ROUND, "min",
         "At the start of Step B the whole fleet breeds. Grids start once a free GA score beats the lowest "
         "climber, or when this many minutes (scaled) have passed. One wave: usually 2-4 minutes, 6-8 on a long region."),
    Knob("far.coarse", 25, 3, 500, LIVE, "count",
         "Step C: rungs tried beyond the nearest tried value on each side of a peak, at least (more when the fleet has room)."),
    Knob("far.second_margin", 0.1, 0.0, 2.0, LIVE, "score",
         "Step C: a second high point on a setting's profile within this of the peak is refined and climbed."),
    Knob("far.landings", 3, 0, 20, LIVE, "count", "Step C: second high points per peak that get a short climb."),
    Knob("far.landing_iterations", 2, 1, 10, LIVE, "count", "Step C: grids in a short climb."),
    Knob("excess.lanes", 0, 0, 4, LIVE, "count",
         "Step B climbers that climb on core above the run's core-vs-FP curve, then on the score (search/excess.py). "
         "0 (default): off, the engine is v9.18.1."),
    Knob("excess.start_far", 0.2, 0.0, 1.0, LIVE, "distance",
         "An excess lane starts at least this far in settings from the best config and from the other lanes."),
    Knob("excess.band_low", 0.95, 0.5, 1.0, LIVE, "x FP",
         "Excess measure: a config below this x the best score's FP is compared with the curve at this FP."),
    Knob("excess.band_high", 1.05, 1.0, 2.0, LIVE, "x FP",
         "Excess measure: a config above this x the best score's FP is compared with the curve here, minus the score's "
         "own price of the extra FP (x over_scale)."),
    Knob("excess.over_scale", 1.0, 0.0, 10.0, LIVE, "x",
         "Excess measure: multiplies the score's own price of FP above band_high (1 = exactly the score's exchange rate)."),
    Knob("excess.convert", 1, 0, 1, LIVE, "flag",
         "1: an excess lane climbs on the score from its J peak; 0: it stops at the J peak."),
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
    Knob("main.resubmit_interval_s", 300.0, 0.0, 3600.0, ROUND, "s", "Least time between two submissions."),
    Knob("main.resubmit_margin", 0.0, 0.0, 5.0, ROUND, "score", "A new best must gain this much to be resubmitted."),
    Knob("main.stop_before_deadline_s", 60.0, 10.0, 900.0, ROUND, "s", "The search stops this long before the deadline."),
    Knob("answer_set.n", 548, 50, 2000, ROUND, "count",
         "How many targets the next round's answer set writes. The local score divides false calls by this. "
         "A change made during a round waits for the next round."),
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
