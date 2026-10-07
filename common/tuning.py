"""Tunable settings of the search, the evaluator and main: defaults here, overrides in a JSON file on main.

The file (`TUNING_FILE`: tuning.json in the repo folder, git-ignored; MINOS_TUNING overrides the path) holds only
overrides, as flat names: {"scout.share": 0.2, "schedule.reduce_to_3_start_min": 30}. The dashboard edits it, and a
person may too. Only main turns the file on (`TUNING.use_file()`); tests, simulations and scripts run on defaults.

Reading: `tune("scout.share")` returns the value in force. The file is re-read when it changes (checked at most once a
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
    "opening": "Opening: wave 0 and the genetic generations",
    "domains": "Domains: the hills climbed at once",
    "schedule": "Schedule: minutes into a 72-minute round, scaled to the real round",
    "merges": "Merges",
    "sharing": "Sharing between domains",
    "scout": "Scout",
    "riders": "Riders",
    "climb": "Climb: one domain's ParamILS search",
    "pick": "Pick: the config main submits",
    "evaluator": "Evaluator: batches and routing on the fleet",
    "main": "Main: submissions",
}

KNOBS_LIST: List[Knob] = [
    Knob("opening.generations", 3, 0, 6, LIVE, "count", "Genetic generations after wave 0."),
    Knob("opening.launch_back", 0.6, 0.1, 1.0, LIVE, "share",
         "A generation launches once this share of the previous wave is back."),
    Knob("opening.launch_idle", 0.25, 0.0, 1.0, LIVE, "share",
         "Or once this share of the workers is idle, with half of the previous wave back."),
    Knob("opening.min_gain", 0.05, 0.0, 1.0, LIVE, "score",
         "A generation that gains less than this ends the opening."),
    Knob("opening.found_back", 0.8, 0.1, 1.0, LIVE, "share",
         "After the last generation, domains are founded once this share of it is back."),
    Knob("opening.max_min", 15.0, 2.0, 30.0, LIVE, "min", "Hard limit of the opening."),
    Knob("opening.islands", 1, 0, 1, LIVE, "count",
         "1: the generations breed on separate islands (one per domain, never mating across, a child that lands on "
         "another island's best is redrawn, equal children per island) and founding takes one domain per island. "
         "0: soft niches."),
    Knob("opening.steady", 1, 0, 1, LIVE, "count",
         "1: no generations - children are bred into idle workers as they free up (opening.generations x the workers "
         "in all); 0: generations."),
    Knob("opening.portfolio_share", 0.4, 0.0, 1.0, ROUND, "share",
         "Share of wave 0 taken by portfolio seeds (the legacy portfolio, then ours); the rest is a hypercube."),
    Knob("domains.max", 12, 1, 12, ROUND, "count",
         "Candidate centres the opening identifies; a small fleet gets one per 8 workers at most. Six of them climb. "
         "The rest are held and take a seat, one at a time, when a first scan gains under climb.decisive."),
    Knob("domains.stall_min", 10.0, 1.0, 40.0, LIVE, "min",
         "After its first scan, a basin whose best has not moved for this long is retired and the next candidate takes its seat."),
    Knob("domains.found_min_far", 6, 4, 12, LIVE, "count",
         "Domains founded from the opening are at least this many settings apart (4 = the hill rule)."),
    Knob("domains.worth_a_hill", 0.5, 0.0, 5.0, LIVE, "score",
         "No domain is founded more than this below the best score."),
    Knob("domains.min_early", 4, 1, 12, LIVE, "count",
         "Until the scheduled merges begin, fewer domains than this are refilled from the best off-hill finds."),
    Knob("domains.min_late", 2, 1, 6, LIVE, "count",
         "After the opening, a valley find refills up to this many basins. The basins stay separate until the final."),
    Knob("domains.importance_share", 0.0, 0.0, 2.0, LIVE, "x",
         "The hill rule counts only settings at least this times as important as the mean (the ranker's forests); 0: all."),
    Knob("domains.replace_margin", 0.05, 0.0, 1.0, LIVE, "score",
         "With every slot full, a find this much better than the weakest domain may replace it (after a valley test)."),
    Knob("schedule.reduce_to_3_start_min", 34.0, 0.0, 72.0, LIVE, "min",
         "By this minute, basins whose first scan gained under climb.decisive are replaced. Sooner if every first scan is already back."),
    Knob("schedule.reduce_to_3_by_min", 40.0, 0.0, 72.0, LIVE, "min",
         "Kept so the schedule stays in order. The replacement itself happens at reduce_to_3_start_min."),
    Knob("schedule.reduce_to_1_start_min", 50.0, 0.0, 72.0, LIVE, "min",
         "Kept so the schedule stays in order. The basins merge once, when the final starts, not at this minute."),
    Knob("schedule.reduce_to_1_by_min", 62.0, 0.0, 72.0, LIVE, "min",
         "The final starts here at the latest."),
    Knob("schedule.final_min", 10.0, 1.0, 30.0, LIVE, "min",
         "The final starts when this much time is left, if it has not started yet."),
    Knob("schedule.no_new_before_final_min", 5.0, 0.0, 30.0, LIVE, "min",
         "For this long before the final, a stalled basin retires and its workers join the basins already climbing. "
         "No new centre is founded."),
    Knob("schedule.pick_s", 90.0, 10.0, 600.0, LIVE, "s",
         "The highest-scoring config is named this many seconds before the search ends (not scaled)."),
    Knob("merges.on_evidence", 0, 0, 1, LIVE, "count",
         "1: two domains on one hill (by distance or a midpoint run) merge at any time and a frozen domain merges into a "
         "live one on its hill. 0: domains merge only when the schedule goes down to one, a flat basin is replaced once, "
         "a frozen domain retires, and a "
         "scout find is recombined with the nearest domain while it keeps climbing (merges.find_children)."),
    Knob("merges.find_children", 60, 1, 1000, LIVE, "count",
         "Children of a scout find recombined with the nearest domain's centre (merges.on_evidence 0)."),
    Knob("merges.retire", 1, 0, 1, LIVE, "count",
         "1: a merge retires the weaker domain, hands its differing settings to the stronger as moves (each alone, "
         "pairs of the first four), and queues one capped set of recombination children behind that domain's scan; "
         "0: the older full recombination, which pauses both domains."),
    Knob("merges.resolve_back", 0.8, 0.1, 1.0, LIVE, "share", "A merge resolves once this share of its children is back."),
    Knob("merges.close", 0.02, 0.0, 1.0, LIVE, "score",
         "Two nearby domains this close in score get a midpoint test."),
    Knob("merges.midpoint_retest_min", 10.0, 1.0, 60.0, LIVE, "min", "A pair's midpoint test is repeated after this long."),
    Knob("merges.exhaustive_max", 8, 1, 10, LIVE, "count",
         "Pairing tries every combination when the centres differ on this many settings or fewer."),
    Knob("merges.random_combos", 256, 0, 2048, LIVE, "count", "Otherwise this many random combinations."),
    Knob("merges.mutated", 50, 0, 500, LIVE, "count", "Mutated children added to each merge."),
    Knob("merges.children_min", 24, 1, 1000, LIVE, "count", "Children evaluated per merge: at least this many."),
    Knob("merges.children_share", 0.25, 0.01, 1.0, LIVE, "share", "Or this share of the workers, whichever is more."),
    Knob("sharing.every_min", 5.0, 1.0, 60.0, LIVE, "min", "How often domains offer each other their winning settings."),
    Knob("sharing.settings", 3, 0, 10, LIVE, "count", "Winning settings each domain offers."),
    Knob("scout.share", 0.0, 0.0, 0.6, LIVE, "share",
         "Share of the workers the scout sampler keeps on the whole space. 0: the sampler is off."),
    Knob("scout.late_share", 0.0, 0.0, 0.6, LIVE, "share",
         "The scout sampler's share once scheduled merges have begun. 0: the sampler is off."),
    Knob("scout.kappa", 1.0, 0.0, 5.0, LIVE, "x",
         "Exploration weight: candidates rank by model mean plus this times its spread."),
    Knob("riders.wave0", 4, 0, 8, ROUND, "count", "Threshold riders on each wave 0 point."),
    Knob("riders.global", 2, 0, 8, LIVE, "count", "Threshold riders on each generation and scout point."),
    Knob("riders.siblings", 2, 0, 2, LIVE, "count", "Siblings per climbing move: one genotyping, one PairHMM."),
    Knob("climb.accept_min", 0.005, 0.0, 0.1, LIVE, "score", "A gain smaller than this is noise."),
    Knob("climb.decisive", 0.05, 0.0, 1.0, LIVE, "score", "A gain this large moves the centre at once."),
    Knob("climb.reapply_full", 1, 0, 1, LIVE, "count",
         "1: a sweep's other winners are re-applied on the new centre alone, added one by one and in pairs of the top "
         "four (v3.2); 0: alone only."),
    Knob("climb.back_frac", 0.8, 0.1, 1.0, LIVE, "share", "A sweep counts as back at this share."),
    Knob("climb.equal_accepts", 2, 0, 10, LIVE, "count", "Better-or-equal acceptances per perturbation cycle."),
    Knob("climb.shake_k", [2, 3, 5, 8], 1, 20, LIVE, "count", "Settings shaken per perturbation, in turn."),
    Knob("climb.shake_points", 12, 1, 200, LIVE, "count", "Shaken points per perturbation, at least."),
    Knob("climb.scan_points", 9, 2, 30, LIVE, "count",
         "A new domain's first sweep: this many points over each setting's whole range."),
    Knob("climb.wide_max", 600, 0, 3000, LIVE, "count",
         "The wide round (all pairs of the strongest settings) after a coarse sweep without a gain: at most this many."),
    Knob("climb.max_chains", 300, 0, 300, LIVE, "count",
         "Parallel ParamILS chains per domain, on the workers its main line leaves idle."),
    Knob("climb.eager_wide_domains", 2, 0, 12, LIVE, "count",
         "With this many live domains or fewer, a domain's first coarse sweep includes the wide all-pairs round (one "
         "domain cannot fill the fleet otherwise); 0: the wide round only after a coarse sweep without a gain."),
    Knob("climb.chain_local", 12, 1, 100, LIVE, "count", "Local moves around each chain's shaken point."),
    Knob("climb.chain_restart", 0.2, 0.0, 1.0, LIVE, "share",
         "Chance that a chain restarts from the domain's centre instead of its own incumbent."),
    Knob("pick.band", 0.05, 0.0, 1.0, LIVE, "score", "Configs within this of the best are pick candidates."),
    Knob("pick.near", 2, 1, 6, LIVE, "count", "Neighbours are at most this many settings apart."),
    Knob("pick.min", 3, 1, 20, LIVE, "count", "A candidate needs this many scored neighbours."),
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
]
KNOBS: Dict[str, Knob] = {k.name: k for k in KNOBS_LIST}
assert len(KNOBS) == len(KNOBS_LIST) and all(k.group in GROUPS for k in KNOBS_LIST)
DEFAULTS: Dict[str, Any] = {k.name: (list(k.default) if k.kind == "list" else k.default) for k in KNOBS_LIST}

# Checks across settings: (the settings involved, test on the would-be values, why it fails).
CROSS: List[Tuple[Sequence[str], Callable[[Dict[str, Any]], bool], str]] = [
    (("schedule.reduce_to_3_start_min", "schedule.reduce_to_3_by_min", "schedule.reduce_to_1_start_min",
      "schedule.reduce_to_1_by_min"),
     lambda v: v["schedule.reduce_to_3_start_min"] <= v["schedule.reduce_to_3_by_min"]
     <= v["schedule.reduce_to_1_start_min"] <= v["schedule.reduce_to_1_by_min"],
     "the merge times must be in order: reduce to 3 start <= by <= reduce to 1 start <= by"),
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
