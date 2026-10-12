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
    "rf": "v10 RF_UCB: the forest, its proposers, QLock, Restart/Reinitialize, the endgame",
    "hill": "v11: several hills kept and climbed at once (RF_UCB per hill), cross-hill recombination",
    "tr": "v12: trust regions (TuRBO/CASMOPOLITAN style) searched at once, asynchronous Thompson sampling with the forest",
    "evaluator": "Evaluator: batches and routing on the fleet",
    "main": "Main: submissions",
    "answer_set": "Answer set: the targets the local score divides by",
}

KNOBS_LIST: List[Knob] = [
    Knob("candidates.portfolio_share", 0.4, 0.0, 1.0, ROUND, "share",
         "Share of wave 0 taken by portfolio seeds (the round's read-length profile first, never the round being "
         "replayed); the rest is a hypercube."),
    Knob("candidates.riders", 4, 0, 6, ROUND, "count", "Calling-threshold riders on each wave 0 point."),
    Knob("candidates.generations", 3, 0, 6, LIVE, "count",
         "Genetic generations after wave 0. The last one crosses the groups (core x fp, both x the others)."),
    Knob("candidates.launch_back", 0.8, 0.1, 1.0, LIVE, "share",
         "A generation launches once this share of the previous wave is back."),
    Knob("candidates.mutate_frac", 0.3, 0.0, 1.0, LIVE, "share",
         "Share of the children that get candidates.mutate_axes settings redrawn anywhere on their ladder."),
    Knob("candidates.mutate_axes", 2, 1, 10, LIVE, "count", "Settings redrawn in a mutated child."),
    Knob("candidates.group_size", 20, 2, 200, LIVE, "count", "Parents per group (both, core, fp), chosen far apart."),
    Knob("candidates.core_band", 0.01, 0.0, 0.2, LIVE, "score",
         "The core and fp groups take configs whose core is within this of the best core."),
    Knob("candidates.max_min", 30.0, 2.0, 50.0, LIVE, "min",
         "A safety limit of Step A (scaled to the round), not its schedule: every generation runs and the candidates "
         "are picked once the last one is candidates.launch_back back. On 20261008T175300Z a 15-minute limit cut "
         "generation 2 short and skipped generation 3 (real GATK waves took 4.5-10 minutes)."),
    Knob("candidates.workers_per_candidate", 300, 8, 5000, ROUND, "count",
         "One climbing candidate per this many workers (at least one)."),
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
    Knob("climb.wait_riders", 0, 0, 1, LIVE, "count",
         "1: a wave decides once its riders are back too. 0: once its grid points are back; riders that come back "
         "later count at the next decision."),
    Knob("climb.quorum", 1.0, 0.5, 1.0, LIVE, "share", "A wave counts as back at this share."),
    Knob("climb.end_min", 50.0, 10.0, 72.0, LIVE, "min", "Step B ends here (scaled to the round); Step C follows."),
    Knob("far.coarse", 25, 3, 500, LIVE, "count",
         "Step C: rungs tried beyond the nearest tried value on each side of a peak, at least (more when the fleet has room)."),
    Knob("far.second_margin", 0.1, 0.0, 2.0, LIVE, "score",
         "Step C: a second high point on a setting's profile within this of the peak is refined and climbed."),
    Knob("far.landings", 3, 0, 20, LIVE, "count", "Step C: second high points per peak that get a short climb."),
    Knob("far.landing_iterations", 2, 1, 10, LIVE, "count", "Step C: grids in a short climb."),
    Knob("far.pick_s", 90.0, 10.0, 600.0, LIVE, "s",
         "The highest-scoring config is named this many seconds before the search ends (not scaled)."),
    Knob("rf.trees", 200, 10, 1000, LIVE, "count", "Trees in the forest."),
    Knob("rf.min_leaf", 2, 1, 50, LIVE, "count", "Fewest configs in a leaf of a tree."),
    Knob("rf.kappa", 1.0, 0.0, 10.0, LIVE, "x",
         "UCB's exploration weight, in units of the spread of the predicted means among the 5% best predicted."),
    Knob("rf.kappa_dist", 1.0, 0.0, 10.0, LIVE, "x",
         "UCB's novelty weight: distance to the nearest evaluated config, in units of the spread of the predicted means "
         "among the 5% best predicted. A forest does not know where it has not looked; this term does."),
    Knob("rf.riders", 3, 0, 8, LIVE, "count", "Cheap settings moved one rung each, riding in every pick's GATK run."),
    Knob("rf.refit_per_fleet", 17, 1, 200, LIVE, "count",
         "A refit starts when the previous one has finished and workers/this many new results are in (legacy: 17)."),
    Knob("rf.fit_threads", 14, 1, 64, LIVE, "count", "Threads fitting the forest (main has 16 cores and uses few)."),
    Knob("rf.tree_sample", 0.3, 0.05, 1.0, LIVE, "share",
         "Each tree is grown on this random share of all the round's configs (1.0: every config, slower)."),
    Knob("rf.pool", 6000, 200, 50000, LIVE, "count", "Candidates QGA, GA and Inject propose per refit (Local adds its own)."),
    Knob("rf.share_qga", 0.55, 0.0, 1.0, LIVE, "share", "Share of rf.pool proposed by QGA (legacy: 0.55)."),
    Knob("rf.share_ga", 0.35, 0.0, 1.0, LIVE, "share", "Share of rf.pool proposed by GA-crossover (legacy: 0.35)."),
    Knob("rf.share_mutate", 0.20, 0.0, 1.0, LIVE, "share",
         "Share of the picks that are mutations (1-3 settings, 1-3 rungs) of the best-ranked proposals, unranked: the "
         "top of a forest's ranking is often near-duplicates from one leaf's box."),
    Knob("rf.qga_picks", 0.10, 0.0, 1.0, LIVE, "share",
         "Share of the picks guaranteed to QGA's best-ranked proposals (in one ranking its many-setting draws never "
         "reached the top)."),
    Knob("rf.share_sweep", 0.10, 0.0, 1.0, LIVE, "share",
         "Share of the picks from the sweep: every setting's line through the actual best (gaps and rungs beyond the "
         "nearest tried value, nearest first), unranked, rebuilt at every refit. Legacy: the sweeps beside the model "
         "produced every winning config after the opening."),
    Knob("rf.mutate_top", 50, 1, 5000, LIVE, "count", "The best-ranked proposals mutations are drawn from."),
    Knob("rf.centres", 4, 1, 20, LIVE, "count",
         "Distinct best configs (far apart by callset distance) QGA rotates toward and Local works around."),
    Knob("rf.share_inject", 0.10, 0.0, 1.0, LIVE, "share", "Share of the picks that are fresh hypercube configs (unranked)."),
    Knob("rf.qga_bins", 32, 2, 1000, LIVE, "count", "QGA bins per setting (a short ladder gets one bin per rung)."),
    Knob("rf.qga_theta", 0.10, 0.0, 1.0, LIVE, "share", "QGA rotation toward the elites at the start of the loop."),
    Knob("rf.qga_theta_end", 0.02, 0.0, 1.0, LIVE, "share", "QGA rotation by the endgame."),
    Knob("rf.qga_mut", 0.02, 0.0, 1.0, LIVE, "share", "Chance a setting's QGA vector swaps two bins at a rotation."),
    Knob("rf.qga_floor", 0.01, 0.0, 0.1, LIVE, "share", "Lowest probability of a QGA bin (per bin at 8 bins)."),
    Knob("rf.elites", 64, 4, 1000, LIVE, "count", "Best configs GA breeds from (plus the core/FP front)."),
    Knob("rf.front_window", 1.0, 0.0, 10.0, LIVE, "score",
         "The core/false-positive front among configs within this of the best joins the elites."),
    Knob("rf.ga_mut_rate", 0.10, 0.0, 1.0, LIVE, "share", "Chance each setting of a GA child mutates."),
    Knob("rf.ga_mut_rungs", 3, 1, 50, LIVE, "count", "A GA mutation moves a setting 1 to this many rungs."),
    Knob("rf.local_moves", 200, 0, 5000, LIVE, "count", "Random 2-3 setting moves around each of those elites."),
    Knob("rf.qlock_min_obs", 60, 8, 100000, LIVE, "count", "QLock waits for this many results."),
    Knob("rf.qlock_resolution", 0.003, 0.0, 1.0, LIVE, "score",
         "QLock freezes a setting whose sweep at the best moves the forest's prediction (with its spread) less than this."),
    Knob("rf.restart_waves", 2.5, 0.5, 20.0, LIVE, "x",
         "Restart after this many fleets' worth of results without a gain of rf.restart_min_gain; a second stall in a "
         "row Reinitializes."),
    Knob("rf.restart_min_gain", 0.005, 0.0, 1.0, LIVE, "score", "A gain that resets the stall count."),
    Knob("rf.portfolio_share", 0.15, 0.0, 1.0, ROUND, "share",
         "Share of the opening from the portfolio (small, so the forest is not fitted around other rounds' configs)."),
    Knob("rf.density_radius", 0.02, 0.0, 1.0, LIVE, "x",
         "Configs within this mean ladder distance of each other share their weight in the fit (0: equal weights)."),
    Knob("rf.opening_back", 0.5, 0.1, 1.0, LIVE, "share",
         "The loop starts once this share of the opening is back (its first refit runs before any pick)."),
    Knob("rf.no_restart_min", 10.0, 0.0, 60.0, LIVE, "min",
         "No Restart or Reinitialize within this many minutes (scaled) of the endgame."),
    Knob("rf.endgame_priority", 10, 0, 1000, LIVE, "count",
         "The endgame's proposals carry this priority, so they run ahead of the loop's leftovers."),
    Knob("rf.converge_gain", 0.02, 0.0, 1.0, LIVE, "score",
         "The loop has converged when its best gains less than this for rf.converge_min minutes."),
    Knob("rf.converge_min", 10.0, 1.0, 72.0, LIVE, "min", "See rf.converge_gain (scaled to the round's length)."),
    Knob("rf.converge_after_min", 20.0, 0.0, 72.0, LIVE, "min",
         "Convergence is not judged before this minute (scaled); when converged, the endgame starts at once."),
    Knob("rf.endgame_riders", 4, 0, 16, LIVE, "count", "Cheap riders per point of an endgame grid (v9's climb sends 16)."),
    Knob("hill.count", 4, 1, 10, LIVE, "count", "v11: hills funded at once (the best by score, see hill.drop)."),
    Knob("hill.callset_far", 0.03, 0.0, 1.0, LIVE, "x",
         "Two configs are on different hills when their callsets differ by at least this (class-weighted 1 - Jaccard "
         "from the sketches). Measured on saved callsets: variants of one hill 0.000-0.016, different trade-offs "
         "(core 0.865 vs 0.859, another PCR model) 0.067-0.095."),
    Knob("hill.setting_far", 0.05, 0.0, 1.0, LIVE, "x",
         "The same when a config has no sketch (a simulator): mean normalised setting difference."),
    Knob("hill.window", 1.5, 0.1, 10.0, LIVE, "score",
         "A config founds a new hill only within this of the best (lower ones join the nearest hill or are ignored)."),
    Knob("hill.track", 24, 4, 100, LIVE, "count", "Hills tracked (founders), of which hill.count are funded."),
    Knob("hill.points", 2000, 200, 50000, LIVE, "count", "The best this many results are grouped into hills."),
    Knob("hill.lead_share", 0.4, 0.0, 1.0, LIVE, "share",
         "Of the hills' picks, the leading hill's share; the other funded hills split the rest equally."),
    Knob("hill.share_cross", 0.10, 0.0, 0.5, LIVE, "share",
         "Picks that recombine two funded hills' bests by whole stage blocks (read, activity, assembly, PairHMM, "
         "genotyping), ranked by the forest."),
    Knob("hill.share_qga", 0.10, 0.0, 0.5, LIVE, "share", "Picks from QGA (rotated toward every funded hill): new hills."),
    Knob("hill.share_sweep", 0.20, 0.0, 1.0, LIVE, "share", "Within a hill: the next of every setting's line through its best."),
    Knob("hill.share_mutate", 0.20, 0.0, 1.0, LIVE, "share", "Within a hill: a mutation of one of its best-ranked proposals."),
    Knob("hill.drop", 0.5, 0.0, 10.0, LIVE, "score",
         "A funded hill this far below the leader with no gain for hill.patience_min loses its funding to the next."),
    Knob("hill.patience_min", 10.0, 1.0, 72.0, LIVE, "min", "See hill.drop (scaled)."),
    Knob("hill.valley_drop", 0.1, 0.0, 5.0, LIVE, "score",
         "Hill-valley test: two hills are one when a config between them scores within this of the lower; just above "
         "the near-peak roughness (one rung of padding moved the score 0.058 at the median on 153000Z)."),
    Knob("hill.valley_min", 5, 1, 100, LIVE, "count", "Configs between two hills, all lower, that confirm a valley."),
    Knob("hill.valley_mid", 0.25, 0.0, 0.5, LIVE, "share",
         "A config is between two hills when closer to each best than they are to each other and at least this share "
         "of that distance from both (not a near copy of one)."),
    Knob("hill.valley_hills", 6, 2, 24, LIVE, "count", "The valley test runs on the pairs among this many best hills."),
    Knob("hill.endgame", 4, 1, 10, LIVE, "count", "Endgame climbs: the best of this many funded hills."),
    Knob("hill.endgame_window", 0.3, 0.0, 10.0, LIVE, "score",
         "Endgame climbs only hills within this of the leader: a climb gains ~0.1-0.25 in the endgame, and every "
         "extra climb slows the others' waves."),
    Knob("tr.regions", 4, 1, 12, LIVE, "count", "v12: trust regions searched at once."),
    Knob("tr.length_init", 0.8, 0.01, 2.0, LIVE, "x", "A new region's length (TuRBO: 0.8)."),
    Knob("tr.length_min", 0.0078125, 0.0001, 0.5, LIVE, "x",
         "Below this a region is finished and restarts elsewhere (TuRBO: 0.5^7)."),
    Knob("tr.length_max", 1.6, 0.1, 4.0, LIVE, "x", "A region's largest length (TuRBO: 1.6)."),
    Knob("tr.step_results", 100, 10, 2000, LIVE, "count",
         "A region's step: this many of its picks back (asynchronous: TuRBO counts per batch)."),
    Knob("tr.step_max_min", 3.0, 0.2, 30.0, LIVE, "min",
         "A region whose step has not completed in this long (scaled) is judged on what it has: a starved or stuck "
         "region fails, shrinks and restarts (a simulated round left one with 28 results in 50 minutes)."),
    Knob("tr.min_share", 0.08, 0.0, 0.25, LIVE, "share",
         "Each region's floor share of the picks (its best candidates by the forest's mean); Thompson sampling decides "
         "the rest."),
    Knob("tr.sh_regions", 16, 1, 48, LIVE, "count",
         "Successive halving: regions opened after the opening, halved down to tr.regions (<= tr.regions: off)."),
    Knob("tr.sh_step_min", 4.0, 0.5, 30.0, LIVE, "min", "Time between halvings (scaled)."),
    Knob("tr.sh_gain_weight", 0.0, 0.0, 10.0, LIVE, "x",
         "A screened region's rank: best + this x its gain since opened. 0: best only (fleet 044200Z, 2026-10-11: at 1, "
         "far starts from ~75 gained +6-7 in 8 min and took 3 of 4 places; the best normal start was cut first)."),
    Knob("tr.sh_floor", 0.8, 0.0, 1.0, LIVE, "share", "While screening, picks shared out in equal turns across regions."),
    Knob("tr.sh_indel_keep", 1, 0, 8, LIVE, "count",
         "Each halving keeps this many regions with the best indel part at the reference FP, whatever their score."),
    Knob("tr.indel_fp_tol", 0.1, 0.0, 1.0, LIVE, "x",
         "Reference FP for the indel part: the best config's region FP x (1 + this)."),
    Knob("tr.sh_far_starts", 4, 0, 24, LIVE, "count",
         "Screening starts taken far in settings from the opening's top (part of tr.sh_regions)."),
    Knob("tr.sh_far_keep", 1, 0, 8, LIVE, "count", "Each halving keeps this many far starts, whatever their score."),
    Knob("tr.sh_far_setting", 0.2, 0.0, 1.0, LIVE, "x",
         "A far start's setting distance from every one of the opening's top tr.sh_home configs."),
    Knob("tr.sh_home", 20, 1, 500, LIVE, "count", "The opening's top configs that define our usual hill for far starts."),
    Knob("tr.share_lines", 0.25, 0.0, 0.8, LIVE, "share",
         "Picks that walk a region's lines: every setting's untried rungs through its centre, coarse to fine."),
    Knob("tr.importance_weights", 0, 0, 1, LIVE, "0/1",
         "Narrow a region's box on settings the forest finds important (off: on 044200Z the setting that needed a "
         "+17-rung move was kept on a short leash)."),
    Knob("tr.spent_min", 6.0, 1.0, 30.0, LIVE, "min",
         "In the last tr.no_restart_min, when no region has gained for this long (scaled), the endgame starts."),
    Knob("tr.reopen_length", 0.1, 0.001, 2.0, LIVE, "x", "Length of the region reopened on a best config left in none."),
    Knob("tr.success_gain", 0.005, 0.0, 1.0, LIVE, "score", "A step is a success when the region's best gained this."),
    Knob("tr.success_tol", 3, 1, 20, LIVE, "count", "Successes in a row that double a region's length (TuRBO: 3)."),
    Knob("tr.failure_tol", 2, 1, 20, LIVE, "count", "Failures in a row that halve a region's length."),
    Knob("tr.candidates", 2500, 200, 20000, LIVE, "count",
         "Candidates per region per refit (TuRBO: min(5000, max(2000, 200 x settings)))."),
    Knob("tr.importance_clip", 5.0, 1.0, 50.0, LIVE, "x",
         "A setting's box weight (1/importance, geometric mean 1) is clipped to [1/x, x]."),
    Knob("tr.new_far", 0.067, 0.0, 1.0, LIVE, "x",
         "A new region starts at least this far (callset distance) from every active and finished centre: the "
         "measured distance of different hills on saved callsets (0.067-0.095)."),
    Knob("tr.dup_far", 0.016, 0.0, 1.0, LIVE, "x",
         "Two active centres closer than this (callset distance; one hill's variants measured 0.000-0.016): the weaker "
         "region restarts."),
    Knob("tr.setting_new_far", 0.10, 0.0, 1.0, LIVE, "x", "tr.new_far without sketches (a simulator): setting distance."),
    Knob("tr.setting_dup_far", 0.02, 0.0, 1.0, LIVE, "x", "tr.dup_far without sketches: setting distance."),
    Knob("tr.no_restart_min", 20.0, 0.0, 72.0, LIVE, "min",
         "No new region within this of the endgame (scaled): a region opened at minute 40 cannot catch up (a simulated "
         "loss, 2026-10-10). Duplicates are still closed; their workers go to the other regions."),
    Knob("tr.endgame", 4, 1, 10, LIVE, "count", "Endgame climbs: at most this many region centres."),
    Knob("tr.endgame_window", 0.3, 0.0, 10.0, LIVE, "score", "Endgame climbs only centres within this of the best."),
    Knob("rf.endgame_min", 55.0, 10.0, 72.0, LIVE, "min",
         "From here (scaled to the round) v9's grid climb with lines, then Step C, on the best distinct configs."),
    Knob("rf.endgame_hills", 3, 1, 10, LIVE, "count",
         "Distinct configs the endgame climbs: the best score, the best core and the best FP, then the best distinct scores."),
    Knob("rf.ends_window", 0.5, 0.0, 5.0, LIVE, "score",
         "The best-core and best-FP configs are taken among configs within this of the best score (no junk ends)."),
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
