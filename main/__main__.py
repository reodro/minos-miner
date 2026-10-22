"""python -m main replay <round_tag> [--minutes 15]
       python -m main live --wallet NAME --hotkey NAME [--demo]
       python -m main review <run_id|round_tag>

replay  the whole miner on a stored round, as if live: real answer set, real GATK (batched), real scoring, the
        search part's searcher, dry-run submissions, a deadline `--minutes` from now. Then reports what it would have
        submitted, releases the fleet, and reviews the run here (main.review): answer-set accuracy, TRUE scores, rank.
        `--no-review` stops after mining so a follow-up `python -m main review` can run while the next round takes
        the fleet.
live    the real miner (needs bittensor, httpx and a registered wallet; untested - try --demo first). Next to it
        main.review checks every 5 min whether a mined round's truth is out / the round is finalized, and reviews
        it then (answer-set accuracy, true scores, rank) - main's spare CPU is for this.
review  review one stored run now (replay does this itself; live does it in the background).

--fleet [INVENTORY] (both commands): main only proposes; GATK and scoring run on the sub-VPSs listed in the
        inventory (default: fleet.json in the repo folder), which main starts and controls over ssh.
"""
from __future__ import annotations

import argparse
import logging
import os
import time

from common.params import ParamRanges
from common.settings import SETTINGS, default_workers
from common.store import ResultStore, clear_round
from common.types import Round


def ranges() -> ParamRanges:
    try:
        return ParamRanges.from_platform()
    except Exception as e:  # noqa: BLE001 - offline: the checked-in snapshot
        logging.warning("live parameter ranges unavailable (%s); using the test fixture", e)
        return ParamRanges.from_file(os.path.join(os.path.dirname(__file__), "..", "tests", "fixtures",
                                                  "parameter_ranges.json"))


def reference_for(chrom: str) -> str:
    return SETTINGS.reference_fasta(chrom)


def _inventory_default():
    from evaluator.fleet.inventory import INVENTORY       # repo-folder fleet.json, or MINOS_FLEET_INVENTORY
    return INVENTORY


_BOARD = {}


def _board():
    """One status board per process (the dashboard reads what it writes)."""
    if "b" not in _BOARD:
        from main.status import StatusBoard
        _BOARD["b"] = StatusBoard()
    return _BOARD["b"]


def baselines(r: ParamRanges) -> list:
    """The known-good configs, best first; [0] is what a round with nothing scored submits. main.generalist_first
    (user, 2026-10-22): generalist g0 (/minos/data/portfolio/generalists_a1.json, mean gap 1.448 to the round's best over
    15 released rounds) before GATK's defaults (4-6 below). Round 03:46 fell back to defaults when its answer set came
    45 min late. A g0 outside the published ranges, or a missing file, leaves the defaults alone."""
    import json
    from common.tuning import tune
    from common.types import Config
    out = [r.defaults("gatk")]
    if int(tune("main.generalist_first")):
        try:
            with open(os.path.join(SETTINGS.data_dir, "portfolio", "generalists_a1.json")) as f:
                g0 = Config("gatk", {**out[0].params, **json.load(f)["configs"][0]})
            if not r.errors(g0):
                out.insert(0, g0)
        except (OSError, ValueError, KeyError, IndexError):
            pass
    return out


def make_parts(rnd: Round, r: ParamRanges, submitter, workers: int, fleet=None, scoring=None):
    """The pieces of one round, each from its own part's factory: main never builds another part's internals."""
    import evaluator
    import search
    import answer_set
    from common.tuning import TUNING, tune
    from main.loop import LoopSettings
    from main.miner import Parts
    TUNING.round_started()                    # settings that wait for a new round take effect now
    executor = fleet if fleet is not None else evaluator.local_executor(workers)
    os.makedirs(os.path.join(SETTINGS.cache_dir, "results"), exist_ok=True)
    store = ResultStore(os.path.join(SETTINGS.cache_dir, "results", f"{rnd.tag}.sqlite"))
    from main.answer_reuse import reusing
    return Parts(reusing(answer_set.builder()), executor, evaluator.scorer(), search.searcher(r, workers), submitter, r,
                 baselines(r), LoopSettings(resubmit_interval_s=tune("main.resubmit_interval_s"),
                                                    resubmit_margin=tune("main.resubmit_margin"),
                                                    stop_before_deadline_s=tune("main.stop_before_deadline_s")), store,
                 status=_board(), scoring=scoring)


def replay(tag: str, minutes: float, workers: int, fleet=None, scoring=None):
    from backtest.data import load
    from main.miner import mine_round
    from main.submitter import DryRunSubmitter
    rnd = load(tag, reference=reference_for(load(tag).chrom))
    if rnd.n_targets is None:                     # as live: N from the rounds revealed before this one started
        import dataclasses
        from main.n_estimate import recent_n
        n = recent_n(tag)
        if n:
            rnd = dataclasses.replace(rnd, n_targets=n, n_source="recent rounds")
    r = ranges()
    sub = DryRunSubmitter()
    parts = make_parts(rnd, r, sub, workers, fleet, scoring)
    t0 = time.time()
    mined = mine_round(rnd, parts, deadline_ts=t0 + minutes * 60)
    print(f"\n{tag}: {mined.mode} {mined.note} | {time.time() - t0:.0f}s")
    out = mined.outcome
    if not sub.sent:
        print("nothing submitted")
        return None
    if out:
        print(f"evaluated {out.evaluated}, failed {out.failed}, submissions {len(out.submissions)}; "
              f"local best {out.best.score if out.best else None}")
    return _board().state.get("round", {}).get("run_id")


def review_now(run_id: str) -> None:
    """Review + rank one run at once (its truth must be on disk or released) and print the result."""
    from main.review import Reviewer
    rv = Reviewer(board=_board())
    row = rv.review(run_id)
    if row is None:
        print(f"{run_id}: truth not released yet"); return
    a, srch = row.get("answer_set") or {}, row["search"]
    if a:
        g, y, n = a["germline_exact"], a["synthetic"], a["n_targets"]
        print(f"answer set vs truth: germline exact P {g['precision']} R {g['recall']} | synthetic P {y['precision']} "
              f"R {y['recall']} | N ours {n['ours']} truth {n['truth']}")
    for c in row["configs"]:
        print(f"  {c['role']:8s} {c['config']}  local {c['local']}  TRUE {c['true']}")
    print(f"final TRUE {srch.get('final_true')} | best true among them {srch.get('best_true')} "
          f"(pick gap {srch.get('pick_gap')}) | spearman {srch.get('spearman')} | rerun here {row['rerun']}")
    rank = rv.rank(run_id)
    if rank:
        print(f"rank {rank['rank']} of {rank['field']} (field best {rank['field_best']}, gap {rank['gap_to_best']})")
    else:
        print("rank: round not finalized yet")


def main() -> None:
    ap = argparse.ArgumentParser(prog="python -m main")
    sp = ap.add_subparsers(dest="cmd", required=True)
    a = sp.add_parser("replay"); a.add_argument("tag"); a.add_argument("--minutes", type=float, default=15)
    a.add_argument("--workers", type=int, default=default_workers())
    a.add_argument("--no-review", action="store_true",
                   help="do not review after mining; the fleet is free for the next round")
    b = sp.add_parser("live"); b.add_argument("--wallet", required=True); b.add_argument("--hotkey", required=True)
    b.add_argument("--demo", action="store_true"); b.add_argument("--workers", type=int, default=default_workers())
    c = sp.add_parser("review"); c.add_argument("run", help="run id (<tag>-<epoch>) or round tag (its newest run)")
    for p in (a, b):
        p.add_argument("--fleet", nargs="?", const=str(_inventory_default()),
                       metavar="INVENTORY",
                       help="run GATK and scoring on the sub-VPSs in the inventory (over ssh) instead of here")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    from common.limits import raise_open_files
    raise_open_files()                            # main's local workers and fleet sessions hold many files at once
    if args.cmd == "review":
        from main.review import find_run
        review_now(find_run(args.run)); return
    _board().kind = args.cmd
    from common.tuning import TUNING
    TUNING.use_file()                             # tuning.json overrides (the dashboard edits it); re-read on change
    import evaluator
    import signal
    import sys
    for sig in (signal.SIGTERM, signal.SIGHUP):   # kill / a closed terminal: run the finally below, which stops the
        signal.signal(sig, lambda signum, _: sys.exit(128 + signum))   # node agents (Python skips it by default)
    fleet = evaluator.fleet_executor(args.fleet) if args.fleet else None
    scoring = None if fleet is not None else evaluator.scoring_pool(args.workers)   # the nodes score on a fleet
    run_id = None
    try:
        if args.cmd == "replay":
            run_id = replay(args.tag, args.minutes, args.workers, fleet, scoring)
        else:
            from main.review import Reviewer
            Reviewer(board=_board()).start()                  # every 5 min: truth out? finalized? -> review / rank
            from main.miner import PlatformSubmitter, UpstreamPlatform, run_forever
            platform = UpstreamPlatform(args.wallet, args.hotkey, demo=args.demo)
            r = ranges()
            run_forever(platform, lambda rnd: make_parts(rnd, r, PlatformSubmitter(platform), args.workers,
                                                         fleet, scoring), reference_for)
    finally:
        if fleet is not None:
            fleet.close()                         # every node agent stops and kills its GATK containers
        if scoring is not None:
            scoring.close()
    if args.cmd == "replay" and run_id and not args.no_review:
        logging.info("fleet released; review runs here")
        review_now(run_id)


if __name__ == "__main__":
    main()
