"""python -m main replay <round_tag> [--minutes 15] [--strategy NAME]      (names: search.STRATEGY_NAMES)
       python -m main live --wallet NAME --hotkey NAME [--demo]
       python -m main review <run_id|round_tag>

replay  the whole miner on a stored round, as if live: real answer set, real GATK (batched), real scoring, the
        chosen strategy, dry-run submissions, a deadline `--minutes` from now. Then reports what it would have
        submitted, and reviews the run (main.review): answer-set accuracy, TRUE scores, rank.
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
from common.store import ResultStore
from common.types import Round
from search.lineup import STRATEGY_NAMES


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


def make_parts(rnd: Round, r: ParamRanges, strategy: str, submitter, workers: int, fleet=None):
    """The pieces of one round, each from its own part's factory: main never builds another part's internals."""
    import evaluator
    import search
    import answer_set
    from main.loop import LoopSettings
    from main.miner import Parts
    executor = fleet if fleet is not None else evaluator.local_executor(workers)
    os.makedirs(os.path.join(SETTINGS.cache_dir, "results"), exist_ok=True)
    store = ResultStore(os.path.join(SETTINGS.cache_dir, "results", f"{rnd.tag}.sqlite"))
    return Parts(answer_set.builder(), executor, evaluator.scorer(), search.strategies(strategy, r, workers), submitter, r,
                 [r.defaults("gatk")], LoopSettings(resubmit_interval_s=300.0, stop_before_deadline_s=60.0), store,
                 status=_board())


def replay(tag: str, minutes: float, strategy: str, workers: int, fleet=None) -> None:
    from backtest.data import load
    from main.miner import mine_round
    from main.submitter import DryRunSubmitter
    rnd = load(tag, reference=reference_for(load(tag).chrom))
    r = ranges()
    sub = DryRunSubmitter()
    parts = make_parts(rnd, r, strategy, sub, workers, fleet)
    t0 = time.time()
    mined = mine_round(rnd, parts, deadline_ts=t0 + minutes * 60)
    print(f"\n{tag}: {mined.mode} {mined.note} | {time.time() - t0:.0f}s")
    out = mined.outcome
    if not sub.sent:
        print("nothing submitted"); return
    if out:
        print(f"evaluated {out.evaluated}, failed {out.failed}, submissions {len(out.submissions)}; "
              f"local best {out.best.score if out.best else None}")
    review_now(_board().state.get("round", {}).get("run_id"))


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
    a.add_argument("--strategy", default=STRATEGY_NAMES[0], choices=STRATEGY_NAMES)
    a.add_argument("--workers", type=int, default=default_workers())
    b = sp.add_parser("live"); b.add_argument("--wallet", required=True); b.add_argument("--hotkey", required=True)
    b.add_argument("--demo", action="store_true"); b.add_argument("--workers", type=int, default=default_workers())
    c = sp.add_parser("review"); c.add_argument("run", help="run id (<tag>-<epoch>) or round tag (its newest run)")
    for p in (a, b):
        p.add_argument("--fleet", nargs="?", const=str(_inventory_default()),
                       metavar="INVENTORY",
                       help="run GATK and scoring on the sub-VPSs in the inventory (over ssh) instead of here")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    if args.cmd == "review":
        from main.review import find_run
        review_now(find_run(args.run)); return
    _board().kind = args.cmd
    import evaluator
    fleet = evaluator.fleet_executor(args.fleet) if args.fleet else None
    try:
        if args.cmd == "replay":
            replay(args.tag, args.minutes, args.strategy, args.workers, fleet)
        else:
            from main.review import Reviewer
            Reviewer(board=_board()).start()                  # every 5 min: truth out? finalized? -> review / rank
            from main.miner import PlatformSubmitter, UpstreamPlatform, run_forever
            platform = UpstreamPlatform(args.wallet, args.hotkey, demo=args.demo)
            r = ranges()
            run_forever(platform, lambda rnd: make_parts(rnd, r, "rf_ucb", PlatformSubmitter(platform), args.workers,
                                                         fleet), reference_for)
    finally:
        if fleet is not None:
            fleet.close()                         # every node agent stops and kills its GATK containers


if __name__ == "__main__":
    main()
