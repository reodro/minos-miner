"""python -m main replay <round_tag> [--minutes 15] [--strategy rf_ucb|random]
       python -m main live --wallet NAME --hotkey NAME [--demo]

replay  the whole miner on a stored round, as if live: real answer set, real GATK (batched), real scoring, the
        chosen strategy, dry-run submissions, a deadline `--minutes` from now. Then reports what it would have
        submitted, that config's TRUE score (released truth) and its rank among the real submissions.
live    the real miner (needs bittensor, httpx and a registered wallet; untested - try --demo first).

--fleet [INVENTORY] (both commands): main only proposes; GATK and scoring run on the sub-VPSs listed in the
        inventory (default ~/.minos/fleet.json), which main starts and controls over ssh.
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


def ranges() -> ParamRanges:
    try:
        return ParamRanges.from_platform()
    except Exception as e:  # noqa: BLE001 - offline: the checked-in snapshot
        logging.warning("live parameter ranges unavailable (%s); using the test fixture", e)
        return ParamRanges.from_file(os.path.join(os.path.dirname(__file__), "..", "tests", "fixtures",
                                                  "parameter_ranges.json"))


def reference_for(chrom: str) -> str:
    return os.path.join(SETTINGS.data_dir, "reference", chrom, f"{chrom}.fa")


_FLEET = {}


def fleet_executor(inventory: str):
    """One coordinator per process, started on first use: main only proposes; sub-VPS nodes run and score."""
    if "ex" not in _FLEET:
        from pathlib import Path
        from evaluator.fleet.coordinator import FleetExecutor
        from evaluator.gatk import DEFAULT_PATCH_JAR, GatkRunner
        planner = GatkRunner(patch_jar=DEFAULT_PATCH_JAR)          # only groups and prices families; never runs
        _FLEET["ex"] = FleetExecutor(planner, inventory=Path(inventory).expanduser()).start()
        logging.info("fleet: nodes from %s", inventory)
    return _FLEET["ex"]


def make_parts(rnd: Round, r: ParamRanges, strategy: str, submitter, workers: int, fleet=None):
    from answer_set import Provider
    from evaluator import DecisionCacheScorer, GatkRunner, LocalExecutor
    from evaluator.gatk import DEFAULT_PATCH_JAR
    from main.loop import LoopSettings
    from main.miner import Parts
    from search import random_stream, rf_ucb
    if fleet is not None:
        executor = fleet
    else:
        runner = GatkRunner(threads=1, memory_gb=4,
                            patch_jar=DEFAULT_PATCH_JAR if os.path.exists(DEFAULT_PATCH_JAR) else None)
        executor = LocalExecutor(runner, workers=workers)
    strat = rf_ucb(r, in_flight=workers) if strategy == "rf_ucb" else random_stream(r, in_flight=workers)
    os.makedirs(os.path.join(SETTINGS.cache_dir, "results"), exist_ok=True)
    store = ResultStore(os.path.join(SETTINGS.cache_dir, "results", f"{rnd.tag}.sqlite"))
    return Parts(Provider(), executor, DecisionCacheScorer(), [strat], submitter, r,
                 [r.defaults("gatk")], LoopSettings(resubmit_interval_s=300.0, stop_before_deadline_s=60.0), store)


def replay(tag: str, minutes: float, strategy: str, workers: int, fleet=None) -> None:
    from answer_set import OracleProvider
    from backtest.data import field_scores, load
    from backtest.kpi import rank_in_field
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
    final = sub.sent[-1][1]
    cs = out.callsets.get(final) if out else None
    if cs is None or not os.path.exists(cs.vcf):
        print("final config's callset is not on this machine (fleet): rerun it locally to get its TRUE score"); return
    true = parts.scorer.score(cs, OracleProvider().build(rnd)).score
    field = field_scores(tag)
    print(f"final submission {final}: TRUE score {true:.4f} -> rank {rank_in_field(true, field)} of {len(field) + 1} "
          f"(field best {max(field) if field else None})")


def main() -> None:
    ap = argparse.ArgumentParser(prog="python -m main")
    sp = ap.add_subparsers(dest="cmd", required=True)
    a = sp.add_parser("replay"); a.add_argument("tag"); a.add_argument("--minutes", type=float, default=15)
    a.add_argument("--strategy", default="rf_ucb", choices=["rf_ucb", "random"])
    a.add_argument("--workers", type=int, default=default_workers())
    b = sp.add_parser("live"); b.add_argument("--wallet", required=True); b.add_argument("--hotkey", required=True)
    b.add_argument("--demo", action="store_true"); b.add_argument("--workers", type=int, default=default_workers())
    for p in (a, b):
        p.add_argument("--fleet", nargs="?", const="~/.minos/fleet.json", metavar="INVENTORY",
                       help="run GATK and scoring on the sub-VPSs in the inventory (over ssh) instead of here")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    fleet = fleet_executor(args.fleet) if args.fleet else None
    try:
        if args.cmd == "replay":
            replay(args.tag, args.minutes, args.strategy, args.workers, fleet)
        else:
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
