"""The evaluator pieces main uses, built here so runner settings and the scorer choice stay the evaluator's.
Owner: evaluator - main calls these and never constructs GatkRunner / LocalExecutor / scorers itself.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional

_FLEET: dict = {}


def local_executor(workers: Optional[int] = None, batching: bool = True):
    """GATK on this machine: 1 thread and 4 GB per job, the batching patch when built; `workers` defaults to the
    machine rule (common.settings.default_workers). batching=False runs every config solo (a batch runs its variants
    one after another, so a handful of configs finish sooner solo and in parallel)."""
    from evaluator.gatk import DEFAULT_PATCH_JAR, GatkRunner
    from evaluator.local import LocalExecutor
    runner = GatkRunner(threads=1, memory_gb=4,
                        patch_jar=DEFAULT_PATCH_JAR if os.path.exists(DEFAULT_PATCH_JAR) else None)
    return LocalExecutor(runner, workers=workers, max_batch=32 if batching else 1)


def callset_dir(round_tag: str) -> Path:
    """Where this machine's GATK runs write a round's VCFs (<key>_full.vcf.gz)."""
    from common.settings import SETTINGS
    return Path(SETTINGS.cache_dir) / "callsets" / round_tag


def scorer():
    """Exact scoring: upstream's hap.py path, one run per callset."""
    from evaluator.gatk import HappyScorer
    return HappyScorer()


def fleet_executor(inventory: str):
    """One fleet coordinator per process, started on first use: the sub-VPSs in `inventory` run GATK and score;
    main only proposes."""
    if "ex" not in _FLEET:
        from evaluator.fleet.coordinator import FleetExecutor
        from evaluator.gatk import DEFAULT_PATCH_JAR, GatkRunner
        planner = GatkRunner(patch_jar=DEFAULT_PATCH_JAR)          # only groups and prices families; never runs
        _FLEET["ex"] = FleetExecutor(planner, inventory=Path(inventory).expanduser()).start()
        logging.info("fleet: nodes from %s", inventory)
    return _FLEET["ex"]
