"""Paths and endpoints. Override any field with an environment variable of the same name, upper-cased."""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, fields


@dataclass(frozen=True)
class Settings:
    data_dir: str = "/minos/data"                      # rounds/, giab/, panel/, cache/ - never in git
    minos_subnet_dir: str = "/minos/minos_subnet"      # upstream repo (scoring formula, templates)
    platform_url: str = "https://api.theminos.ai"

    @property
    def rounds_dir(self) -> str:
        return os.path.join(self.data_dir, "rounds")

    @property
    def kg_dir(self) -> str:
        """1000G data from scripts/fetch_reference_data.sh: raw/, phased/, cram_index/."""
        return os.path.join(self.data_dir, "1kg")

    @property
    def cache_dir(self) -> str:
        return os.path.join(self.data_dir, "cache")


def default_workers(cores: int = 0) -> int:
    """GATK workers per machine (1 thread each): floor(cores * 0.9 - 2). The other ~10% of cores score
    (default_score_threads) and run the system. Same rule on every fleet machine; override with MINOS_WORKERS."""
    env = os.environ.get("MINOS_WORKERS")
    if env:
        return int(env)
    cores = cores or os.cpu_count() or 1
    return max(1, math.floor(cores * 0.9 - 2))


def default_score_threads(cores: int = 0) -> int:
    """Scorer processes per fleet node: ceil(cores * 0.1), rounded up (36 cores -> 4): a scorer short stalls GATK, a
    spare one costs little. Override with MINOS_SCORE_THREADS."""
    env = os.environ.get("MINOS_SCORE_THREADS")
    if env:
        return int(env)
    cores = cores or os.cpu_count() or 1
    return max(1, math.ceil(cores * 0.1))


def _load() -> Settings:
    kw = {}
    for f in fields(Settings):
        v = os.environ.get(f.name.upper())
        if v is not None:
            kw[f.name] = v
    return Settings(**kw)


SETTINGS = _load()
