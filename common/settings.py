"""Paths and endpoints. Override any field with an environment variable of the same name, upper-cased."""
from __future__ import annotations

import os
from dataclasses import dataclass, fields


@dataclass(frozen=True)
class Settings:
    data_dir: str = "/minos/data"                      # rounds/, giab/, panel/, cache/ - never in git
    minos_subnet_dir: str = "/minos/minos_subnet"      # upstream repo (scoring formula, templates)
    platform_url: str = "https://api.theminos.ai"
    files_url: str = "https://rounds.minosproxy.com"   # released round files

    @property
    def rounds_dir(self) -> str:
        return os.path.join(self.data_dir, "rounds")

    @property
    def kg_dir(self) -> str:
        """1000G data from scripts/fetch_reference_data.sh: raw/, phased/, mask/."""
        return os.path.join(self.data_dir, "1kg")

    @property
    def cache_dir(self) -> str:
        return os.path.join(self.data_dir, "cache")


def _load() -> Settings:
    kw = {}
    for f in fields(Settings):
        v = os.environ.get(f.name.upper())
        if v is not None:
            kw[f.name] = v
    return Settings(**kw)


SETTINGS = _load()
