"""Paths and endpoints. Override any field with an environment variable of the same name, upper-cased
(DATA_DIR, CACHE_DIR, STATUS_DIR, MINOS_SUBNET_DIR, PLATFORM_URL)."""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, fields
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Settings:
    """Where things live. Three kinds, never mixed:
    data_dir    INPUTS, shareable between clones on one machine: reference/, 1kg/, panel/, giab/, rounds/,
                truths/ (released truth + mutations VCFs of rounds whose BAM was never stored)
    cache_dir   rebuildable caches of THIS clone (git-ignored, never synced to nodes): callsets/, results/, parent/,
                truth_slices/, answer_set/, selftest/ (the fleet self-test's BAM slice)
    status_dir  what THIS clone's miner writes for the dashboard and the review (git-ignored, never synced)"""
    data_dir: str = "/minos/data"
    cache_dir: str = str(REPO / "cache")
    status_dir: str = str(REPO / "status")
    minos_subnet_dir: str = "/minos/minos_subnet"      # upstream repo (scoring formula, templates)
    platform_url: str = "https://api.theminos.ai"

    @property
    def rounds_dir(self) -> str:
        """Round inputs: <tag>/input.bam (+ truth.vcf.gz, mutations.vcf.gz, scores.json once released)."""
        return os.path.join(self.data_dir, "rounds")

    @property
    def kg_dir(self) -> str:
        """1000G data from scripts/fetch_reference_data.sh: raw/, phased/, cram_index/."""
        return os.path.join(self.data_dir, "1kg")

    @property
    def giab_dir(self) -> str:
        """GIAB NIST truth sets: HG00x.vcf.gz (+ .bed)."""
        return os.path.join(self.data_dir, "giab")

    @property
    def reference_dir(self) -> str:
        """GRCh38 per chromosome: <chrom>/<chrom>.fa (+ .fai, .dict, .sdf) from scripts/fetch_genome.sh."""
        return os.path.join(self.data_dir, "reference")

    def reference_fasta(self, chrom: str) -> str:
        return os.path.join(self.reference_dir, chrom, f"{chrom}.fa")


def default_workers(cores: int = 0) -> int:
    """GATK workers per machine (1 thread each): floor(cores * 0.9). Same rule on every fleet machine; override with
    MINOS_WORKERS."""
    env = os.environ.get("MINOS_WORKERS")
    if env:
        return int(env)
    cores = cores or os.cpu_count() or 1
    return max(1, math.floor(cores * 0.9))


def _load() -> Settings:
    kw = {}
    for f in fields(Settings):
        v = os.environ.get(f.name.upper())
        if v is not None:
            kw[f.name] = v
    return Settings(**kw)


SETTINGS = _load()
