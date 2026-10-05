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
    data_dir    INPUTS, shareable between clones on one machine: reference/, 1kg/, giab/, rounds/
    cache_dir   rebuildable caches of THIS clone (git-ignored, never synced to nodes): callsets/, results/,
                truth_slices/, dcache/ (scoring decision cache), rootfs/ (unpacked images), fleet/ (a node's rounds),
                answer_set/, selftest/ (the fleet self-test's BAM slice)
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
        """1000G data: phased/ (the phased panel, scripts/fetch_reference_data.sh) and the answer-set builder's
        gatk_sites/, index/, masks/ (answer_set/k2/util.py)."""
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


GB_PER_WORKER = 1.3     # measured: a 1-thread HaplotypeCaller JVM uses ~1.0-1.2 GB resident (its -Xmx4g heap is
                        # a cap it does not reach); + headroom for hap.py and the node agent (fleet, 2026-10-03)


def _ram_worker_cap() -> int:
    """Workers this machine's available RAM can hold: floor(MemAvailable / GB_PER_WORKER)."""
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    gb = int(line.split()[1]) / (1024 * 1024)
                    return max(1, math.floor(gb / GB_PER_WORKER))
    except OSError:
        pass
    return 10 ** 9


def default_workers(cores: int = 0) -> int:
    """GATK workers per machine (1 thread each): floor(cores * 0.9), and when the count is this machine's own,
    also floor(available RAM / GB_PER_WORKER). Same rule on every fleet machine; override with MINOS_WORKERS.

    An explicit `cores` argument is the CPU rule only, so a caller sizing a remote node is not capped by this
    machine's memory.
    """
    env = os.environ.get("MINOS_WORKERS")
    if env:
        return int(env)
    auto = not cores
    cores = cores or os.cpu_count() or 1
    workers = max(1, math.floor(cores * 0.9))
    if auto:
        workers = min(workers, _ram_worker_cap())
    return workers


def _load() -> Settings:
    kw = {}
    for f in fields(Settings):
        v = os.environ.get(f.name.upper())
        if v is not None:
            kw[f.name] = v
    return Settings(**kw)


SETTINGS = _load()
