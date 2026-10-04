"""Designed cheap families: genotyping and PairHMM variants that can share one GATK run.

Both strategies use this. The family is chosen for its score, then the executor batches it; we never drop a
parameter because it would be inconvenient to run, only because it cannot score or it is not in the cheap tiers.
"""
from __future__ import annotations

from typing import Callable, List

from common.types import Config
from search.space import Space


def designed_family(space: Space, base: Config, limit: int, accept: Callable[[Config], bool],
                    genotyping_only: bool = False) -> List[Config]:
    """Up to `limit` fresh cheap variants of `base`. Genotyping first; a third of the slots stay for PairHMM
    unless `genotyping_only` (the last minutes, when a PairHMM variant is too slow to finish)."""
    if limit <= 0:
        return []
    geno = [c for c in space.genotyping_variants(base) if accept(c)]
    if genotyping_only:
        return geno[:limit]
    pair = [c for c in space.pairhmm_variants(base) if accept(c)]
    n_pair = min(len(pair), max(1, limit // 3)) if pair else 0
    if n_pair >= limit:
        return pair[:limit]
    return geno[:limit - n_pair] + pair[:n_pair]
