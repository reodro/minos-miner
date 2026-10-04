"""The strategies main runs in a round. Owner: search - main only passes a name from STRATEGY_NAMES.

Adding or retuning a strategy (or running several side by side) is a change here, never in main.
"""
from __future__ import annotations

from typing import List

from common.params import ParamRanges
from search.base import Strategy, random_stream
from search.basin import basin_race
from search.regions import gradient_regions
from search.rf_ucb import rf_ucb
from search.tpe import tpe

STRATEGY_NAMES = ("rf_ucb", "pattern", "tpe", "regions", "random")   # the first one is the default


def strategies(name: str, ranges: ParamRanges, workers: int) -> List[Strategy]:
    """The lineup for `name`; `workers` is the executor's size on this machine (main's local workers)."""
    if name == "rf_ucb":
        # Tiered RF-UCB is the main lane. Multivariate TPE explores a different set of expensive keys beside it.
        return [rf_ucb(ranges, in_flight=workers, seed=0), tpe(ranges, in_flight=max(2, workers // 4), seed=1)]
    if name == "pattern":
        # Basin race: structurally different skeletons, cheap genotyping family on each traversal.
        # A submission sticks unless a different basin leads by the local-score band.
        return [basin_race(ranges, in_flight=workers, seed=0)]
    if name == "tpe":
        return [tpe(ranges, in_flight=workers)]
    if name == "regions":
        return [gradient_regions(ranges, in_flight=workers)]
    if name == "random":
        return [random_stream(ranges, in_flight=workers)]
    raise ValueError(f"unknown strategy {name!r}; known: {', '.join(STRATEGY_NAMES)}")
