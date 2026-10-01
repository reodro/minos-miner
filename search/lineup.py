"""The strategies main runs in a round. Owner: search - main only passes a name from STRATEGY_NAMES.

Adding or retuning a strategy (or running several side by side) is a change here, never in main.
"""
from __future__ import annotations

from typing import List

from common.params import ParamRanges
from search.base import Strategy, random_stream
from search.rf_ucb import rf_ucb

STRATEGY_NAMES = ("rf_ucb", "random")          # the first one is the default


def strategies(name: str, ranges: ParamRanges, workers: int) -> List[Strategy]:
    """The lineup for `name`; `workers` is the executor's size on this machine (main's local workers)."""
    if name == "rf_ucb":
        return [rf_ucb(ranges, in_flight=workers)]
    if name == "random":
        return [random_stream(ranges, in_flight=workers)]
    raise ValueError(f"unknown strategy {name!r}; known: {', '.join(STRATEGY_NAMES)}")
