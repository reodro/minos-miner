"""Which answer-set builder runs. Owner: answer_set - main only calls builder(); adding or switching builders is a
change here, never in main.

ANSWER_SET_BUILDER=<name> picks another builder for one process without changing the default, e.g. to replay a
round guided by the released truth:  ANSWER_SET_BUILDER=oracle python -m main replay <tag>
"""
from __future__ import annotations

import dataclasses
import os
from typing import Optional

from answer_set.base import AnswerSetProvider

BUILDER_NAMES = ("k2", "oracle")     # the first one is the default; "oracle" is for released rounds only


class SearchableOracle:
    """Test only: the released truth as the answer set main searches against (released rounds only).
    OracleProvider marks its answer set Regime.UNKNOWN (main would then submit a baseline without searching); this
    wrapper marks it GIAB so a replay searches exactly as with a real builder - the regime is only a gate in main."""

    def build(self, round):
        from answer_set.oracle import OracleProvider
        from common.types import Regime
        a = OracleProvider().build(round)
        return dataclasses.replace(a, regime=Regime.GIAB, notes=dict(a.notes, regime="placeholder (oracle test builder)"))

    def update(self, answer_set, new_callsets):
        return answer_set


def builder(name: Optional[str] = None) -> AnswerSetProvider:
    """The answer-set builder `name` (None: $ANSWER_SET_BUILDER, else the default)."""
    name = name or os.environ.get("ANSWER_SET_BUILDER") or BUILDER_NAMES[0]
    if name == "k2":
        from answer_set.k2 import K2Provider
        return K2Provider()
    if name == "oracle":
        return SearchableOracle()
    raise ValueError(f"unknown answer-set builder {name!r}; known: {', '.join(BUILDER_NAMES)}")
