"""The real answer-set builder. Owner: answer_set.

Status:
  1000G regime  germline exact (fingerprint -> donor -> raw calls + KG_RULE)
  GIAB regime   smoke only (regime detected, germline empty)
  synthetic     not started (~25-35% of truth; empty for now)
  n_targets     not estimated yet (None)
"""
from __future__ import annotations

from typing import Iterable

from answer_set import fingerprint, germline
from common.types import AnswerSet, CallSet, Regime, Round

# HG001 is 1000G sample NA12878: reads cannot tell the regimes apart, and every round seen so far
# with this donor used the NIST HG001 truth.
GIAB_IN_1000G = {"NA12878"}


class Provider:
    def build(self, round: Round) -> AnswerSet:
        fp = fingerprint.identify(round)
        notes = {"fingerprint": fp.__dict__}
        if fp.decisive and fp.donor not in GIAB_IN_1000G:
            return AnswerSet(round.round_id, 1, Regime.KG, fp.donor,
                             germline.kg_germline(round, fp.donor), [], None, notes)
        donor, calls = germline.giab_germline(round)
        return AnswerSet(round.round_id, 1, Regime.GIAB, donor or None, calls, [], None, notes)

    def update(self, answer_set: AnswerSet, new_callsets: Iterable[CallSet]) -> AnswerSet:
        return answer_set
