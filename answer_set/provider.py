"""The real answer-set builder. Owner: answer_set. Not implemented yet.

Plan (from the 2026-09-28 investigation):
  1. regime + donor: genotype a few hundred known sites from the BAM (pileup, no
     GATK) and match against the 7 GIAB NIST v4.2.1 sets and the 3,202-sample
     1000G 30x panel. Observed margin: best 0.83-0.93 Jaccard vs runner-up <=0.46.
  2. germline: that donor's calls in the region. GIAB rounds matched 1.000;
     1000G rounds 0.83-0.93 against the phased panel - exact source still to pin
     (hypothesis: NYGC raw per-sample genotypes).
  3. synthetic (~5K, 25-35% of truth, 97% not in the panel): candidates from
     callers, each labelled with p by a classifier trained on released rounds
     (read evidence, mappability, panel AF, caller agreement).
  4. n_targets: estimate N (mutations.vcf size, observed 434-779).
  5. update(): label newly surfaced sites with the same classifier, bump version.
"""
from __future__ import annotations

from typing import Iterable

from common.types import AnswerSet, CallSet, Round


class Provider:
    def build(self, round: Round) -> AnswerSet:
        raise NotImplementedError("answer_set.Provider.build")

    def update(self, answer_set: AnswerSet, new_callsets: Iterable[CallSet]) -> AnswerSet:
        raise NotImplementedError("answer_set.Provider.update")
