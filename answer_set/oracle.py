"""Oracle provider: the platform's released truth, as an AnswerSet.

Two uses: (1) lets main and evaluator develop before the real builder exists;
(2) is the upper bound in answer-set backtests (oracle gap = 0 by definition).
Only valid for rounds whose files have been released.
"""
from __future__ import annotations

import gzip
import os
from typing import Iterable

from common.types import AnswerSet, CallSet, Regime, Round, SyntheticCandidate
from common.vcf import records


class OracleProvider:
    def build(self, round: Round) -> AnswerSet:
        d = round.released_dir
        if not d or not os.path.exists(os.path.join(d, "truth.vcf.gz")):
            raise FileNotFoundError(f"no released truth for {round.round_id} (released_dir={d!r})")
        germline, synthetic = [], []
        for v, info in records(os.path.join(d, "truth.vcf.gz")):
            if "SYNTHETIC" in info.split(";"):
                synthetic.append(SyntheticCandidate(v, 1.0))
            else:
                germline.append(v)
        with gzip.open(os.path.join(d, "mutations.vcf.gz"), "rt") as f:
            n_targets = sum(1 for line in f if not line.startswith("#"))
        return AnswerSet(round.round_id, 1, Regime.UNKNOWN, None, germline, synthetic, n_targets,
                         notes={"source": "oracle", "truth_vcf": os.path.join(d, "truth.vcf.gz"),
                                "mutations_vcf": os.path.join(d, "mutations.vcf.gz"), "region": round.region,
                                "reference": round.reference})

    def update(self, answer_set: AnswerSet, new_callsets: Iterable[CallSet]) -> AnswerSet:
        return answer_set
