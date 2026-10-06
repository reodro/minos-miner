"""K2Provider: the model-based answer-set builder (answer_set.k2) behind the AnswerSetProvider interface.

  candidates  bcftools mpileup on the round BAM (every allele 2+ reads show; no calling)
  format      GIAB donor list support >= 0.95 -> giab track; else 1000G identify -> kg track; else unrecognised
  models      the track's ACTIVE version (<data_dir>/answer_set/models/<track>/ACTIVE): gradient-boosted trees
              decide P(in truth) and P(hom) per candidate; GIAB rounds copy the donor's NIST list and add model
              calls inside its BED
  output      truth.vcf.gz, mutations.vcf.gz (N=548 entries), scores.tsv.gz (P per candidate), report.json,
              in <cache_dir>/answer_set/k2/<round tag>/v1/

Scoring reads notes['truth_vcf'] / notes['mutations_vcf'], i.e. exactly the files the builder wrote. germline /
synthetic below are parsed from the same truth file (for status and review); they never feed the score.

Backtests and replays (round.released_dir set) use the leave-out fold model that never saw the round; live rounds
are in no training set and use the full active model.

Runtime needs: Docker (the bcftools 1.20 image the validator's scoring uses), internet (DRAGEN on S3, EBI GATK
calls, Google NovaSeq originals), and the inputs under SETTINGS.data_dir listed in answer_set/k2/util.py.
"""
from __future__ import annotations

import gzip
import json
import logging
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from answer_set.k2.build import build
from answer_set.k2.util import CACHE, trim
from common.types import AnswerSet, CallSet, Regime, Round, SyntheticCandidate, Variant
from common.vcf import records

log = logging.getLogger(__name__)


class K2Provider:
    def __init__(self, work_dir: str = "", threads: int = 8):
        self.work_dir = Path(work_dir) if work_dir else CACHE
        self.threads = threads

    def build(self, round: Round) -> AnswerSet:
        out = self.work_dir / round.tag / "v1"
        if out.exists():                       # never mix files of an earlier (possibly partial) build
            shutil.rmtree(out)
        build(round.bam, round.region, out, self.threads,
              exclude_round=round.tag if round.released_dir else None)
        answer = read_output(round, out)
        for w in answer.notes["warnings"]:
            log.warning("%s answer set: %s", round.tag, w)
        return answer

    def update(self, answer_set: AnswerSet, new_callsets: Iterable[CallSet]) -> AnswerSet:
        return answer_set


def read_output(round: Round, out: Path) -> AnswerSet:
    """AnswerSet from a build directory (truth.vcf.gz, mutations.vcf.gz, scores.tsv.gz, report.json)."""
    report = json.loads((out / "report.json").read_text())
    regime, donor = regime_of(report["mode"])
    p_of = read_scores(out / "scores.tsv.gz")
    germline: List[Variant] = []
    synthetic: List[SyntheticCandidate] = []
    for v, info in records(str(out / "truth.vcf.gz")):
        if "SYNTHETIC" in info.split(";"):
            ps = [p_of.get(trim(v.pos, v.ref, alt)) for alt in v.alt.split(",")]
            ps = [p for p in ps if p is not None]
            synthetic.append(SyntheticCandidate(v, min(ps) if ps else 1.0))
        else:
            germline.append(v)
    mutations = str(out / "mutations.vcf.gz")
    with gzip.open(mutations, "rt") as f:
        n_file = sum(1 for line in f if not line.startswith("#"))
    notes = {"source": "k2", "truth_vcf": str(out / "truth.vcf.gz"), "mutations_vcf": mutations,
             "scores_tsv": str(out / "scores.tsv.gz"), "region": round.region, "reference": round.reference,
             "model": report.get("model"), "threshold": report.get("threshold"), "mode": report["mode"],
             "n_targets_source": "platform" if round.n_targets else "builder (median)", "report": report,
             "warnings": warnings_of(report)}
    return AnswerSet(round.round_id, 1, regime, donor, germline, synthetic, round.n_targets or n_file, notes)


def warnings_of(report: dict) -> List[str]:
    """What made this build weaker than it should be (empty when nothing did)."""
    out = []
    if report.get("inputs", "ok") != "ok":
        out.append(f"inputs differ from the manifest: {report['inputs']}")
    rd = report.get("readdiff") or {}
    if rd.get("error"):
        out.append(f"read comparison failed, built without it (weaker on inserted indels): {rd['error']}")
    elif rd.get("failed_tries"):
        out.append(f"read comparison needed {len(rd['failed_tries']) + 1} tries: {rd['failed_tries']}")
    if (report.get("person_match") or {}).get("error"):
        out.append(f"1000G DRAGEN/panel download failed: {report['person_match']['error']}")
    if report.get("gatk_person_error"):
        out.append(f"1000G per-person GATK download failed: {report['gatk_person_error']}")
    return out


def regime_of(mode: str) -> Tuple[Regime, Optional[str]]:
    """report.json 'mode': 'GIAB:HG004', '1000G:HG03559[ WITHOUT DRAGEN ...]', or 'unrecognised (reads only)'."""
    if mode.startswith("GIAB:"):
        return Regime.GIAB, mode[5:].split()[0]
    if mode.startswith("1000G:"):
        return Regime.KG, mode[6:].split()[0]
    return Regime.UNKNOWN, None


def read_scores(path: Path) -> Dict[Tuple[int, str, str], float]:
    """{(pos, ref, alt) minimal representation: P(in truth)} from scores.tsv.gz."""
    out = {}
    with gzip.open(path, "rt") as f:
        next(f, None)
        for line in f:
            pos, ref, alt, p = line.rstrip("\n").split("\t")
            out[(int(pos), ref, alt)] = float(p)
    return out
