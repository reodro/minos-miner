"""The real answer-set builder. Owner: answer_set.

  donor      1000G: panel fingerprint (fingerprint.py); GIAB: genotype match against the 7 NIST sets (giab.py).
             NA12878 is HG001: every round seen with this donor used the NIST truth, so it goes the GIAB way.
  germline   exact lookup - 1000G: raw calls + KG_RULE (36/36 rounds exact); GIAB: NIST calls (24/24 exact)
  synthetic  round-BAM candidates that the donor's own reads (parent) show as reference (synthetic.py, parent.py)
  n_targets  estimated: N_PER_SYNTHETIC x the synthetic calls put in the truth (Round.n_targets wins if ever given)
The truth and target VCFs the scorers read are written to the cache and named in notes. mutations.vcf only feeds N
(the v2 core uses every truth record), so it just has to be non-empty.
Degraded paths: no donor -> Regime.UNKNOWN with no germline (main then submits a static baseline, so nothing else is
built); no parent reads -> every candidate goes into the truth unfiltered (real donor variants included: degraded,
but scorable). build() may still raise on missing data; main falls back to the static baseline then.

Status (2026-09-30), measured with test/answer_set_eval.py and test/answer_set_kpi.py:
  germline   exact: 1000G 36/36 rounds, GIAB 24/24 (donor ID 0.99+ vs runner-up <= 0.72)
  synthetic  1000G: SNP recall/precision ~0.998/0.997, indel ~0.99/0.985; GIAB HG005: SNP 0.995/0.985, indel 0.96/0.95
  KPI        same callsets scored under ours vs released truth: Spearman 0.98-0.998 and oracle gap 0 on 4 of 5
             rounds; 081700Z: Spearman 0.84, gap 0.21 points
  build      ~30-60 s for 1000G, ~110-140 s for GIAB (300x parent pileups dominate classify)
Scoring note: the v2 core uses EVERY truth record (germline + all synthetic); targets only give N.
Open work:
  - find the truth differences behind the 081700Z gap (which records, which difficulty classes)
  - GIAB indels: lower precision/recall than 1000G (different sequencing in the parent)
  - p=0.5 candidates (thin/ambiguous parent) are left out of the truth; decide per evidence
  - N estimate: median error 13%; measure its effect on rankings, try better predictors
  - speed: GIAB classify (pileups at 300x) and parent fetch; build() must stay within a few minutes
  - update(): live refinement during the round (parked by the user)
"""
from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Iterable, List

import pysam

from answer_set import fingerprint, germline, giab, parent, synthetic
from common.settings import SETTINGS
from common.types import AnswerSet, CallSet, Regime, Round, Variant

GIAB_IN_1000G = {"NA12878"}
TRUTH_MIN_P = 1.0          # synthetic calls at or above this go into the truth and target VCFs
# N (mutations.vcf entries) is unknown live and only sets the fp_per_target denominator. Fitted on 60 rounds:
# N ~ 0.1069 x synthetic pool size, median error 13% (per-Mb of region: 29%).
N_PER_SYNTHETIC = 0.1069           # fitted against the confident (parent-filtered) synthetic count


class Provider:
    def __init__(self, work_dir: str = ""):
        self.work_dir = Path(work_dir or os.path.join(SETTINGS.cache_dir, "answer_set"))

    def build(self, round: Round) -> AnswerSet:
        t, notes = {}, {}
        t0 = time.time()
        fp = fingerprint.identify(round)
        notes["fingerprint"] = fp.__dict__
        if fp.decisive and fp.donor not in GIAB_IN_1000G:
            regime, donor = Regime.KG, fp.donor
            germ = germline.kg_germline(round, donor)
        else:
            m = giab.identify(round)
            notes["giab_match"] = m.__dict__
            if m.decisive:
                regime, donor = Regime.GIAB, m.donor
                germ = giab.nist_germline(round, donor)
            else:
                t["donor+germline"] = time.time() - t0
                notes.update(seconds={k: round_(v) for k, v in t.items()}, region=round.region, reference=round.reference)
                return AnswerSet(round.round_id, 1, Regime.UNKNOWN, None, [], [], None, notes)   # main goes static
        t["donor+germline"] = time.time() - t0

        t0 = time.time()
        par = None
        url = parent.source(regime, donor) if donor else None
        if url:
            try:
                par = parent.fetch(round, url)
            except Exception as e:  # noqa: BLE001 - a missing parent degrades the answer set, never kills the round
                notes["parent_error"] = str(e)[:300]
        t["parent"] = time.time() - t0
        t0 = time.time()
        cands = synthetic.call_candidates(round)
        t["candidates"] = time.time() - t0
        t0 = time.time()
        syn, stats = synthetic.classify(round, cands, par, germ)
        if par is not None:
            par.close()
        t["classify"] = time.time() - t0
        notes.update(synthetic_stats=stats, seconds={k: round_(v) for k, v in t.items()}, region=round.region,
                     reference=round.reference, parent=url)

        d = self.work_dir / round.tag / "v1"
        d.mkdir(parents=True, exist_ok=True)
        confident = [s.variant for s in syn if s.p >= TRUTH_MIN_P]
        if par is None:                        # no parent: nothing is confident, but an empty target set is unscorable
            confident = [s.variant for s in syn]
            notes["degraded"] = "no parent reads: synthetic pool unfiltered (includes real donor variants)"
        notes["truth_vcf"] = write_vcf(d / "truth.vcf.gz", round, germ, confident)
        notes["mutations_vcf"] = write_vcf(d / "mutations.vcf.gz", round, [], confident)
        n_targets = round.n_targets or max(1, int(N_PER_SYNTHETIC * len(confident) + 0.5))
        notes["n_targets_source"] = "platform" if round.n_targets else "estimate"
        return AnswerSet(round.round_id, 1, regime, donor, germ, syn, n_targets, notes)

    def update(self, answer_set: AnswerSet, new_callsets: Iterable[CallSet]) -> AnswerSet:
        return answer_set


def round_(x: float) -> float:
    return round(x, 1)


def write_vcf(path: Path, round: Round, germ: List[Variant], syn: List[Variant]) -> str:
    """bgzipped + indexed VCF: germline records, then synthetic ones flagged SYNTHETIC, sorted by position."""
    length = None
    if os.path.exists(round.reference + ".fai"):
        with open(round.reference + ".fai") as f:
            for line in f:
                name, ln = line.split("\t")[:2]
                if name == round.chrom:
                    length = int(ln)
    header = ["##fileformat=VCFv4.2", '##FILTER=<ID=PASS,Description="All filters passed">',
              f"##contig=<ID={round.chrom}" + (f",length={length}>" if length else ">"),
              '##INFO=<ID=SYNTHETIC,Number=0,Type=Flag,Description="Injected variant (reconstructed)">',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"]
    rows = [(v.pos, f"{v.chrom}\t{v.pos}\t.\t{v.ref}\t{v.alt}\t50\tPASS\t.\tGT\t{v.gt}") for v in germ]
    rows += [(v.pos, f"{v.chrom}\t{v.pos}\t.\t{v.ref}\t{v.alt}\t50\tPASS\tSYNTHETIC\tGT\t{v.gt}") for v in syn]
    rows.sort(key=lambda r: r[0])
    plain = path.with_name(path.name[:-3])
    plain.write_text("\n".join(header + [r for _, r in rows]) + "\n")
    pysam.tabix_compress(str(plain), str(path), force=True)
    plain.unlink()
    pysam.tabix_index(str(path), preset="vcf", force=True)
    return str(path)
