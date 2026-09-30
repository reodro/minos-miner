"""Shared contracts between answer_set, evaluator and main.

Changing anything in this file changes the interface for all three parts:
it goes through a PR reviewed by every owner.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Dict, List, Optional

# The six v2 difficulty classes, in the validator's order.
CLASSES = ("snp_hom", "snp_het", "indel_1bp", "indel_2_3bp", "indel_4_7bp", "indel_8bp")


class Regime(str, Enum):
    """Where the round's germline truth comes from."""
    GIAB = "giab"        # NIST v4.2.1 calls of one of HG001-HG007
    KG = "1000g"         # one 1000 Genomes 30x sample
    UNKNOWN = "unknown"  # no donor identified; germline must be reconstructed


@dataclass(frozen=True)
class Round:
    round_id: str               # platform id, e.g. "2026-09-28T08:17:00+00:00"
    chrom: str
    start: int                  # 1-based inclusive, as in the platform region string
    end: int
    bam: str                    # local path to input.bam (indexed)
    reference: str              # local path to the chromosome FASTA (indexed)
    deadline: Optional[str] = None       # submission_end_time, ISO 8601
    released_dir: Optional[str] = None   # dir holding truth/mutations once revealed (backtest only)
    n_targets: Optional[int] = None      # platform's num_mutations (fp_per_target denominator), known live

    @property
    def region(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"

    @property
    def tag(self) -> str:
        """Directory-safe id used by the platform's file host: 20260928T081700Z."""
        return self.round_id[:19].replace("-", "").replace(":", "") + "Z"

    @staticmethod
    def parse_region(region: str):
        chrom, span = region.split(":")
        start, end = span.replace(",", "").split("-")
        return chrom, int(start), int(end)


@dataclass(frozen=True)
class Config:
    """A submittable variant-caller configuration."""
    tool: str                   # "gatk" | "deepvariant" | "bcftools"
    params: Dict[str, object] = field(default_factory=dict)

    @property
    def key(self) -> str:
        """Stable content hash: equal configs always share a key."""
        blob = json.dumps({"tool": self.tool, "params": self.params}, sort_keys=True)
        return hashlib.sha256(blob.encode()).hexdigest()[:16]

    def to_submission(self) -> dict:
        """The shape the platform accepts (accepted_top_level_keys)."""
        return {"tool": self.tool, f"{self.tool}_options": dict(self.params)}

    def to_dict(self) -> dict:
        return {"tool": self.tool, "params": dict(self.params)}

    @classmethod
    def from_dict(cls, d: dict) -> "Config":
        return cls(tool=d["tool"], params=dict(d["params"]))

    def __hash__(self):
        return hash(self.key)


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: str                     # normalised "0/1", "1/1", "1/2" ...


@dataclass(frozen=True)
class SyntheticCandidate:
    variant: Variant
    p: float                    # probability this is a real injected variant (1.0 = certain)


@dataclass
class AnswerSet:
    """Our best reconstruction of the platform's truth for one round.

    `version` increases on every update during a live round; every ScoreResult
    records which version it was scored against, so main can tell stale scores.
    """
    round_id: str
    version: int
    regime: Regime
    donor: Optional[str]                        # "HG002", "NA18878", or None
    germline: List[Variant]                     # treated as exact
    synthetic: List[SyntheticCandidate]         # uncertain, each with p
    n_targets: Optional[int]                    # N = entries in mutations.vcf (fp_per_target denominator); estimated live
    notes: Dict[str, object] = field(default_factory=dict)


FULL = "full"   # fidelity of a full-region evaluation; anything else is a sub-region like "chr20:1-1000000"


@dataclass(frozen=True)
class EvalJob:
    """A request to evaluate one config on one round. Search decides every field; the executor only runs it."""
    job_id: str
    round: Round
    config: Config
    fidelity: str = FULL
    priority: int = 0           # higher runs first
    tag: str = ""               # opaque to the executor: strategy lane, generation, ...


@dataclass(frozen=True)
class CallSet:
    """Output of running one config on one round."""
    round_id: str
    config_key: str
    vcf: str                    # path to the callset VCF(.gz)
    sha256: str
    runtime_s: float
    ok: bool = True
    error: Optional[str] = None
    job_id: str = ""
    fidelity: str = FULL
    score: Optional["ScoreResult"] = None   # set when the machine that ran it also scored it (fleet nodes); then
                                            # nobody else needs the VCF, which stays on that machine


@dataclass
class ClassCounts:
    tp: int = 0
    fp: int = 0
    fn: int = 0


@dataclass
class ScoreResult:
    round_id: str
    config_key: str
    answer_set_version: int
    counts: Dict[str, ClassCounts]              # keyed by CLASSES
    region_fp: int                              # all FPs in the challenge region (SNP + indel)
    n_targets: int
    ratios: Dict[str, float] = field(default_factory=dict)   # titv_/hethom_ truth/query for the plausibility gate
    core: Optional[float] = None
    score: Optional[float] = None               # 0-100, v2 formula; expected value when the answer set is uncertain
    gate_passed: Optional[bool] = None
    score_sd: Optional[float] = None            # spread of score over answer-set uncertainty (None = exact truth)
    fidelity: str = FULL
    job_id: str = ""

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "ScoreResult":
        d = dict(d)
        d["counts"] = {k: ClassCounts(**v) for k, v in d["counts"].items()}
        return cls(**d)
