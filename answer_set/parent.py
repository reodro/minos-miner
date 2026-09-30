"""The donor's own sequencing reads ("parent") over the round region, fetched remotely by region.

Round reads are real sequencing of the donor (downsampled, renamed), so a caller finds real donor variants that
are not in the truth; synthetic (injected) variants are absent from any real sequencing of the donor. Checked on
2 Mb of 20260928T081700Z: synthetic truth SNPs 0/1112 show the alt in the parent, germline 3189/3192 do.

Parents:
  1000G  the sample's NYGC 30x CRAM (1000 Genomes high-coverage, S3 mirror; ~3 s per 0.5 Mb)
  GIAB   the donor's NIST Illumina 100-300x BAM on NCBI (~8 s per 0.5 Mb at 300x)
The region is fetched as independent chunks in parallel; a site is looked up in the chunk containing it (every
read overlapping the chunk is in it, so coverage at the site is complete).
"""
from __future__ import annotations

import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

from common.settings import SETTINGS
from common.types import Regime, Round

SAMTOOLS_IMAGE = "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
KG_S3 = "https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage"
KG_INDEXES = {"1000G_2504_high_coverage.sequence.index": "data",
              "1000G_698_related_high_coverage.sequence.index": "additional_698_related/data"}
GIAB_BAMS = {
    "HG001": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam",
    "HG002": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam",
    "HG003": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.GRCh38.300x.bam",
    "HG004": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.GRCh38.300x.bam",
    "HG005": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/HG005.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam",
    "HG006": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG006_NA24694-huCA017E_father/NA24694_Father_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG006.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam",
    "HG007": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG007_NA24695-hu38168_mother/NA24695_Mother_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams/HG007.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam",
}


def kg_cram_url(sample: str) -> Optional[str]:
    """S3 URL of a 1000G sample's 30x CRAM, from the sequence indexes (scripts/fetch_reference_data.sh)."""
    d = Path(SETTINGS.kg_dir) / "cram_index"
    for name, sub in KG_INDEXES.items():
        path = d / name
        if not path.exists():
            continue
        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                url = line.split("\t", 1)[0]
                fname = url.rsplit("/", 1)[1]
                if fname == f"{sample}.final.cram":
                    run = url.rsplit("/", 2)[1]
                    return f"{KG_S3}/{sub}/{run}/{fname}"
    return None


def source(regime: Regime, donor: str) -> Optional[str]:
    if regime == Regime.KG:
        return kg_cram_url(donor)
    if regime == Regime.GIAB:
        return GIAB_BAMS.get(donor)
    return None


class Parent:
    """Parent reads over a region, as chunk BAMs; pileup counts at single sites."""

    def __init__(self, chunks: List[Tuple[int, int, str]], chrom: str):
        self.chunks, self.chrom = chunks, chrom            # (start, end, bam path), 1-based inclusive
        self._open: Dict[str, pysam.AlignmentFile] = {}

    def _bam(self, pos: int) -> Optional[pysam.AlignmentFile]:
        for a, b, path in self.chunks:
            if a <= pos <= b:
                if path not in self._open:
                    self._open[path] = pysam.AlignmentFile(path)
                return self._open[path]
        return None

    def counts(self, pos: int, ref: str, alt: str, min_bq: int = 13, min_mq: int = 20) -> Tuple[int, int]:
        """(reads supporting alt, reads covering) at a VCF-style site. SNP: base at pos. Indel (left-anchored,
        ref[0] == alt[0]): reads with an indel of that length right after pos."""
        bam = self._bam(pos)
        if bam is None:
            return 0, 0
        snp = len(ref) == 1 and len(alt) == 1
        want = len(alt) - len(ref)
        a = n = 0
        for col in bam.pileup(self.chrom, pos - 1, pos, truncate=True, min_base_quality=min_bq if snp else 0,
                              min_mapping_quality=min_mq, stepper="nofilter", ignore_orphans=False):
            for pr in col.pileups:
                aln = pr.alignment
                if aln.is_secondary or aln.is_supplementary or aln.is_duplicate or aln.is_qcfail:
                    continue
                if snp:
                    if pr.is_del or pr.is_refskip:
                        continue
                    n += 1
                    a += aln.query_sequence[pr.query_position] == alt
                else:
                    n += 1
                    a += pr.indel == want
        return a, n

    def close(self):
        for b in self._open.values():
            b.close()


def fetch(round: Round, url: str, out_dir: Optional[str] = None, chunk_bp: int = 500_000, workers: int = 8) -> Parent:
    """Parent reads over the round region (cached per round). Raises if any chunk fails."""
    d = Path(out_dir or os.path.join(SETTINGS.cache_dir, "parent", round.tag))
    d.mkdir(parents=True, exist_ok=True)
    spans = [(a, min(a + chunk_bp - 1, round.end)) for a in range(round.start, round.end + 1, chunk_bp)]
    ref = Path(round.reference).resolve()

    def one(span):
        a, b = span
        path = d / f"{round.chrom}_{a}_{b}.bam"
        if not (path.exists() and Path(str(path) + ".bai").exists()):
            tmp = path.with_suffix(".part.bam")
            cmd = ["docker", "run", "--rm", "-v", f"{d}:/out", "-v", f"{ref.parent}:/ref:ro", SAMTOOLS_IMAGE, "sh", "-c",
                   f"samtools view -T /ref/{ref.name} -b -o /out/{tmp.name} '{url}' {round.chrom}:{a}-{b} "
                   f"&& samtools index /out/{tmp.name}"]
            for attempt in range(4):                       # remote servers drop connections now and then
                p = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
                if p.returncode == 0:
                    break
                time.sleep(2 * (attempt + 1))
            else:
                raise RuntimeError(f"parent fetch {round.chrom}:{a}-{b} failed: {p.stderr[-300:]}")
            tmp.rename(path)
            Path(str(tmp) + ".bai").rename(str(path) + ".bai")
        return a, b, str(path)

    with ThreadPoolExecutor(workers) as pool:
        chunks = list(pool.map(one, spans))
    return Parent(chunks, round.chrom)
