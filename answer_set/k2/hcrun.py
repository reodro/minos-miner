"""GATK HaplotypeCaller evidence for GIAB answer sets (meta hc_evidence). Owner: answer_set.

The read scan, the read comparison and the haplotype genotyper miss ~10 inserted 1-bp indels per GIAB round; on 12
development rounds (2026-10-01 14:17 .. 2026-10-04 07:05, after g026's training, before the frozen set) a GATK run
found 4.3 of them per round (plus 1.8 2-3 bp indels and 5 SNPs), among ~590 false 1-bp calls - so its calls are
evidence for the model, never answers on their own.

One fixed configuration (GATK defaults, calling threshold 10 for recall), never a configuration being tuned: the
evidence must not follow what the miners tune. The region runs as `pieces` parallel intervals in the team's GATK
image; the calls are split and left-aligned like the candidates (util.normalize_vcf).
"""
import concurrent.futures as cf
import os
import subprocess
import uuid
from pathlib import Path

GATK_IMAGE = "broadinstitute/gatk:4.5.0.0"
VERSION = 1                       # the evidence a model is trained with (meta hc_evidence)
CALL_CONF = 10
PIECE_TIMEOUT_S = 900
FEATS = ["hc_called", "hc_qual", "hc_gq", "hc_af", "hc_hom", "hc_near"]


def _piece(bam, ref, chrom, x, y, out):
    bam, ref, out = Path(bam).resolve(), Path(ref).resolve(), Path(out).resolve()
    name = f"k2-hc-{os.getpid()}-{uuid.uuid4().hex[:10]}"
    cmd = ["docker", "run", "--rm", "--name", name, "--cpus=1", "--memory=4g", "--network", "none",
           "-u", f"{os.getuid()}:{os.getgid()}",
           "-v", f"{bam.parent}:/b:ro", "-v", f"{ref.parent}:/r:ro", "-v", f"{out.parent}:/o",
           GATK_IMAGE, "gatk", "--java-options", "-Xmx3g", "HaplotypeCaller",
           "-R", f"/r/{ref.name}", "-I", f"/b/{bam.name}", "-O", f"/o/{out.name}", "-L", f"{chrom}:{x}-{y}",
           "-stand-call-conf", str(CALL_CONF), "--native-pair-hmm-threads", "1"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=PIECE_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        subprocess.run(["docker", "rm", "-f", name], capture_output=True)
        raise RuntimeError(f"HaplotypeCaller {chrom}:{x}-{y} timed out")
    if r.returncode != 0:
        raise RuntimeError(f"HaplotypeCaller {chrom}:{x}-{y} failed: {r.stderr[-800:]}")
    return out


def run(bam, ref, chrom, a, b, workdir, pieces=8):
    """{(pos, ref, alt): {hc_qual, hc_gq, hc_af, hc_hom}} of the GATK calls in [a, b] (split, left-aligned, trimmed)."""
    import pysam
    from answer_set.k2.util import bcftools, normalize_vcf, trim
    workdir = Path(workdir); workdir.mkdir(parents=True, exist_ok=True)
    step = (b - a + pieces) // pieces
    spans = [(a + i * step, min(b, a + (i + 1) * step - 1)) for i in range(pieces) if a + i * step <= b]
    outs = [workdir / f"hc_p{i}.vcf.gz" for i in range(len(spans))]
    with cf.ThreadPoolExecutor(len(spans)) as ex:
        list(ex.map(lambda i: _piece(bam, ref, chrom, spans[i][0], spans[i][1], outs[i]), range(len(spans))))
    merged = workdir / "hc.vcf.gz"
    bcftools(f"cd /o && bcftools concat -a {' '.join(o.name for o in outs)} -Oz -o {merged.name} "
             f"&& bcftools index -f -t {merged.name}", [(str(workdir.resolve()), "/o")])
    norm = workdir / "hc.norm.vcf.gz"
    normalize_vcf(merged, norm, ref)
    calls = load_calls(norm, a, b)
    for o in outs:
        for f in (o, Path(str(o) + ".tbi")):
            f.unlink(missing_ok=True)
    return calls


def load_calls(norm, a=None, b=None):
    """The calls of a run's split, left-aligned VCF (work/hc.norm.vcf.gz) as run() returns them."""
    import pysam
    from answer_set.k2.util import trim
    calls = {}
    for rec in pysam.VariantFile(str(norm)):
        s = rec.samples[0]
        gt = s.get("GT")
        if not gt or any(g is None for g in gt) or not any(gt):
            continue
        ad, dp = s.get("AD"), s.get("DP")
        for i in {g for g in gt if g}:
            alt = rec.alleles[i]
            if alt == "*" or alt.startswith("<"):
                continue
            k = trim(rec.pos, rec.ref, alt)
            if a is not None and not a <= k[0] <= b:
                continue
            af = (ad[i] / sum(x for x in ad if x is not None)) if ad and ad[i] is not None and sum(
                x for x in ad if x is not None) else float("nan")
            calls[k] = {"hc_qual": float(rec.qual or 0.0), "hc_gq": float(s.get("GQ") or 0.0), "hc_af": af,
                        "hc_hom": float(len(set(gt)) == 1)}
    return calls
