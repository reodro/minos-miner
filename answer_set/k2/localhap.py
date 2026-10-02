"""Local haplotype check: which version of the local DNA does each read support?

The read scan (bcftools mpileup) counts letters column by column. Where two variants sit close
together, or an indel lies in a repeat, the same reads can be written several ways and the scan
often reports the wrong allele (a SNP instead of an indel, an indel shifted or of the wrong length)
and undercounts indel reads (the aligner clips or misplaces them).

Here, candidates within `gap` bp form a cluster; for each cluster we build the reference haplotype,
one haplotype per candidate allele and one per compatible pair of alleles, and align every read
that spans the cluster (whole read, soft-clipped ends included) to each haplotype (edit distance,
edlib). A read supports the haplotype it matches best (ties: no support). No variant caller is
used, so the answer set stays independent of the GATK configs being tuned.
"""
from collections import defaultdict
from itertools import combinations

import edlib
import pysam

FEATS = ["hap_alt", "hap_ref", "hap_n", "hap_frac", "hap_rank", "hap_margin", "hap_cluster"]


def _clusters(keys, gap):
    ks = sorted(keys)
    out, cur, end = [], [], -1
    for k in ks:
        s, e = k[0], k[0] + len(k[1]) - 1
        if cur and s - end > gap:
            out.append(cur)
            cur = []
        cur.append(k)
        end = max(end if cur[:-1] else e, e)
    if cur:
        out.append(cur)
    return out


def _apply(seq, w0, alleles):
    """seq = reference from 1-based position w0; alleles sorted, non-overlapping."""
    out, i = [], 0
    for pos, ref, alt in sorted(alleles):
        j = pos - w0
        out.append(seq[i:j])
        out.append(alt)
        i = j + len(ref)
    out.append(seq[i:])
    return "".join(out)


def hap_support(bam_path, fasta_path, chrom, keys, only=None, gap=10, flank=170, max_alleles=6,
                min_mapq=1):
    """{key: {feature: value}} for every key in `only` (default: all keys). `keys` = every candidate
    (including the person's own variants), used to form clusters and haplotypes."""
    only = set(keys) if only is None else set(only)
    bam = pysam.AlignmentFile(str(bam_path))
    fa = pysam.FastaFile(str(fasta_path))
    res = {}
    for cl in _clusters(keys, gap):
        if not any(k in only for k in cl):
            continue
        # keep the most relevant alleles: those asked about first, then the rest (bounded)
        cl = sorted(cl, key=lambda k: (k not in only, k))[:max_alleles]
        c0 = min(k[0] for k in cl)
        c1 = max(k[0] + len(k[1]) - 1 for k in cl)
        w0, w1 = max(1, c0 - flank), c1 + flank
        seq = fa.fetch(chrom, w0 - 1, w1).upper()
        haps = [((), seq)]
        for k in cl:
            haps.append(((k,), _apply(seq, w0, [k])))
        for x, y in combinations(sorted(cl), 2):
            if x[0] + len(x[1]) - 1 < y[0]:          # compatible (non-overlapping) pair
                haps.append(((x, y), _apply(seq, w0, [x, y])))
        support = defaultdict(int)
        n = 0
        for rd in bam.fetch(chrom, max(0, c0 - 2), c1 + 1):
            if rd.is_secondary or rd.is_supplementary or rd.is_duplicate or rd.mapping_quality < min_mapq:
                continue
            if rd.reference_start > c0 - 3 or rd.reference_end < c1 + 2:
                continue                              # must span the whole cluster (+2 bp)
            q = rd.query_sequence
            if not q:
                continue
            d = [edlib.align(q, h, mode="HW", task="distance")["editDistance"] for _, h in haps]
            best = min(d)
            if d.count(best) != 1:
                continue
            n += 1
            support[d.index(best)] += 1
        ref_n = support.get(0, 0)
        per = {k: 0 for k in cl}
        for hi, cnt in support.items():
            for k in haps[hi][0]:
                per[k] += cnt
        ranked = sorted(cl, key=lambda k: -per[k])
        for k in cl:
            if k not in only:
                continue
            others = [per[o] for o in cl if o != k]
            res[k] = {"hap_alt": per[k], "hap_ref": ref_n, "hap_n": n,
                      "hap_frac": per[k] / n if n else float("nan"),
                      "hap_rank": ranked.index(k) + 1,
                      "hap_margin": per[k] - (max(others) if others else 0),
                      "hap_cluster": len(cl)}
    return res
