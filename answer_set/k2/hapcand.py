"""Local-sequence helpers for the haplotype genotyper (hapgeno.py). Owner: answer_set.

The read scan and the read comparison see one allele at a time, in whatever spelling the alignment produced. Where
Minos's edits sit close together or inside a repeat, reads spell the same event differently (shifted, merged with a
neighbour, an insertion + deletion written as mismatches). Here: windows around seed positions, and each read's sequence
between two fixed anchor stretches flanking a window (one mismatch allowed, soft-clipped ends included) - independent of
how the aligner wrote it; local sequences carried by several reads (low-quality bases as wildcards); left-alignment of
events like bcftools norm.
"""
from collections import Counter

import edlib

from answer_set.k2.util import trim

ANCHOR, PAD, MAX_WINDOW = 12, 15, 80
MIN_READS, MIN_SHARE, MIN_ANCHORED = 2, 0.15, 4
MASK_Q = 20
MAX_INDEL = 10


def windows(seeds, cluster_gap=10):
    """Merge seed spans [(start, end), ...] (1-based, inclusive) into windows [start - PAD, end + PAD] of at most
    MAX_WINDOW bases; seeds within cluster_gap of each other share a window."""
    out = []
    for s, e in sorted(seeds):
        if out and s - out[-1][1] <= cluster_gap and max(e, out[-1][1]) - out[-1][0] <= MAX_WINDOW - 2 * PAD:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s - PAD, e + PAD) for s, e in out]


def _find(anchor, q, start=0):
    """First position after `start` where `anchor` matches q with at most one mismatch: (begin, end) or None."""
    a = edlib.align(anchor, q[start:], mode="HW", task="locations", k=1)
    if a["editDistance"] < 0 or not a["locations"]:
        return None
    b, e = a["locations"][0]
    if e is None:
        return None
    if b is None:                                          # edlib may report only the end of a match
        b = max(0, e - len(anchor) + 1)
    return start + b, start + e


def local_sequences(bam, fa, chrom, s, e, min_q=MASK_Q):
    """({sequence between the anchors: reads}, reference sequence of the window, anchored reads).
    Bases below min_q are wildcards: a read whose masked sequence matches one carried by other reads (same length,
    equal wherever it is not masked) supports that sequence - Minos's noise is mostly low quality and would otherwise
    split one real sequence into many versions, each too weak to count."""
    left = fa.fetch(chrom, s - 1 - ANCHOR, s - 1).upper()
    right = fa.fetch(chrom, e, e + ANCHOR).upper()
    ref = fa.fetch(chrom, s - 1, e).upper()
    reads = []
    for rd in bam.fetch(chrom, max(0, s - 1), e):
        if rd.is_secondary or rd.is_supplementary or rd.is_duplicate or not rd.query_sequence:
            continue
        q = rd.query_sequence
        L = _find(left, q)
        if L is None:
            continue
        R = _find(right, q, L[1] + 1)
        if R is None or R[0] - (L[1] + 1) > len(ref) + 2 * MAX_INDEL:
            continue
        i0, i1 = L[1] + 1, R[0]
        qual = rd.query_qualities
        masked = "".join(c if qual is None or qual[i0 + j] >= min_q else "N" for j, c in enumerate(q[i0:i1]))
        reads.append(masked)
    clean = Counter(x for x in reads if "N" not in x)
    seqs = Counter()
    for x in reads:
        if "N" not in x:
            seqs[x] += 1
            continue
        fits = [c for c in clean if len(c) == len(x) and all(a == "N" or a == b for a, b in zip(x, c))]
        if len(fits) == 1:                 # unambiguous: counts for that sequence
            seqs[fits[0]] += 1
        elif not fits and x.count("N") <= 2 and len(x) == len(ref) and all(a == "N" or a == b for a, b in zip(x, ref)):
            seqs[ref] += 1
    return seqs, ref, len(reads)


def left_align(fa, chrom, pos, ref, alt):
    """bcftools-norm-like: trim, then shift an indel left while the shared last base allows (minimal form)."""
    pos, ref, alt = trim(pos, ref, alt)
    while len(ref) != len(alt) and pos > 1 and ref[-1] == alt[-1]:
        b = fa.fetch(chrom, pos - 2, pos - 1).upper()
        ref, alt, pos = b + ref[:-1], b + alt[:-1], pos - 1
    return trim(pos, ref, alt)
