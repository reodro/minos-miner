"""Repeat-aware het/hom for 1000G inserted indels. Owner: answer_set.

In tandem repeats the kg het/hom model misjudges inserted indels (reads with another repeat length look like a second
copy, and the reverse). For a candidate indel of >= MIN_LEN bp inside a 1-6 bp tandem repeat spanning >= MIN_SPAN bp,
a second het/hom model (hom_repeat.joblib) replaces p_hom. It sees the version's own measurements plus:
  repeat unit length, repeat span, copies, event length in units, whether it is a whole number of units;
  reads judged by their DNA against the window with and without the event alone (carrying / reference / other /
  partial, and shares), and CIGAR stutter (spanning reads with an indel of another length in the repeat, distinct
  such lengths, reads with the event's length).
Trained on the true inserted repeat indels of k029's 109 rounds (test/answer_set/kg_hom_feats.py): held-out genotype
errors 393 -> 352 on k029's folds, of them >= 2 bp 258 -> 210; 1 bp got worse (135 -> 142), so 1 bp keeps k029's.
On only for a version whose meta sets hom_repeat; build.py applies it after the GATK genotype rule.
"""
from collections import Counter

import numpy as np

NEW = ["rep_unit", "rep_span", "rep_copies", "ev_units", "ev_whole_units", "rd_carry", "rd_ref", "rd_other", "rd_partial",
       "rd_carry_share", "rd_ref_share", "rd_other_share", "cig_stutter", "cig_stutter_lens", "cig_same"]
FL = 16
MIN_LEN, MIN_SPAN = 2, 8


def repeat_at(fa, chrom, k):
    """(unit length, span bp, copies) of the longest 1-6 bp tandem repeat starting right after the anchor base"""
    s = fa.fetch(chrom, k[0], k[0] + 120).upper()
    best = (0, 0, 0)
    for n in range(1, 7):
        u = s[:n]
        if len(u) < n:
            break
        c = 1
        while s[c * n:(c + 1) * n] == u:
            c += 1
        if c >= 2 and c * n > best[1]:
            best = (n, c * n, c)
    return best


def read_feats(bam, fa, chrom, k, span):
    p, ref, alt = k
    a, b = p - FL, p + len(ref) + FL + span
    refseq = fa.fetch(chrom, a - 1, b).upper()
    j = p - a
    with_ev = refseq[:j] + alt + refseq[j + len(ref):]
    L = len(alt) - len(ref)
    c = Counter(); st = Counter()
    for r in bam.fetch(chrom, a - 1, b):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.query_sequence:
            continue
        if not (r.reference_start + 1 <= a and r.reference_end >= b):
            c["partial"] += 1
            continue
        q, qs, qe, seq, last = r.query_sequence, r.query_alignment_start, r.query_alignment_end, [], None
        for qi, rp in r.get_aligned_pairs():
            if qi is not None and (qi < qs or qi >= qe):
                continue
            if rp is not None:
                if a <= rp + 1 <= b and qi is not None:
                    seq.append(q[qi])
                last = rp + 1
            elif qi is not None and last is not None and a <= last < b:
                seq.append(q[qi])
        s = "".join(seq)
        close = lambda h: len(h) == len(s) and sum(x != y for x, y in zip(h, s)) <= 1   # noqa: E731
        c["carry" if close(with_ev) else "ref" if close(refseq) else "other"] += 1
        rp0, lens = r.reference_start, []                      # CIGAR indels inside the repeat span
        for op, n in r.cigartuples:
            if op in (0, 7, 8):
                rp0 += n
            elif op == 1:
                if p - 1 <= rp0 <= p + span + 1:
                    lens.append(n)
            elif op == 2:
                if p - 1 <= rp0 <= p + span + 1:
                    lens.append(-n)
                rp0 += n
        if any(x != L for x in lens):
            st["stutter"] += 1
            for x in lens:
                if x != L:
                    st[("len", x)] += 1
        if L in lens:
            st["same"] += 1
    n = c["carry"] + c["ref"] + c["other"]
    sh = lambda x: c[x] / n if n else np.nan       # noqa: E731
    return [c["carry"], c["ref"], c["other"], c["partial"], sh("carry"), sh("ref"), sh("other"), st["stutter"],
            len([x for x in st if isinstance(x, tuple)]), st["same"]]


def new_feats(bam, fa, chrom, k):
    """the NEW measurements of candidate k, or None when it is not an indel of >= MIN_LEN bp in a repeat"""
    L = abs(len(k[1]) - len(k[2]))
    if L < MIN_LEN:
        return None
    u, span, cp = repeat_at(fa, chrom, k)
    if span < MIN_SPAN:
        return None
    return [u, span, cp, L / u, float(L % u == 0)] + read_feats(bam, fa, chrom, k, span)


def apply(p_hom, keys, X, person_hom, p_any, thr, bam, fa, chrom, model):
    """p_hom with the repeat model's value for every called candidate (p_any >= thr) that is an indel of >= MIN_LEN bp
    in a repeat and not one of the person's GATK genotypes (person_hom NaN) -> (p_hom, n replaced)"""
    out, n = p_hom.copy(), 0
    for i, k in enumerate(keys):
        if p_any[i] < thr or not np.isnan(person_hom[i]):
            continue
        f = new_feats(bam, fa, chrom, k)
        if f is None:
            continue
        out[i] = model.predict_proba(np.hstack([X[i], f])[None, :])[0, 1]
        n += 1
    return out, n
