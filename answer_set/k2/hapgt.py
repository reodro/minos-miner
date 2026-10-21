"""Haplotype-aware het/hom for 1000G inserted indels. Owner: answer_set.

Most het/hom errors on inserted indels sit next to another variant (62% of k047's wrong calls on its live rounds,
13% of its right ones): a read carrying the indel AND its neighbour does not match "the indel alone", so a hom
indel next to a het neighbour looks het, and the reverse. Here a spanning read carries the event if its DNA in the
window matches any local haplotype holding the event (the event plus any subset of its neighbours), and is
"without" if it matches one without it. The neighbours are what the build knows: the person's own records and the
other candidates it calls.

A het/hom model (hom_hap.joblib) on the version's own measurements + homrep's (repeat shape, event-alone read
judgement) + these replaces p_hom for every called indel candidate that is not one of the person's records.
Held-out on k047's 102 rounds (test/answer_set/gt_study): genotype errors on true inserted indels 544 (k047's
model) -> 510 (+ homrep's measurements) -> 476 (+ haplotype reads); with neighbours 306 -> 242.
On only for a version whose meta sets hom_hap.
"""
import bisect
import itertools

import numpy as np

from answer_set.k2.homrep import read_feats, repeat_at

FL = 16
WIN = 30
MAX_NEIGH = 5
HAP = ["hap_carry", "hap_without", "hap_other", "hap_share", "hap_n_neigh", "hap_n_person", "hap_n_called"]


def _window_seq(r, a, b):
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
    return "".join(seq)


def _edit(refseq, a, vs):
    """refseq (from position a) with the non-overlapping variants vs applied; None if two overlap"""
    out, cur = [], a
    for p, ref, alt in sorted(vs):
        if p < cur:
            return None
        out.append(refseq[cur - a:p - a]); out.append(alt); cur = p + len(ref)
    out.append(refseq[cur - a:])
    return "".join(out)


def hap_counts(bam, fa, chrom, k, span, neigh):
    """{'carry', 'without', 'other'}: spanning reads matching (<= 1 mismatch) a local haplotype with event k, one
    without it, or neither; haplotypes = the event or not, plus any subset of the MAX_NEIGH nearest neighbours"""
    p, ref, alt = k
    neigh = sorted(neigh, key=lambda n: abs(n[0] - p))[:MAX_NEIGH]
    a = min([p - FL] + [n[0] - 2 for n in neigh])
    b = max([p + len(ref) + FL + span] + [n[0] + len(n[1]) + 2 for n in neigh])
    refseq = fa.fetch(chrom, a - 1, b).upper()
    with_h, without_h = set(), set()
    for m in range(len(neigh) + 1):
        for sub in itertools.combinations(neigh, m):
            s1, s0 = _edit(refseq, a, list(sub) + [k]), _edit(refseq, a, list(sub))
            if s1:
                with_h.add(s1)
            if s0:
                without_h.add(s0)
    c = {"carry": 0, "without": 0, "other": 0}
    for r in bam.fetch(chrom, a - 1, b):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.query_sequence:
            continue
        if not (r.reference_start + 1 <= a and r.reference_end >= b):
            continue
        s = _window_seq(r, a, b)
        close = lambda H: any(len(h) == len(s) and sum(x != y for x, y in zip(h, s)) <= 1 for h in H)  # noqa: E731
        c["carry" if close(with_h) else "without" if close(without_h) else "other"] += 1
    return c


def feats(bam, fa, chrom, k, neigh):
    """homrep's measurements (repeat shape + event-alone reads) + the haplotype-aware ones for candidate k;
    neigh: [(key, 'person' | 'called')] in its window"""
    u, span, cp = repeat_at(fa, chrom, k)
    L = abs(len(k[1]) - len(k[2]))
    c = hap_counts(bam, fa, chrom, k, span, [x[0] for x in neigh])
    n = c["carry"] + c["without"]
    return ([u, span, cp, L / u if u else 0.0, float(L % u == 0) if u else 0.0] + read_feats(bam, fa, chrom, k, span)
            + [c["carry"], c["without"], c["other"], c["carry"] / n if n else np.nan, len(neigh),
               sum(x[1] == "person" for x in neigh), sum(x[1] == "called" for x in neigh)])


def neighbours(keys, called, person):
    """sorted [(key, kind)] of the person's records and the called candidates, with their positions"""
    nb = sorted([(k, "person") for k in person] + [(k, "called") for k in called])
    return nb, [x[0][0] for x in nb]


def window_neigh(nb, pos, k, span):
    lo, hi = k[0] - WIN, k[0] + len(k[1]) + span + WIN
    out, j = [], bisect.bisect_left(pos, lo)
    while j < len(nb) and nb[j][0][0] <= hi:
        if nb[j][0] != k:
            out.append(nb[j])
        j += 1
    return out


def apply(p_hom, keys, X, p_any, thr, person, bam, fa, chrom, model):
    """p_hom with the haplotype-aware model's value for every called indel candidate (p_any >= thr) that is not one
    of the person's records -> (p_hom, n replaced)"""
    called = [k for k, p in zip(keys, p_any) if p >= thr and k not in person]
    nb, pos = neighbours(keys, called, person)
    out, n = p_hom.copy(), 0
    for i, k in enumerate(keys):
        if p_any[i] < thr or k in person or len(k[1]) == len(k[2]):
            continue
        _, span, _ = repeat_at(fa, chrom, k)
        f = feats(bam, fa, chrom, k, window_neigh(nb, pos, k, span))
        out[i] = model.predict_proba(np.hstack([X[i], f])[None, :])[0, 1]
        n += 1
    return out, n
