"""Read guard on records the per-variant model did not call (GIAB). Owner: answer_set.

The steps after the per-variant model (cluster decoder, window replacement, emit, run-edge spelling) write inserted
indels the model rejected or never scored; in the monitored g075 rounds nearly half of those were false (106 of 231).
Each such indel is checked against the reads: every read spanning the site is judged by its DNA (its sequence over the
window vs every version of the window OUR records give, the indel taken as het so both versions exist; <= 1 mismatch)
as carrying it or not. Too few carriers -> the record is dropped:
  >= MIN_READS judged reads and < MAX_SHARE of them carry it.
Rule fixed before measuring (2026-10-16); 42 monitored g075 rounds (frozen set): see test/answer_set/guard_eval.py.
On only for a model whose meta sets read_guard; build.py applies it to the final records.
"""
from collections import Counter

FL = 16
MIN_READS, MAX_SHARE = 6, 0.15
MAX_HET = 8                      # windows with more het records than this are left alone


def counts(bam, fa, chrom, key, recs):
    """(carrying, reference) reads at key, each spanning read's DNA vs the versions recs {key: het/hom} give"""
    p, ref, _ = key
    a, b = p - FL, p + len(ref) + FL
    refseq = fa.fetch(chrom, a - 1, b).upper()
    near = sorted(k for k in recs if a <= k[0] and k[0] + len(k[1]) - 1 <= b)
    hom = [k for k in near if recs[k] == "hom"]; het = [k for k in near if recs[k] != "hom"]
    if len(het) > MAX_HET:
        return None

    def applied(evs):
        out, i0 = [], 0
        for pp, rr, xx in sorted(evs):
            j = pp - a
            if j < i0 or refseq[j:j + len(rr)] != rr:
                return None
            out.append(refseq[i0:j]); out.append(xx); i0 = j + len(rr)
        return "".join(out) + refseq[i0:]

    versions = {}
    for mask in range(2 ** len(het)):
        sub = hom + [k for n, k in enumerate(het) if mask >> n & 1]
        h = applied(sub)
        if h is not None:
            versions[h] = versions.get(h, False) or key in sub
    if not versions:
        return None
    c = Counter()
    for r in bam.fetch(chrom, a - 1, b):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.query_sequence:
            continue
        if not (r.reference_start + 1 <= a and r.reference_end >= b):
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
        hit = [has for h, has in versions.items() if len(h) == len(s) and sum(x != y for x, y in zip(h, s)) <= 1]
        if hit:
            c["carries" if any(hit) else "ref"] += 1
    return c["carries"], c["ref"]


def apply(recs, bam, fa, chrom, p_of, thr):
    """records (pos, ref, alt, zyg, synthetic, confidence) without the later steps' inserted indels the reads do not
    carry -> (records, n dropped). p_of {key: p}: the per-variant score after the class cut-offs; a record the model
    called has p >= thr."""
    zyg = {(r[0], r[1], r[2]): r[3] for r in recs}
    out, n = [], 0
    for r in recs:
        k = (r[0], r[1], r[2])
        if not r[4] or len(r[1]) == len(r[2]) or p_of.get(k, 0.0) >= thr:
            out.append(r); continue
        c = counts(bam, fa, chrom, k, {**zyg, k: "het"})
        if c and sum(c) >= MIN_READS and c[0] / sum(c) < MAX_SHARE:
            n += 1
            continue
        out.append(r)
    return out, n
