"""Local haplotype genotyper: which two versions of the local DNA does the person carry? Owner: answer_set.

Inserted indels in repeats are in the reads (2026-10-02, 28 rounds: 91% of the inserted 1-bp indels that never became
candidates have >= 2 anchored reads closer to the truth's local DNA than to the reference; het ~52% of the reads,
hom ~94%), but the aligner spells them differently read by read, and Minos's noise splits one real sequence into many
versions - so neither the column-wise read scan nor counting identical read sequences (hapcand.local_sequences) finds
them.

Per window (hapcand.windows):
  1. reads: each read's sequence between two fixed anchor stretches flanking the window (hapcand._find: one mismatch
     allowed, soft-clipped ends included), with its base qualities - independent of how the aligner wrote it;
  2. candidate haplotypes: the reference, plus local sequences carried by >= 2 reads (exactly, or with low-quality
     bases as wildcards - hapcand.local_sequences), at most MAX_HAPS;
  3. read likelihoods: each read aligned to each haplotype (edlib, global), scored with its base qualities
     (match log(1-e), mismatch log(e/3), indel GAP_OPEN + GAP_EXT per extra base), floored at FLOOR below the read's
     best haplotype so one odd read cannot veto a haplotype;
  4. genotype: every pair of haplotypes (ref/ref, ref/h het, h/h hom, h1/h2), log-likelihood + a prior per non-reference
     haplotype (PRIOR_ALT), posterior over pairs;
  5. events: each non-reference haplotype of a likely pair broken into simple events against the reference
     (_events_of: Minos-style spelling, see SCHEMES; left-aligned like every other candidate).

Per event: P(present) and P(hom) summed over the genotype posterior, and the reads best explained by a haplotype that
carries it. No variant caller is used, so the answer set stays independent of the GATK configs being tuned.

Edit-aware (GIAB, with the donor's Google originals - giab_origin.py): every round read is a piece of an original
read, edited by Minos where it carries an inserted mutation (or Minos noise). Each anchored read is matched to its
original (same strand, position bin, exact substring or best edit distance) and the original's sequence between the
same anchors is compared: "edited" if it differs. A sequence carried by edited reads from several different originals
is Minos's edit; one that came with the originals is the donor's DNA (or the original run's errors). 28 rounds,
2026-10-02: 481/520 inserted 1-bp indels that never became candidates have >= 2 edited reads carrying exactly the
truth's sequence; at 11,200 repeats without a truth record, >= 2 edited reads sharing one 1-bp change: 0.25%.
"""
import math
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import edlib
import parasail
import pysam

from answer_set.k2 import hapcand

MAX_HAPS = 8              # candidate haplotypes per window, the reference included
MIN_HAP_READS = 2
GAP_OPEN, GAP_EXT = math.log(1e-3), math.log(0.1)
FLOOR = 12.0              # a read never counts more than this (natural log) against a haplotype
PRIOR_ALT = -4.0          # log prior per non-reference haplotype in a pair
MIN_EVENT_POST = 0.05     # events of less likely pairs are not reported
EDIT_FEATS = ["hg_edit", "hg_edit_ori", "hg_unedit", "hg_win_edit_other", "hg_edit_q"]
FEATS = ["hg_post", "hg_hom", "hg_reads", "hg_anchored", "hg_lr", "hg_n_haps"] + EDIT_FEATS
# Genotyper version a model is trained with (meta hap_genotyper; build.py runs that one). 1: g009. (Tried 2026-10-03
# and dropped: denoised candidate sequences - truth's sequence absent 64 -> 62 of 275 missed sites, +40 s per round.)
VERSION = 1
BIN = 50                  # originals grouped by (strand, start // BIN); a round read looks in its bin and both neighbours


def _between(q, left, right, max_len):
    L = hapcand._find(left, q)
    if L is None:
        return None
    R = hapcand._find(right, q, L[1] + 1)
    if R is None or R[0] - (L[1] + 1) > max_len:
        return None
    return L[1] + 1, R[0]


def _originals(ob, chrom, s, e):
    """{(is_reverse, start // BIN): [original reads]} around the window."""
    by = defaultdict(list)
    for o in ob.fetch(chrom, max(0, s - 400), e + 400):
        if o.is_unmapped or o.is_secondary or o.is_supplementary or not o.query_sequence:
            continue
        by[(o.is_reverse, o.reference_start // BIN)].append(o)
    return by


ORIGIN_MAX_EDITS = 4      # a round read differs from its original by Minos's edits only (78% exact, 13% one edit)


def _origin(rd, q, by):
    """[(edit distance, original)] of the originals closest to a round read (within ORIGIN_MAX_EDITS), best first."""
    b = rd.reference_start // BIN
    cands = by.get((rd.is_reverse, b - 1), []) + by.get((rd.is_reverse, b), []) + by.get((rd.is_reverse, b + 1), [])
    exact = [(0, o) for o in cands if q in o.query_sequence]
    if exact:
        return exact
    out = []
    for o in cands:
        d = edlib.align(q, o.query_sequence, mode="HW", task="distance", k=ORIGIN_MAX_EDITS)["editDistance"]
        if d >= 0:
            out.append((d, o))
    out.sort(key=lambda t: t[0])
    return [t for t in out if t[0] == out[0][0]] if out else []


def _anchored_reads(bam, fa, chrom, s, e, ob=None):
    """[(sequence, qualities, status, original name, edit quality)] of reads anchored on both sides of [s, e]
    (1-based, inclusive). status: "edit" (differs from its original there), "orig" (as in its original), "?" (no
    original / the original does not span the anchors / no originals given); edit quality: mean base quality of the
    read where it differs from its original (Minos's noise is mostly low quality), NaN unless edited."""
    left = fa.fetch(chrom, s - 1 - hapcand.ANCHOR, s - 1).upper()
    right = fa.fetch(chrom, e, e + hapcand.ANCHOR).upper()
    max_len = (e - s + 1) + 2 * hapcand.MAX_INDEL
    by = _originals(ob, chrom, s, e) if ob is not None else None
    out = []
    for rd in bam.fetch(chrom, max(0, s - 1), e):
        if rd.is_secondary or rd.is_supplementary or rd.is_duplicate or rd.is_unmapped or not rd.query_sequence:
            continue
        q = rd.query_sequence
        span = _between(q, left, right, max_len)
        if span is None:
            continue
        i0, i1 = span
        qual = rd.query_qualities
        status, oname = "?", ""
        if by is not None:
            # the closest originals must agree on the window, and the read may differ from them only inside it -
            # otherwise the read's own original is missing and a look-alike (e.g. the donor's other copy) was found
            locals_ = set()
            best = _origin(rd, q, by)
            for d, o in best:
                os_ = _between(o.query_sequence, left, right, max_len)
                if os_ is None:
                    locals_.add(None)
                    continue
                inner = o.query_sequence[os_[0]:os_[1]]
                d_in = 0 if inner == q[i0:i1] else edlib.align(q[i0:i1], inner, task="distance")["editDistance"]
                locals_.add(inner if d - d_in <= 1 else None)
            if len(locals_) == 1 and None not in locals_:
                inner = locals_.pop()
                status = "orig" if inner == q[i0:i1] else "edit"
                oname = best[0][1].query_name
        rq = list(qual[i0:i1]) if qual is not None else [30] * (i1 - i0)
        eq = math.nan
        if status == "edit" and q[i0:i1] and inner:
            a = edlib.align(q[i0:i1], inner, task="path")
            nice = edlib.getNiceAlignment(a, q[i0:i1], inner)
            qs, qi = [], 0
            for ch, m in zip(nice["query_aligned"], nice["matched_aligned"]):
                if m != "|":
                    qs.append(rq[min(qi, len(rq) - 1)])
                if ch != "-":
                    qi += 1
            eq = sum(qs) / len(qs) if qs else math.nan
        out.append((q[i0:i1], rq, status, oname, eq))
    return out


def _read_ll(seq, qual, hap):
    """log P(read | haplotype): global alignment, base-quality-aware costs."""
    if not seq or not hap:
        return GAP_OPEN + GAP_EXT * max(len(seq), len(hap))
    a = edlib.align(seq, hap, mode="NW", task="path")
    ll, qi = 0.0, 0
    num = ""
    for ch in a["cigar"]:
        if ch.isdigit():
            num += ch
            continue
        n = int(num)
        num = ""
        if ch == "=":
            for k in range(qi, qi + n):
                ll += math.log1p(-10 ** (-qual[k] / 10))
            qi += n
        elif ch == "X":
            for k in range(qi, qi + n):
                ll += math.log(10 ** (-qual[k] / 10) / 3)
            qi += n
        elif ch == "I":                                   # bases in the read, not in the haplotype
            ll += GAP_OPEN + GAP_EXT * (n - 1)
            qi += n
        else:                                             # "D": haplotype bases the read lacks
            ll += GAP_OPEN + GAP_EXT * (n - 1)
    return ll


def genotype_window(bam, fa, chrom, s, e, ob=None):
    """(anchored reads, haplotypes, [(posterior, i, j)] best first, read log-likelihoods) for window [s, e]."""
    ref = fa.fetch(chrom, s - 1, e).upper()
    reads = _anchored_reads(bam, fa, chrom, s, e, ob)
    if len(reads) < hapcand.MIN_ANCHORED:
        return reads, [ref], [], []
    raw = defaultdict(int)
    for q, *_ in reads:
        raw[q] += 1
    masked, _, _ = hapcand.local_sequences(bam, fa, chrom, s, e)
    support = defaultdict(int)
    for x, m in list(raw.items()) + list(masked.items()):
        support[x] = max(support[x], m)
    haps = [ref] + [x for x, m in sorted(support.items(), key=lambda t: -t[1])
                    if x != ref and m >= MIN_HAP_READS and abs(len(x) - len(ref)) <= hapcand.MAX_INDEL][:MAX_HAPS - 1]
    L = []
    for q, qual, *_ in reads:
        row = [_read_ll(q, qual, h) for h in haps]
        best = max(row)
        L.append([max(x, best - FLOOR) for x in row])
    pairs = []
    for i in range(len(haps)):
        for j in range(i, len(haps)):
            ll = sum(math.log(0.5 * math.exp(r[i] - max(r[i], r[j])) + 0.5 * math.exp(r[j] - max(r[i], r[j])))
                     + max(r[i], r[j]) for r in L)
            prior = PRIOR_ALT * len({x for x in (i, j) if x != 0})
            pairs.append((ll + prior, i, j))
    top = max(p[0] for p in pairs)
    z = sum(math.exp(p[0] - top) for p in pairs)
    post = sorted(((math.exp(p[0] - top) / z, p[1], p[2], p[0]) for p in pairs), reverse=True)
    return reads, haps, post, L


# Decomposition: Minos writes its edits as SNPs and indels (<= 10 bp). A global alignment where indels are cheap next
# to mismatches gives its spelling (an insertion + a deletion a few bases apart, not a run of mismatches; two 1-bp
# insertions, not one 2-bp insertion + a SNP): on the test round, 8/11 windows the genotyper fixes come out exactly as
# the truth's records with either scheme, 10/11 with one of them, 0/11 with the prototype's affine decomposition (mismatch 4, gap 6+1).
SCHEMES = ((6, 3, 1), (8, 2, 1))          # (mismatch, gap open, gap extend); match +1
_MATRICES = {mm: parasail.matrix_create("ACGTN", 1, -mm) for mm, _, _ in SCHEMES}


def _decompose(ref, hap, start, mm, go, ge):
    tb = parasail.nw_trace_scan_16(hap, ref, go, ge, _MATRICES[mm]).get_traceback("|", ".", " ")
    qa, ta = tb.query, tb.ref
    if len(qa) != len(ta):
        return None
    ev, ri, i, n = [], 0, 0, len(qa)
    while i < n:
        if qa[i] == ta[i]:
            ri += 1
            i += 1
            continue
        gap_q, gap_t = qa[i] == "-", ta[i] == "-"
        j = i
        while j < n and qa[j] != ta[j] and (qa[j] == "-") == gap_q and (ta[j] == "-") == gap_t:
            j += 1
        o, a = ta[i:j].replace("-", ""), qa[i:j].replace("-", "")
        if len(o) == len(a):
            ev += [(start + ri + k, o[k], a[k]) for k in range(len(o)) if o[k] != a[k]]
        elif ri == 0:
            return None                                    # an indel at the window's first base: no anchor base
        else:
            ev.append((start + ri - 1, ref[ri - 1] + o, ref[ri - 1] + a))
        ri += len(o)
        i = j
    return ev


def _spellings_of(fa, chrom, s, ref, hap, known=frozenset()):
    """The haplotype as simple events, left-aligned: every distinct decomposition of the schemes, preferred first
    (fewest events, then most events already known - other candidates -, then scheme order). Only the preferred one
    becomes candidates: offering the other too (2026-10-03) added ~47,600 candidates per round for ~4 truth records,
    and hap.py, which compares DNA, counts the truth's records found under either spelling."""
    out = {}
    for rank, (mm, go, ge) in enumerate(SCHEMES):
        ev = _decompose(ref, hap, s, mm, go, ge)
        if not ev or any(abs(len(x[1]) - len(x[2])) > hapcand.MAX_INDEL for x in ev):
            continue
        ev = tuple(sorted({hapcand.left_align(fa, chrom, *x) for x in ev}))
        key = (len(ev), -sum(k in known for k in ev), rank)
        if ev not in out or key < out[ev]:
            out[ev] = key
    return [list(ev) for ev, _ in sorted(out.items(), key=lambda t: t[1])]


def _events_of(fa, chrom, s, ref, hap, known=frozenset()):
    """The haplotype's preferred spelling (see _spellings_of), or None."""
    sp = _spellings_of(fa, chrom, s, ref, hap, known)
    return sp[0] if sp else None


def _windows_job(args):
    bam_path, fasta, chrom, wins, known, orig_path = args
    bam, fa = pysam.AlignmentFile(bam_path), pysam.FastaFile(fasta)
    ob = pysam.AlignmentFile(orig_path) if orig_path else None
    out = {}
    for s, e in wins:
        reads, haps, post, L = genotype_window(bam, fa, chrom, s, e, ob)
        if not post:
            continue
        ref = haps[0]
        spell = {h: _spellings_of(fa, chrom, s, ref, haps[h], known)[:1] for h in range(1, len(haps))}
        ev_of = {h: (v[0] if v else None) for h, v in spell.items()}
        best_ll = post[0][3]
        assigned = defaultdict(int)                       # reads whose best haplotype is h
        edit, unedit, edit_ori, edit_q = defaultdict(int), defaultdict(int), defaultdict(set), defaultdict(list)
        for r, (_, _, st, oname, eq) in zip(L, reads):
            m = max(r)
            for h in [h for h, x in enumerate(r) if x == m][:1]:
                assigned[h] += 1
                if st == "edit":
                    edit[h] += 1
                    edit_ori[h].add(oname)
                    if eq == eq:
                        edit_q[h].append(eq)
                elif st == "orig":
                    unedit[h] += 1
        p0, i0, j0, _ = post[0]
        win_edit_other = sum(n for h, n in edit.items() if h not in (i0, j0))
        has_origin = ob is not None
        per = defaultdict(lambda: {"hg_post": 0.0, "hg_hom": 0.0, "hg_reads": 0, "best_without": -math.inf})
        keys_in = set()
        for h, ev in ev_of.items():
            for k in ev or ():
                keys_in.add(k)
        for p, i, j, ll in post:
            carried = {i: set(ev_of.get(i) or ()) if i else set(), j: set(ev_of.get(j) or ()) if j else set()}
            a, b = carried[i], carried[j]
            for k in keys_in:
                if k in a or k in b:
                    per[k]["hg_post"] += p
                    if k in a and k in b:
                        per[k]["hg_hom"] += p
                else:
                    per[k]["best_without"] = max(per[k]["best_without"], ll)
        for k in keys_in:
            if per[k]["hg_post"] < MIN_EVENT_POST:
                continue
            hs = [h for h, ev in ev_of.items() if ev and k in ev]
            reads_k = sum(assigned[h] for h in hs)
            out[k] = {"hg_post": per[k]["hg_post"], "hg_hom": per[k]["hg_hom"], "hg_reads": float(reads_k),
                      "hg_anchored": float(len(reads)), "hg_lr": min(50.0, best_ll - per[k]["best_without"]),
                      "hg_n_haps": float(len(haps) - 1),
                      "hg_edit": float(sum(edit[h] for h in hs)) if has_origin else math.nan,
                      "hg_edit_ori": float(len(set().union(*[edit_ori[h] for h in hs]))) if has_origin else math.nan,
                      "hg_unedit": float(sum(unedit[h] for h in hs)) if has_origin else math.nan,
                      "hg_win_edit_other": float(win_edit_other) if has_origin else math.nan,
                      "hg_edit_q": (sum(x for h in hs for x in edit_q[h]) / max(1, sum(len(edit_q[h]) for h in hs))
                                    if any(edit_q[h] for h in hs) else math.nan)}
        called = []                                       # the best pair's non-reference sequences, for window decisions
        for h in sorted({i0, j0} - {0}):
            called.append({"events": ev_of.get(h) or [], "hom": i0 == j0, "reads": assigned[h],
                           "edit": edit[h] if has_origin else math.nan,
                           "edit_ori": len(edit_ori[h]) if has_origin else math.nan,
                           "unedit": unedit[h] if has_origin else math.nan,
                           "edit_q": sum(edit_q[h]) / len(edit_q[h]) if edit_q[h] else math.nan})
        out[("window", s, e)] = {"hg_anchored": float(len(reads)), "hg_n_haps": float(len(haps) - 1),
                                 "hg_win_edit_other": float(win_edit_other) if has_origin else math.nan,
                                 "ref": ref, "pair": (haps[i0], haps[j0]), "post": p0,
                                 "post2": post[1][0] if len(post) > 1 else 0.0,
                                 "lr2": (post[0][3] - post[1][3]) if len(post) > 1 else 50.0,
                                 "ref_reads": assigned[0], "called": called}
    return out


def genotype(bam_path, fasta, chrom, seeds, threads=8, known=(), orig_bam=None):
    """{event key: features} for the events of likely haplotype pairs, plus {("window", s, e): {...}} per genotyped
    window (so a candidate inside a window that no likely pair carries can be told apart from an unexamined one).
    known: events found otherwise (read scan, read comparison) - preferred when two spellings are equally simple.
    orig_bam: the donor's original reads over the region (giab_origin.fetch); without it the edit features are NaN."""
    wins = hapcand.windows(seeds)
    if not wins:
        return {}
    parts = [wins[i::max(1, threads)] for i in range(max(1, threads))]
    known = frozenset(known)
    out = {}
    with ProcessPoolExecutor(max(1, threads)) as ex:
        for part in ex.map(_windows_job, [(bam_path, fasta, chrom, p, known, str(orig_bam) if orig_bam else None) for p in parts if p]):
            for k, v in part.items():
                if k[0] != "window" and k in out and out[k]["hg_post"] >= v["hg_post"]:
                    continue                              # an event seen in two windows: keep the stronger
                out[k] = v
    return out
