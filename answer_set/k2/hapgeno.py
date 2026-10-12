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
import bisect
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
# 2: repeat-aware window edges (_repeat_aware). A window edge inside a tandem repeat puts its 12-bp anchor in the repeat,
# where reads that start or end in the repeat cannot be anchored and the anchor may match at a shifted place.
# 28 development rounds (2026-10-10), clustered inserted 1-bp indels g026 missed, windows seeded from the neighbouring
# records: truth's DNA called 59 -> 69 of 112; 120 clustered sites g026 found: 70 -> 70.
REPEAT_MAX_EXTEND = 30    # bases an edge may move out of a repeat
# 3: tight windows in repeats (_tight_windows). A window is padded by hapcand.PAD (15) on each side and its 12-bp anchors
# lie beyond that, so a read must span the repeat + ~27 bp on each side; the read comparison needs ~8. In long repeats
# few reads qualify (g042, 23 development rounds: windows that called another sequence than the truth's anchored 6 reads
# (median), 8.5 where right) - 34% of the missed 1-3 bp indels had >= 2 reads showing the truth's edit in the read
# comparison while the window called something else. Version 3: where the seeds touch a tandem repeat, the window is the
# seeds + the repeats they touch + TIGHT_PAD, so the anchors sit next to the repeat.
# 4: version 3 + a called sequence that no edited read's ORIGINAL confirms (_confirmed == 0: carried only by reads with
# no original found / unedited reads) is dropped and the window re-decided, unless the call is clear-cut (log-likelihood
# margin >= V4_MARGIN). Version 3's extra reads also brought noise (repeats with no truth record within 50 bp, 112
# sampled: a non-reference pair called 12 -> 28 times); development rounds: the truth's sequences v3 called 190 -> 183
# kept, false non-reference sequences 29 -> 11.
# 5: version 2's windows, each with fallbacks (version 3's tight window, then version 1's) tried when it has fewer than
# hapcand.MIN_ANCHORED anchored reads; and every candidate inside a window judged by DNA, whatever its spelling (_carries:
# a sequence carries an event when the event lies on an optimal affine-gap alignment path from the reference to it), plus the 1-bp
# indels a likely sequence carries (_indels_carried) as new events. Versions 1-4 credit an event only when the
# genotyper's own spelling of a sequence contains it. 11 live g042 rounds (2026-10-11), 85 missed 1-bp indels: 23 called
# with the wrong genotype, 22 scored too low, 17 dropped at entry - most with the reads' DNA clearly carrying them
# (median 86% of spanning reads where the read comparison counted 43%); 15 in a version-2 window with < 4 anchored reads
# (v1 4, v3 3), both latest-round misses among them.
TIGHT_PAD = 3
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


V4_MARGIN = 8.0           # version 4: a called sequence no original confirms needs this log-likelihood margin


def _confirmed(reads, x):
    """Edited reads carrying exactly x, cut from distinct originals."""
    return len({oname for q, _, st, oname, _ in reads if q == x and st == "edit"})


def genotype_window(bam, fa, chrom, s, e, ob=None, version=1):
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
    post = _pair_posterior(haps, L)
    if version == 4 and ob is not None and any(st in ("edit", "orig") for _, _, st, _, _ in reads):
        for _ in range(3):           # drop a called sequence no original confirms, unless the call is clear-cut
            margin = (post[0][3] - post[1][3]) if len(post) > 1 else 50.0
            weak = [h for h in {post[0][1], post[0][2]} if h != 0 and _confirmed(reads, haps[h]) == 0]
            if not weak or margin >= V4_MARGIN:
                break
            keep = [k for k in range(len(haps)) if k not in weak]
            haps = [haps[k] for k in keep]
            L = [[row[k] for k in keep] for row in L]
            post = _pair_posterior(haps, L)
    return reads, haps, post, L


def _pair_posterior(haps, L):
    """[(posterior, i, j, log-likelihood + prior)] of every haplotype pair, best first."""
    pairs = []
    for i in range(len(haps)):
        for j in range(i, len(haps)):
            ll = sum(math.log(0.5 * math.exp(r[i] - max(r[i], r[j])) + 0.5 * math.exp(r[j] - max(r[i], r[j])))
                     + max(r[i], r[j]) for r in L)
            prior = PRIOR_ALT * len({x for x in (i, j) if x != 0})
            pairs.append((ll + prior, i, j))
    top = max(p[0] for p in pairs)
    z = sum(math.exp(p[0] - top) for p in pairs)
    return sorted(((math.exp(p[0] - top) / z, p[1], p[2], p[0]) for p in pairs), reverse=True)


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


COST_MM, COST_OPEN, COST_EXT = 4, 3, 1   # version 5's DNA judgement: a SNP (4) beats an insertion + deletion (6), which
_COST_MATRIX = parasail.matrix_create("ACGTN", 0, -COST_MM)   # beats two SNPs (8); one long indel beats its pieces


def _ed(a, b):
    """Alignment cost of b against a (COST_*: mismatch 4, gap 3 + 1 per extra base)."""
    if a == b:
        return 0
    return -parasail.nw_scan_16(b, a, COST_OPEN, COST_EXT, _COST_MATRIX).score


def _apply_key(ref, s, k):
    """ref (window from 1-based s) with event k applied, or None if k does not fit inside."""
    j = k[0] - s
    if j < 0 or j + len(k[1]) > len(ref) or ref[j:j + len(k[1])] != k[1]:
        return None
    return ref[:j] + k[2] + ref[j + len(k[1]):]


def _carries(ref, s, hap, d_hap, k):
    """Does sequence hap (cost d_hap from ref, _ed) carry event k? True when k lies on an optimal alignment path from
    the reference to hap - independent of how anyone spells hap."""
    if d_hap == 0:
        return False
    rk = _apply_key(ref, s, k)
    if rk is None:
        return False
    c = _ed(ref, rk)
    return c > 0 and c + _ed(rk, hap) == d_hap


def _indels_carried(fa, chrom, ref, s, hap, d_hap):
    """{left-aligned 1-bp indel hap carries (_carries): its spelling inside the window}, one per resulting sequence."""
    out, seen = {}, set()
    for i in range(1, len(ref)):
        edits = [(s + i - 1, ref[i - 1:i + 1], ref[i - 1])]                       # delete ref[i]
        edits += [(s + i - 1, ref[i - 1], ref[i - 1] + b) for b in "ACGT"]         # insert b after ref[i-1]
        for k in edits:
            rk = _apply_key(ref, s, k)
            if rk is None or rk in seen:
                continue
            if _ed(ref, rk) + _ed(rk, hap) == d_hap:
                seen.add(rk)
                out.setdefault(hapcand.left_align(fa, chrom, *k), k)
    return out


def _windows_job(args):
    bam_path, fasta, chrom, wins, known, orig_path, version = args
    bam, fa = pysam.AlignmentFile(bam_path), pysam.FastaFile(fasta)
    ob = pysam.AlignmentFile(orig_path) if orig_path else None
    known_sorted = sorted(known) if version >= 5 else []
    out = {}
    for alts in wins:
        if isinstance(alts[0], int):
            alts = (alts,)
        for s, e in alts:                                 # version 5: fallbacks while too few reads are anchored
            reads, haps, post, L = genotype_window(bam, fa, chrom, s, e, ob, version)
            if post:
                break
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
        car_of = ev_of
        if version >= 5:                                  # every event in the window, judged by DNA
            d_of = {h: _ed(ref, haps[h]) for h in range(1, len(haps))}
            local = {k: k for k in keys_in}               # event -> its spelling inside the window (judged)
            likely = {h for p_, i_, j_, _ in post if p_ >= MIN_EVENT_POST for h in (i_, j_) if h}
            for h in likely:
                for k, kw in _indels_carried(fa, chrom, ref, s, haps[h], d_of[h]).items():
                    local.setdefault(k, kw)
            jk = bisect.bisect_left(known_sorted, (s,))
            while jk < len(known_sorted) and known_sorted[jk][0] <= e:
                k = known_sorted[jk]
                if k[0] + len(k[1]) - 1 <= e:
                    local.setdefault(k, k)
                jk += 1
            keys_in = set(local)
            car_of = {h: {k for k, kw in local.items() if _carries(ref, s, haps[h], d_of[h], kw)}
                      for h in range(1, len(haps))}
        for p, i, j, ll in post:
            carried = {i: set(car_of.get(i) or ()) if i else set(), j: set(car_of.get(j) or ()) if j else set()}
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
            hs = [h for h, ev in car_of.items() if ev and k in ev]
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


def _repeat_span(seq, i):
    """[a, b) of the longest tandem repeat (period 1-6, >= 2 copies, span >= 6) covering seq[i], or None."""
    best = None
    for per in range(1, 7):
        for o in range(per):
            st = i - o
            if st < 0 or st + per > len(seq):
                continue
            unit = seq[st:st + per]
            a, b = st, st + per
            while a - per >= 0 and seq[a - per:a] == unit:
                a -= per
            while b + per <= len(seq) and seq[b:b + per] == unit:
                b += per
            while b < len(seq) and seq[b] == unit[(b - a) % per]:
                b += 1
            if b - a >= max(6, 2 * per) and (best is None or b - a > best[1] - best[0]):
                best = (a, b)
    return best


def _repeat_aware(wins, fasta, chrom):
    """Windows whose edges (with their anchors) touch a tandem repeat moved out of it, by at most REPEAT_MAX_EXTEND."""
    fa = pysam.FastaFile(fasta)
    out = []
    for s, e in wins:
        lo = max(1, s - hapcand.ANCHOR - REPEAT_MAX_EXTEND - 10)
        seq = fa.fetch(chrom, lo - 1, e + hapcand.ANCHOR + REPEAT_MAX_EXTEND + 10).upper()
        ns = s
        for x in range(s - hapcand.ANCHOR - 1, s + 1):          # the left anchor and the window's first base
            sp = _repeat_span(seq, x - lo)
            if sp:
                ns = min(ns, sp[0] + lo - 1)
        ne = e
        for x in range(e, e + hapcand.ANCHOR + 2):              # the window's last base and the right anchor
            sp = _repeat_span(seq, x - lo)
            if sp:
                ne = max(ne, sp[1] + lo + 1)
        out.append((max(ns, s - REPEAT_MAX_EXTEND), min(ne, e + REPEAT_MAX_EXTEND)))
    return out


def _tight_windows(seeds, fasta, chrom, cluster_gap=10):
    """Version 3 windows: seeds merged as hapcand.windows does; a merged span touching a tandem repeat becomes the span +
    every repeat it touches + TIGHT_PAD (at most hapcand.MAX_WINDOW wide); spans touching no repeat keep the version-2
    window."""
    fa = pysam.FastaFile(fasta)
    merged = []
    for s, e in sorted(seeds):
        if merged and s - merged[-1][1] <= cluster_gap and max(e, merged[-1][1]) - merged[-1][0] <= hapcand.MAX_WINDOW - 2 * hapcand.PAD:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    out = []
    for s, e in merged:
        lo = max(1, s - REPEAT_MAX_EXTEND - 10)
        seq = fa.fetch(chrom, lo - 1, e + REPEAT_MAX_EXTEND + 10).upper()
        a, b, hit = s, e, False
        for x in range(s - 1, e + 2):
            sp = _repeat_span(seq, x - lo)
            if sp:
                hit = True
                a, b = min(a, sp[0] + lo), max(b, sp[1] + lo - 1)
        ws, we = a - TIGHT_PAD, b + TIGHT_PAD
        if not hit or we - ws + 1 > hapcand.MAX_WINDOW:
            out.extend(_repeat_aware([(s - hapcand.PAD, e + hapcand.PAD)], fasta, chrom))
        else:
            out.append((ws, we))
    return out


def genotype(bam_path, fasta, chrom, seeds, threads=8, known=(), orig_bam=None, version=1):
    """{event key: features} for the events of likely haplotype pairs, plus {("window", s, e): {...}} per genotyped
    window (so a candidate inside a window that no likely pair carries can be told apart from an unexamined one).
    known: events found otherwise (read scan, read comparison) - preferred when two spellings are equally simple.
    orig_bam: the donor's original reads over the region (giab_origin.fetch); without it the edit features are NaN.
    version: the VERSION a model was trained with (2: repeat-aware window edges; 3: tight windows in repeats; 4: 3 +
    original-confirmed sequences; 5: 2 with fallback windows + events judged by DNA)."""
    wins = hapcand.windows(seeds)
    if version >= 5:                                  # version 2's window, then version 3's, then version 1's
        v1, v2, v3 = wins, _repeat_aware(wins, fasta, chrom), _tight_windows(seeds, fasta, chrom)
        if len(v3) != len(v1):                        # (the same seed clusters; never expected to differ)
            v3 = v1
        wins = [tuple(dict.fromkeys(x)) for x in zip(v2, v3, v1)]
    elif version >= 3:                                # 4: same windows, original-supported sequences only
        wins = _tight_windows(seeds, fasta, chrom)
    elif version >= 2:
        wins = _repeat_aware(wins, fasta, chrom)
    if not wins:
        return {}
    parts = [wins[i::max(1, threads)] for i in range(max(1, threads))]
    known = frozenset(known)
    out = {}
    with ProcessPoolExecutor(max(1, threads)) as ex:
        for part in ex.map(_windows_job, [(bam_path, fasta, chrom, p, known, str(orig_bam) if orig_bam else None, version) for p in parts if p]):
            for k, v in part.items():
                if k[0] != "window" and k in out and out[k]["hg_post"] >= v["hg_post"]:
                    continue                              # an event seen in two windows: keep the stronger
                out[k] = v
    return out
