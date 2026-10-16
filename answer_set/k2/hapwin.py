"""Window decisions: where the genotyper's whole-window answer replaces the per-variant answer. Owner: answer_set.

The local haplotype genotyper (hapgeno.py) gets the person's two local sequences right where the per-variant model
misses (28 rounds, 2026-10-03: right in 437 of 619 windows holding an inserted indel the read scan never made a
candidate, the per-variant answer in 104), but handing its sequences to the per-variant model as single events loses
most of it: they are spelled differently from Minos's records, and each event is judged alone among ~1,500 genotyper
events per round of which ~12 are real. Used as a blanket override it breaks far more than it fixes (outside the
confident region, the person's own GIAB records, Google's and Minos's noise).

Here the decision is per window, on DNA: among windows inside the person's confident region whose called pair differs
from the per-variant answer's, a model (trained on released rounds: is the called pair the truth's DNA?) decides
where to replace the answer's records by the called sequences' events. The person's own GIAB records are never
removed: a window qualifies only if the called pair carries every one of them. hap.py compares DNA, so the spelling
of the events does not matter for the score.

emit_from_called: when the called pair is the donor haplotypes plus leftover 1-3 bp Minos events (the original
reads vs the round), and the window has no synthetic records yet, write those leftovers if they are already
candidates with p >= EMIT_MIN_P (below the 1-bp cut, so the per-variant model declined them).
"""
from __future__ import annotations

import gzip
import math
import pickle

import numpy as np

FEATS = ["post", "post2", "lr2", "anchored", "n_haps", "edit_other", "ref_reads",
         "n_called", "called_hom", "reads_min", "reads_max", "edit_min", "edit_max", "edit_ori_min", "edit_ori_max",
         "unedit_max", "edit_q_min", "n_called_events", "n_answer", "kind", "called_len_max", "answer_len_max",
         "called_indel", "answer_indel", "called_p_min", "called_p_max", "called_not_cand", "answer_p_min",
         "answer_p_max", "hp_max", "donor_in", "width"]
KINDS = {"adds": 0, "removes": 1, "changes": 2}
# window models trained with the GATK evidence (meta window_hc; hcrun.py): GATK's calls inside the window - how many,
# whether their DNA pair is the genotyper's called pair or the answer's, their QUAL. A third opinion on the same DNA.
FEATS_HC = ["hc_n", "hc_agree_called", "hc_agree_answer", "hc_qual_min", "hc_qual_max"]
# window models trained with lr_answer (meta window_lr): how much better the called pair explains the window's anchored
# reads than the per-variant answer's records do (log-likelihood, hapgeno's read model; best phasing of the answer's
# het records). Where the genotyper is right and the answer wrong it is large (108 training rounds: median 30.4) - not
# one read may fit the answer (2026-10-11 16:42, chr13:28,078,730: the answer's G>T, read as a substitution in 6 reads
# that carry an insertion + a deletion of equal length) - and smaller where the answer is right (median 13.1). Held-out
# windows, same folds: fixed / broken / net at every cut-off better (best net 17,609 -> 17,987 SNP-error units).
FEATS_LR = ["lr_answer"]


def answer_lr(bam, fa, chrom, s, e, w, A):
    """log-likelihood of the window's anchored reads under the called pair minus under the best phasing of the answer's
    records A [(key, zyg)], capped at +-50; NaN when it cannot be computed."""
    from answer_set.k2 import hapgeno
    try:
        reads = hapgeno._anchored_reads(bam, fa, chrom, s, e, None)
        called = tuple(sorted(w["pair"]))
        aps = pairs_of(w["ref"], s, A) or set()
        if not reads or not aps:
            return math.nan
        haps = sorted({h for pr in aps for h in pr} | set(called))
        L = {h: [hapgeno._read_ll(q, qual, h) for q, qual, *_ in reads] for h in haps}
        best = [max(L[h][i] for h in haps) for i in range(len(reads))]

        def pair_ll(a, b):
            tot = 0.0
            for i in range(len(reads)):
                x, y = max(L[a][i], best[i] - hapgeno.FLOOR), max(L[b][i], best[i] - hapgeno.FLOOR)
                mx = max(x, y)
                tot += mx + math.log(0.5 * math.exp(x - mx) + 0.5 * math.exp(y - mx))
            return tot
        return min(50.0, max(-50.0, pair_ll(*called) - max(pair_ll(*pr) for pr in aps)))
    except Exception:  # noqa: BLE001  (one odd window must not stop the build)
        return math.nan


def save_windows(path, hg):
    """The genotyper's windows (hapgeno.genotype output) to path, for window decisions."""
    wins = [(k[1], k[2], v) for k, v in hg.items() if k[0] == "window"]
    with gzip.open(path, "wb") as f:
        pickle.dump(wins, f)


def load_windows(path):
    with gzip.open(path, "rb") as f:
        return pickle.load(f)


def _apply(seq, s, alle):
    out, i = [], 0
    for p, r, a in sorted(alle):
        j = p - s
        if j < i or seq[j:j + len(r)] != r:
            return None
        out.append(seq[i:j]); out.append(a); i = j + len(r)
    return "".join(out) + seq[i:]


def pairs_of(ref, s, recs):
    """Every (sorted) pair of window sequences the records can give (hom on both copies, het on one); None if too
    many het records to enumerate."""
    hom = [k for k, z in recs if z == "hom"]
    het = [k for k, z in recs if z != "hom"]
    if len(het) > 10:
        return None
    out = set()
    for mask in range(2 ** max(0, len(het) - 1) if het else 1):
        a = hom + [k for b, k in enumerate(het) if not mask >> b & 1]
        b_ = hom + [k for b, k in enumerate(het) if mask >> b & 1]
        ha, hb = _apply(ref, s, a), _apply(ref, s, b_)
        if ha is not None and hb is not None:
            out.add(tuple(sorted((ha, hb))))
    return out


def _inside(sorted_keys, pos, s, e):
    """records fully inside [s, e]; None if one crosses an edge."""
    import bisect
    i0, i1 = bisect.bisect_left(pos, s - 60), bisect.bisect_right(pos, e)
    out = []
    for k in sorted_keys[i0:i1]:
        k0, k1 = k[0], k[0] + len(k[1]) - 1
        if k1 < s or k0 > e:
            continue
        if k0 < s or k1 > e:
            return None
        out.append(k)
    return out


def pair_spellings(fa, chrom, s, ref, hap, max_len=10):
    """Every insertion + deletion pair (each <= max_len bp, the insertion's bases taken from hap) that turns ref
    (window starting at 1-based s) into hap exactly, left-aligned like bcftools norm. Minos often makes an insertion
    and a deletion 10-20 bp apart in a repeat: the DNA between them is only shifted, so an aligner puts the two
    indels elsewhere (2026-10-03: all 42 replaced windows whose spelling no alignment scheme reached were such pairs)."""
    from answer_set.k2 import hapcand
    d = len(hap) - len(ref)
    out = set()
    for a in range(0, max_len + 1):
        b = a - d
        if b < 1 or b > max_len or a < 1:
            continue                                       # single indels: the alignment schemes cover them
        for i in range(1, len(ref)):                       # insertion of hap[i:i+a] before ref[i]
            if hap[:i] != ref[:i]:
                break
            ins = hap[i:i + a]
            for j in range(i, len(ref) - b):               # deletion of ref[j:j+b] after the insertion
                if ref[i:j] + ref[j + b:] == hap[i + a:]:
                    ev = {(s + i - 1, ref[i - 1], ref[i - 1] + ins), (s + j - 1, ref[j - 1:j + b], ref[j - 1])}
                    out.add(tuple(sorted(hapcand.left_align(fa, chrom, *x) for x in ev)))
        for j in range(1, len(ref) - b):                   # deletion first, insertion after it
            if hap[:j] != ref[:j]:
                break
            for i in range(j + b, len(ref) + 1):
                ins = hap[i - b:i - b + a]
                if ref[j + b:i] == hap[j:i - b] and ref[i:] == hap[i - b + a:]:
                    ev = {(s + j - 1, ref[j - 1:j + b], ref[j - 1]), (s + i - 1, ref[i - 1], ref[i - 1] + ins)}
                    out.add(tuple(sorted(hapcand.left_align(fa, chrom, *x) for x in ev)))
    return [list(x) for x in sorted(out)]                 # a fixed order: the spelling chooser breaks ties by position


def _hp_max(seq):
    best = run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def candidates(windows, answer, p_of, donor, bed, hc=None, lr_fn=None):
    """Windows inside the confident region whose called pair differs from the answer's, with their features.
    answer: {key: (zyg, synthetic)} (the per-variant answer, the person's records included); p_of: {key: model
    score} of the per-variant candidates; donor: the person's GIAB records {key: zyg}; hc: the GATK run's calls
    (hcrun.load_calls) for window models trained with them (FEATS + FEATS_HC). -> [(s, e, window,
    answer records inside [(key, zyg)], called events {key: zyg}, feature vector)]."""
    ak = sorted(answer)
    ap = [k[0] for k in ak]
    if hc is not None:
        hk = sorted(hc)
        hp = [k[0] for k in hk]
    out = []
    for s, e, w in windows:
        if "called" not in w or not bed.contains(s) or not bed.contains(e) or bed.covered_bp(s, e) != e - s + 1:
            continue
        A = _inside(ak, ap, s, e)
        if A is None:
            continue
        ref = w["ref"]
        called = tuple(sorted(w["pair"]))
        aps = pairs_of(ref, s, [(k, answer[k][0]) for k in A])
        if aps is None or called in aps:
            continue
        ev = {}
        for c in w["called"]:
            for k in c["events"]:
                ev[k] = "hom" if c["hom"] or k in ev else "het"
        if any(k not in ev for k in A if k in donor):
            continue                                       # never remove the person's own GIAB records
        cl = w["called"]
        nan = math.nan
        f = {"post": w["post"], "post2": w.get("post2", nan), "lr2": w.get("lr2", nan), "anchored": w["hg_anchored"],
             "n_haps": w["hg_n_haps"], "edit_other": w["hg_win_edit_other"], "ref_reads": w.get("ref_reads", nan),
             "n_called": len(cl), "called_hom": float(any(c["hom"] for c in cl)),
             "reads_min": min([c["reads"] for c in cl], default=nan), "reads_max": max([c["reads"] for c in cl], default=nan),
             "edit_min": min([c["edit"] for c in cl], default=nan), "edit_max": max([c["edit"] for c in cl], default=nan),
             "edit_ori_min": min([c["edit_ori"] for c in cl], default=nan),
             "edit_ori_max": max([c["edit_ori"] for c in cl], default=nan),
             "unedit_max": max([c["unedit"] for c in cl], default=nan),
             "edit_q_min": min([c["edit_q"] for c in cl if c["edit_q"] == c["edit_q"]], default=nan),
             "n_called_events": len(ev), "n_answer": len(A),
             "kind": KINDS["adds" if not A else "removes" if not ev else "changes"],
             "called_len_max": max([abs(len(k[1]) - len(k[2])) for k in ev], default=0),
             "answer_len_max": max([abs(len(k[1]) - len(k[2])) for k in A], default=0),
             "called_indel": float(any(len(k[1]) != len(k[2]) for k in ev)),
             "answer_indel": float(any(len(k[1]) != len(k[2]) for k in A)),
             "called_p_min": min([p_of[k] for k in ev if k in p_of], default=nan),
             "called_p_max": max([p_of[k] for k in ev if k in p_of], default=nan),
             "called_not_cand": float(sum(k not in p_of for k in ev)),
             "answer_p_min": min([p_of.get(k, 2.0) for k in A], default=nan),
             "answer_p_max": max([p_of.get(k, 2.0) for k in A], default=nan),
             "hp_max": _hp_max(ref), "donor_in": float(sum(k in donor for k in A)), "width": e - s + 1}
        names = FEATS
        if hc is not None:
            names = FEATS + FEATS_HC
            H = _inside(hk, hp, s, e)
            hps = pairs_of(ref, s, [(k, "hom" if hc[k]["hc_hom"] else "het") for k in H]) if H is not None else None
            q = [hc[k]["hc_qual"] for k in H] if H else []
            f.update(hc_n=nan if H is None else float(len(H)),
                     hc_agree_called=nan if hps is None else float(called in hps),
                     hc_agree_answer=nan if hps is None else float(bool(hps & aps)),
                     hc_qual_min=min(q, default=nan), hc_qual_max=max(q, default=nan))
        if lr_fn is not None:                              # window models trained with lr_answer (FEATS_LR)
            names = names + FEATS_LR
            f["lr_answer"] = lr_fn(s, e, w, [(k, answer[k][0]) for k in A])
        out.append((s, e, w, [(k, answer[k][0]) for k in A], ev, np.array([f[n] for n in names], dtype=float)))
    return out


def label(cands, truth):
    """Per candidate window: 1 if the called pair is the truth's DNA, and whether the answer's pair was right.
    truth: {key: (zyg, synthetic)} (released truth)."""
    tk = sorted(truth)
    tp = [k[0] for k in tk]
    y, answer_ok = [], []
    for s, e, w, A, ev, _ in cands:
        T = _inside(tk, tp, s, e)
        ref = w["ref"]
        tps = pairs_of(ref, s, [(k, truth[k][0]) for k in T]) if T is not None else None
        if not tps:
            y.append(0); answer_ok.append(True)            # truth record crossing the edge: leave the window alone
            continue
        y.append(int(tuple(sorted(w["pair"])) in tps))
        aps = pairs_of(ref, s, A)
        answer_ok.append(bool(aps and aps & tps))
    return np.array(y, dtype=np.int8), answer_ok


# Spelling of replaced windows. The validator's scorer gives each truth record a class (SNP, 1-bp indel, ...) from
# its spelling, so where our answer set stands in for the truth, a window written as one SNP where Minos wrote a
# deletion + an insertion is scored in the wrong class. Options: the decompositions of several alignment schemes and
# every insertion + deletion pair (pair_spellings); a model trained on released rounds (train_spelling.py) picks the
# one Minos most likely wrote. 2026-10-03, windows the window decision fixes, 108 rounds held out: exactly Minos's
# records 59.9% (fewest events) -> 72.0% (learned); reachable 83.3%.
SPELL_SCHEMES = [(6, 3, 1), (8, 2, 1), (4, 6, 1), (12, 1, 1), (20, 1, 1), (10, 2, 0), (5, 4, 1), (16, 3, 1)]
SPELL_FEATS = ["n_ev", "n_indel", "n_snp", "n_known", "p_sum", "p_max", "p_min", "scheme", "n_ins", "n_del",
               "len_max", "span", "del_ins_pair", "unknown_indel", "n_options", "rank_fewest"]


def _clean(ev, s, ref, hap):
    """A usable spelling: events that do not overlap one another after left-alignment and give exactly hap.
    Overlapping records (an insertion inside a deletion, two events at one position) are the same DNA but hap.py
    cannot match them (2026-10-03: 4 of 28 rounds scored up to 0.05 lower with them)."""
    ev = sorted(ev)
    for x, y in zip(ev, ev[1:]):
        if y[0] <= x[0] + len(x[1]) - 1:
            return False
    return _apply(ref, s, ev) == hap


def spelling_options(fa, chrom, s, ref, hap):
    """[(events, scheme index; 99 = insertion + deletion pair)] - every distinct usable spelling of hap."""
    import parasail
    from answer_set.k2 import hapcand, hapgeno
    out, seen = [], set()
    for i, (mm, go, ge) in enumerate(SPELL_SCHEMES):
        hapgeno._MATRICES.setdefault(mm, parasail.matrix_create("ACGTN", 1, -mm))
        ev = hapgeno._decompose(ref, hap, s, mm, go, ge)
        if not ev or any(abs(len(x[1]) - len(x[2])) > hapcand.MAX_INDEL for x in ev):
            continue
        ev = tuple(sorted({hapcand.left_align(fa, chrom, *x) for x in ev}))
        if ev not in seen and _clean(ev, s, ref, hap):
            seen.add(ev); out.append((ev, i))
    for ev in pair_spellings(fa, chrom, s, ref, hap):
        ev = tuple(ev)
        if ev not in seen and _clean(ev, s, ref, hap):
            seen.add(ev); out.append((ev, 99))
    return out


def spelling_features(opts, p_of):
    """Feature rows (SPELL_FEATS) for the options of one sequence."""
    order = sorted(range(len(opts)), key=lambda j: (len(opts[j][0]), opts[j][1]))
    rows = []
    for j, (ev, i) in enumerate(opts):
        ind = [k for k in ev if len(k[1]) != len(k[2])]
        ps = [p_of[k] for k in ev if k in p_of]
        ins = [k for k in ind if len(k[2]) > len(k[1])]
        dels = [k for k in ind if len(k[1]) > len(k[2])]
        rows.append([len(ev), len(ind), len(ev) - len(ind), len(ps), sum(ps), max(ps, default=math.nan),
                     min(ps, default=math.nan), i, len(ins), len(dels),
                     max([abs(len(k[1]) - len(k[2])) for k in ind], default=0), max(k[0] for k in ev) - min(k[0] for k in ev),
                     float(any(abs(a[0] - b[0]) <= 20 for a in ins for b in dels)), float(sum(k not in p_of for k in ind)),
                     len(opts), order.index(j)])
    return np.array(rows, dtype=float)


def respell(fa, chrom, s, w, p_of, chooser):
    """The window's called pair as events in the spelling the chooser prefers {key: zyg}; None if a sequence has no
    usable spelling."""
    ref = w["ref"]
    picked = []
    for h in sorted(set(w["pair"]) - {ref}):
        opts = spelling_options(fa, chrom, s, ref, h)
        if not opts:
            return None
        q = chooser.predict_proba(spelling_features(opts, p_of))[:, 1]
        picked.append(set(opts[int(np.argmax(q))][0]))
    hom = w["pair"][0] == w["pair"][1]
    out = {}
    for ev in picked:
        for k in ev:
            out[k] = "hom" if hom or k in out else "het"
    return out


def window_group(c):
    """what a candidate window changes: "1-3 bp indel" (any), "SNP only", or "longer indel" - for per-group cut-offs"""
    L = [abs(len(k[1]) - len(k[2])) for k in list(c[4]) + [k for k, _ in c[3]]]
    return "1-3 bp indel" if any(1 <= x <= 3 for x in L) else "SNP only" if all(x == 0 for x in L) else "longer indel"


def select_key_changes(cands, p, thr, answer):
    """Windows at or above thr whose synthetic alleles differ from the answer.

    A window that only flips het/hom on the same alleles is not returned: on the rounds
    where that flip was checked, the per-variant genotype was the true one. answer maps
    a key to (zyg, synthetic)."""
    donor_keys = {k for k, (z, syn) in answer.items() if not syn}
    t_of = (lambda i: thr.get(window_group(cands[i]), thr.get("default", 0.5))) if isinstance(thr, dict) \
        else (lambda i: thr)
    keep, pk = [], []
    for i, c in enumerate(cands):
        if p[i] < t_of(i):
            continue
        syn_a = {k for k, _z in c[3] if k not in donor_keys}
        syn_e = {k for k in c[4] if k not in donor_keys}
        if syn_a != syn_e:
            keep.append(c)
            pk.append(float(p[i]))
    return keep, np.array(pk, dtype=float)


def apply(recs, cands, p, thr, respell_fn=None):
    """recs (assemble_records output: (pos, ref, alt, zyg, synthetic, confidence)) with the windows scored >= thr
    replaced by their called events (most confident first; overlapping windows skipped). thr: one cut-off, or
    {window_group: cut-off} (meta window_thr_group; "default" for a group not listed). respell_fn(cand) -> events in
    another spelling (or None to keep the genotyper's); windows holding the person's own records keep theirs."""
    order = sorted(range(len(cands)), key=lambda i: -p[i])
    t_of = (lambda i: thr.get(window_group(cands[i]), thr.get("default", 0.5))) if isinstance(thr, dict) else (lambda i: thr)
    taken, drop, add = [], set(), []
    own = {(r[0], r[1], r[2]) for r in recs if not r[4]}
    for i in order:
        if p[i] < t_of(i):
            continue
        s, e, w, A, ev, _ = cands[i]
        if any(not (e < a or s > b) for a, b in taken):
            continue
        taken.append((s, e))
        drop.update(k for k, _ in A)
        if respell_fn is not None and not any(k in own for k in ev):
            ev = respell_fn(cands[i]) or ev
        add.extend((k[0], k[1], k[2], z, True, float(p[i])) for k, z in ev.items())
    # the person's own GIAB records (synthetic False) stay exactly as published, het/hom included
    kept = [r for r in recs if (r[0], r[1], r[2]) not in drop or not r[4]]
    have = {(r[0], r[1], r[2]) for r in kept}
    return kept + [r for r in add if (r[0], r[1], r[2]) not in have], len(taken)


def _donor_assignments(donor_recs):
    """Each way to place het donor records on one haplotype (hom always on). None if too many hets."""
    hom = [k for k, z in donor_recs if z == "hom"]
    het = [k for k, z in donor_recs if z != "hom"]
    if len(het) > 8:
        return None
    return [hom + [k for i, k in enumerate(het) if mask >> i & 1] for mask in range(2 ** len(het))]


def _applied_gpos(seq, s, alle):
    """1-based genomic pos of each base in _apply(seq, s, alle) (inserted bases keep the event POS)."""
    out, i = [], 0
    for p, r, a in sorted(alle):
        j = p - s
        if j < i or seq[j:j + len(r)] != r:
            return None
        out.extend(s + k for k in range(i, j))
        for t in range(len(a)):
            out.append(p + t if t < len(r) else p + len(r) - 1)
        i = j + len(r)
    out.extend(s + k for k in range(i, len(seq)))
    return out


def _events_src_to_dst(fa, chrom, src, dst, gpos):
    """Simple events turning src into dst, positions mapped from src indices to genomic via gpos."""
    import parasail
    from answer_set.k2 import hapcand, hapgeno
    if src == dst:
        return []
    best = None
    for mm, go, ge in hapgeno.SCHEMES:
        hapgeno._MATRICES.setdefault(mm, parasail.matrix_create("ACGTN", 1, -mm))
        raw = hapgeno._decompose(src, dst, 1, mm, go, ge)
        if not raw or any(abs(len(x[1]) - len(x[2])) > hapcand.MAX_INDEL for x in raw):
            continue
        ev = []
        ok = True
        for p, r, a in raw:
            if p < 1 or p > len(gpos):
                ok = False
                break
            ev.append(hapcand.left_align(fa, chrom, gpos[p - 1], r, a))
        if not ok:
            continue
        key = (len(ev),)
        if best is None or key < best[0]:
            best = (key, ev)
    return None if best is None else best[1]


def _leftover_on_hap(fa, chrom, s, ref, hap, used):
    """Minos events on top of donor assignment `used` that produce hap, or None."""
    base = _apply(ref, s, used)
    if base is None:
        return None
    if base == hap:
        return []
    gpos = _applied_gpos(ref, s, used)
    if gpos is None or len(gpos) != len(base):
        return None
    extra = _events_src_to_dst(fa, chrom, base, hap, gpos)
    if extra is None:
        return None
    merged = list(used) + list(extra)
    if _apply(ref, s, merged) != hap:
        return None
    if not _clean(tuple(sorted(set(merged))), s, ref, hap):
        return None
    return extra


def spell_called(fa, chrom, s, w, donor_recs):
    """Called pair as {key: zyg}, donor records kept as their own events. None if a hap cannot be spelled that way.

    Pair each called haplotype with a donor haplotype (the originals) and take the leftover Minos events;
    a homopolymer can look like the donor insertion on both copies, so assignments are chosen so the two
    leftovers agree as well as they can, then checked to reproduce the called pair.
    """
    ref = w["ref"]
    asg = _donor_assignments(donor_recs)
    if asg is None:
        return None
    used_of = {}
    for used in asg:
        seq = _apply(ref, s, used)
        if seq is None:
            return None
        used_of[seq] = used
    dpair = pairs_of(ref, s, donor_recs) if donor_recs else {(ref, ref)}
    if not dpair or len(dpair) != 1:
        return None
    d1, d2 = next(iter(dpair))
    c1, c2 = tuple(sorted(w["pair"]))
    best = None
    for x, y in ((d1, d2), (d2, d1)):
        if x not in used_of or y not in used_of:
            continue
        e1 = _leftover_on_hap(fa, chrom, s, ref, c1, used_of[x])
        e2 = _leftover_on_hap(fa, chrom, s, ref, c2, used_of[y])
        if e1 is None or e2 is None:
            continue
        s1, s2 = set(e1), set(e2)
        key = (len(s1 | s2), len(s1 ^ s2))
        if best is None or key < best[0]:
            best = (key, s1, s2)
    if best is None:
        return None
    _, a, b = best
    out = {k: z for k, z in donor_recs}
    for k in a | b:
        if k not in out:
            out[k] = "hom" if k in a and k in b else "het"
    got = pairs_of(ref, s, list(out.items()))
    if not got or tuple(sorted(w["pair"])) not in got:
        return None
    return out


def _short_indel_extras(ev, donor_keys):
    """Leftover events worth writing: 1-3 bp indels, no SNPs, disjoint from the person's records."""
    extra = [k for k in ev if k not in donor_keys]
    if not extra or len(extra) > 3:
        return None
    if any(k in donor_keys for k in extra):
        return None
    L = [abs(len(k[1]) - len(k[2])) for k in extra]
    if any(x == 0 or x > 3 for x in L):
        return None
    return extra


EMIT_MIN_P = 0.3                 # leftover must already be a per-variant candidate this strong (below the 1-bp cut)


def emit_from_called(recs, windows, donor, bed, fa, chrom, p_of=None):
    """Where the genotyper's pair is the donor haplotypes plus leftover 1-3 bp Minos events, the
    window has no synthetic records yet, and each leftover is already a candidate with p >= EMIT_MIN_P,
    write those leftovers (original → round). The person's GIAB records are never removed.
    -> (records, n_respelled, n_added)."""
    p_of = p_of or {}
    ak = sorted((r[0], r[1], r[2]) for r in recs)
    ap = [k[0] for k in ak]
    by = {(r[0], r[1], r[2]): r for r in recs}
    taken, n_sp, n_add = [], 0, 0
    add, drop = [], set()
    order = sorted((w for w in windows if "pair" in w[2] and "ref" in w[2]), key=lambda t: -t[2].get("post", 0.0))
    for s, e, w in order:
        if not bed.contains(s) or not bed.contains(e) or bed.covered_bp(s, e) != e - s + 1:
            continue
        if any(not (e < a or s > b) for a, b in taken):
            continue
        if w.get("post", 0.0) < 0.95:
            continue
        if w["pair"][0] == w["ref"] and w["pair"][1] == w["ref"]:
            continue
        dlen = [abs(len(h) - len(w["ref"])) for h in w["pair"]]
        has_indel = any(len(k[1]) != len(k[2]) for c in w.get("called") or () for k in c.get("events") or ())
        if max(dlen) == 0 and not has_indel:
            continue
        A = _inside(ak, ap, s, e)
        if A is None:
            continue
        donor_recs = [(k, donor[k]) for k in A if k in donor]
        ev = spell_called(fa, chrom, s, w, donor_recs)
        if ev is None:
            continue
        if any(k not in ev for k, _ in donor_recs):
            continue
        ref = w["ref"]
        called = tuple(sorted(w["pair"]))
        aps = pairs_of(ref, s, [(k, by[k][3]) for k in A])
        eps = pairs_of(ref, s, list(ev.items()))
        if not eps or called not in eps:
            continue
        extra = _short_indel_extras(ev, set(donor))
        if extra:
            extra = [k for k in extra if p_of.get(k, 0.0) >= EMIT_MIN_P]
        syn_in = [k for k in A if k not in donor]
        # only ADD into windows that have no synthetic records yet (donor / empty).
        if syn_in or not extra or (aps and called in aps):
            continue
        taken.append((s, e))
        hom_pos = {k[0] for k in extra if ev[k] == "hom"}
        for k in A:
            if k not in donor and k[0] in hom_pos:
                drop.add(k)
        for k in extra:
            if k in donor or (k in by and not by[k][4]):
                continue
            add.append((k[0], k[1], k[2], ev[k], True, float(w["post"])))
        n_add += 1
    kept = [r for r in recs if (r[0], r[1], r[2]) not in drop]
    have = {(r[0], r[1], r[2]) for r in kept}
    return kept + [r for r in add if (r[0], r[1], r[2]) not in have], n_sp, n_add


LONG_IMF = 0.20          # a true long indel the window left out still had at least 23% of the reads (g075 rounds)


def drop_unneeded(recs, windows, imf=None):
    """Drop a synthetic indel whose removal leaves a confident window's two haplotypes unchanged.

    Two synthetic alleles at one position: the one the haplotypes do not use is a second spelling.
    A single indel of 8 bp or more is the same, when fewer than one read in five carries it.
    Exactly one such allele per position is removed. A donor record stays, and a window that
    already matches is left alone."""
    from answer_set.k2.util import trim
    imf = {trim(*k): v for k, v in (imf or {}).items()}
    by = {(r[0], r[1], r[2]): r for r in recs}
    ak = sorted(by)
    ap = [k[0] for k in ak]
    syn_at = {}
    for r in recs:
        if r[4]:
            syn_at.setdefault(r[0], []).append((r[0], r[1], r[2]))
    named = []
    for s, e, w in windows:
        if w.get("post", 0.0) < 0.98:
            continue
        if w["pair"][0] == w["ref"] and w["pair"][1] == w["ref"]:
            continue
        A = _inside(ak, ap, s, e)
        if not A:
            continue
        inside = [by[k] for k in A]
        if sum(1 for r in inside if r[3] != "hom") > 8:
            continue
        if not any(r[4] and len(r[1]) != len(r[2]) and (
                len(syn_at.get(r[0], ())) >= 2
                or abs(len(r[1]) - len(r[2])) >= 8 and imf.get(trim(r[0], r[1], r[2]), 1.0) < LONG_IMF)
                for r in inside):
            continue
        called = tuple(sorted(w["pair"]))
        cur = pairs_of(w["ref"], s, [((r[0], r[1], r[2]), r[3]) for r in inside])
        if not cur or called in cur:
            continue
        cands = []
        for r in inside:
            if not r[4] or len(r[1]) == len(r[2]):
                continue
            key = (r[0], r[1], r[2])
            sib = len(syn_at.get(r[0], ())) >= 2
            long = abs(len(r[1]) - len(r[2])) >= 8 and imf.get(trim(*key), 1.0) < LONG_IMF
            if not sib and not long:
                continue
            rest = [x for x in inside if (x[0], x[1], x[2]) != key]
            if sib and not any(x[4] and x[0] == r[0] and (x[0], x[1], x[2]) != key for x in rest):
                if not long:
                    continue
            if sum(1 for x in rest if x[3] != "hom") > 8:
                continue
            pr = pairs_of(w["ref"], s, [((x[0], x[1], x[2]), x[3]) for x in rest])
            if pr and called in pr:
                cands.append(key)
        if len(cands) == 1:
            named.append(cands[0])
    drop = set()
    by_pos = {}
    for k in named:
        by_pos.setdefault(k[0], set()).add(k)
    for keys in by_pos.values():
        if len(keys) == 1:
            drop.add(next(iter(keys)))
    if not drop:
        return recs
    return [r for r in recs if (r[0], r[1], r[2]) not in drop]
