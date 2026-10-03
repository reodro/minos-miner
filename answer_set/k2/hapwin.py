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


def candidates(windows, answer, p_of, donor, bed):
    """Windows inside the confident region whose called pair differs from the answer's, with their features.
    answer: {key: (zyg, synthetic)} (the per-variant answer, the person's records included); p_of: {key: model
    score} of the per-variant candidates; donor: the person's GIAB records {key: zyg}. -> [(s, e, window,
    answer records inside [(key, zyg)], called events {key: zyg}, feature vector)]."""
    ak = sorted(answer)
    ap = [k[0] for k in ak]
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
        out.append((s, e, w, [(k, answer[k][0]) for k in A], ev, np.array([f[n] for n in FEATS], dtype=float)))
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


def apply(recs, cands, p, thr, respell_fn=None):
    """recs (assemble_records output: (pos, ref, alt, zyg, synthetic, confidence)) with the windows scored >= thr
    replaced by their called events (most confident first; overlapping windows skipped). respell_fn(cand) -> events in
    another spelling (or None to keep the genotyper's); windows holding the person's own records keep theirs."""
    order = sorted(range(len(cands)), key=lambda i: -p[i])
    taken, drop, add = [], set(), []
    own = {(r[0], r[1], r[2]) for r in recs if not r[4]}
    for i in order:
        if p[i] < thr:
            break
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
