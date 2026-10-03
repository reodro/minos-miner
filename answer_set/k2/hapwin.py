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


def apply(recs, cands, p, thr):
    """recs (assemble_records output: (pos, ref, alt, zyg, synthetic, confidence)) with the windows scored >= thr
    replaced by their called events (most confident first; overlapping windows skipped)."""
    order = sorted(range(len(cands)), key=lambda i: -p[i])
    taken, drop, add = [], set(), []
    for i in order:
        if p[i] < thr:
            break
        s, e, w, A, ev, _ = cands[i]
        if any(not (e < a or s > b) for a, b in taken):
            continue
        taken.append((s, e))
        drop.update(k for k, _ in A)
        add.extend((k[0], k[1], k[2], z, True, float(p[i])) for k, z in ev.items())
    # the person's own GIAB records (synthetic False) stay exactly as published, het/hom included
    kept = [r for r in recs if (r[0], r[1], r[2]) not in drop or not r[4]]
    have = {(r[0], r[1], r[2]) for r in kept}
    return kept + [r for r in add if (r[0], r[1], r[2]) not in have], len(taken)
