"""Train ONE track (giab or kg) as its own version, with honest held-out testing within the track.

  python -m answer_set.k2.train_track --track giab                          # every prepared GIAB round
  python -m answer_set.k2.train_track --track giab --exclude-rounds test.txt  # keep a test set out
  python -m answer_set.k2.train_track --track giab --same-rounds-as g006_20261001T1057
  python -m answer_set.k2.train_track --track giab --promote
Owner: answer_set. Prepared rounds: answer_set/k2/train.py; model versions: answer_set/k2/tracks.py.

giab track: decision model sees only candidates NOT on the person's GIAB list (the rest is copied);
            het/hom model is trained on inserted mutations only; rule: nothing outside the person's
            GIAB confident region.
kg track:   decision and het/hom models see every candidate of 1000 Genomes-person rounds.
"""
import argparse
import io
import json
import sys
import time
from collections import Counter, defaultdict
from contextlib import redirect_stdout

import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier

from answer_set.k2 import tracks
from answer_set.k2.features import FEATURES
from answer_set.k2.train import PREPARED, assemble, cls, load, prf, truth_of

FOLDS = 5
IN_DONOR = FEATURES.index("in_mode_donor")
IN_BED = FEATURES.index("in_donor_bed")
# The validator's v2 core weighs classes (SNP het 0.18 + hom 0.02, indel 1 bp 0.40, 2-3 bp 0.16, 4-7 bp 0.10, >= 8 bp
# 0.14); per class, one error costs ~ weight / (2 x truth records of the class). giab track: training rows are
# weighted by that cost (POINT_ALPHA: 0 = every row alike) and each class gets its own cut-off, chosen on held-out
# folds for weight x F1 - points, not counts: a missed 1-bp indel costs ~20 SNP errors. kg track: the same, and its
# cut-offs also count each extra record's cost in the false-positive part (~3.75 / N points; 1000G answers carry
# ~110 extras per round, ~0.7 points).
CLASS_W = {"snp": 0.20, "ind1": 0.40, "ind2_3": 0.16, "ind4_7": 0.10, "ind8+": 0.14}
POINT_ALPHA = 1.0
# kg track with Minos's exact rule for the person's part (kg.gatk_person_own, meta kg_own_rule): that part is copied
# and the person's other GATK calls are never written, so the model learns and is judged only on the rest
KG_OWN_RULE = False
# kg track with the round-level measurements (tracks.ROUND_FEATS, meta round_feats)
ROUND_FEATS = False
THR_GRID = [round(0.1 + 0.05 * i, 2) for i in range(17)]


def round_track(r):
    if r["info"]["mode_giab"]:
        return "giab"
    return "kg" if r["info"].get("n_dragen") else None      # unrecognised rounds belong to no track


def model():
    return HistGradientBoostingClassifier(max_iter=400, learning_rate=0.06, max_leaf_nodes=63,
                                          l2_regularization=1.0, random_state=0)


def rows(track, r, for_hom=False):
    if track == "giab":
        m = r["X"][:, IN_DONOR] == 0                  # candidates the person's list does not explain
    elif KG_OWN_RULE:
        m = np.nan_to_num(r["X"][:, FEATURES.index("gatk_person_called")]) != 1
    else:
        m = np.ones(len(r["X"]), dtype=bool)
    if for_hom:
        m = m & (r["y_any"] == 1)
    return m


def classes_of(r):
    if "_cls" not in r:
        r["_cls"] = np.array([cls(k) for k in zip(r["pos"].tolist(), r["ref"].tolist(), r["alt"].tolist())])
    return r["_cls"]


def class_cost(rs):
    """Relative cost of one error per class (SNP = 1): validator weight / truth records of the class."""
    n = Counter()
    for r in rs:
        for k in truth_of(r):
            n[cls(k)] += 1
    raw = {c: CLASS_W[c] / max(1, n[c]) for c in CLASS_W}
    return {c: raw[c] / raw["snp"] for c in CLASS_W}


def fit(track, rs, cols, cost=None):
    def w(r, m):
        if not cost or POINT_ALPHA == 0:
            return None
        return np.array([cost[c] ** POINT_ALPHA for c in classes_of(r)[m]])
    X = np.vstack([r["X"][rows(track, r)][:, cols] for r in rs])
    y = np.concatenate([r["y_any"][rows(track, r)] for r in rs])
    Xh = np.vstack([r["X"][rows(track, r, True)][:, cols] for r in rs])
    yh = np.concatenate([r["y_hom"][rows(track, r, True)] for r in rs])
    sw = [w(r, rows(track, r)) for r in rs]
    swh = [w(r, rows(track, r, True)) for r in rs]
    sw = np.concatenate(sw) if sw[0] is not None else None
    swh = np.concatenate(swh) if swh[0] is not None else None
    return model().fit(X, y, sample_weight=sw), model().fit(Xh, yh, sample_weight=swh)


def class_points(rs, oof, truths, c, t, track="giab"):
    """Mean over rounds of weight x F1 of class c with cut-off t (held-out scores; giab: the person's own records are
    copied, so they count as found; a record of ours counts only with the right het/hom). kg: minus each extra
    record's false-positive cost (in core units, 3.75 / 580 points / 70)."""
    vals = []
    for r in rs:
        pa, ph = oof[r["name"]]
        m = rows(track, r) & (classes_of(r) == c)
        if track == "giab":
            m = m & (r["X"][:, IN_BED] == 1)
        pred = m & (pa >= t)
        tp = int(np.sum(pred & (r["y_any"] == 1) & ((ph >= 0.5) == (r["y_hom"] == 1))))
        t_ = truths[r["name"]]
        own = (sum(1 for k, (z, syn) in t_.items() if not syn and cls(k) == c) if track == "giab" else
               sum(1 for k in t_ if cls(k) == c and k in r["_own"]) if KG_OWN_RULE else 0)
        n_truth = sum(1 for k in t_ if cls(k) == c)
        v = CLASS_W[c] * 2 * (own + tp) / max(1, own + int(pred.sum()) + n_truth)
        if track == "kg":
            v -= int(np.sum(pred & (r["y_any"] == 0))) * 3.75 / 580 / 70
        vals.append(v)
    return float(np.mean(vals))


def predict(pair, r, cols):
    m, h = pair
    X = r["X"][:, cols]
    return m.predict_proba(X)[:, 1], h.predict_proba(X)[:, 1]


def evaluate(track, name, rs):
    """{round: (precision, recall, f1)} and per-class counts for a track version on rounds."""
    m, h, meta, _ = tracks.load(track, name)
    cols = [tracks.col_index(f) for f in meta["features"]]
    per, agg = {}, defaultdict(Counter)
    for r in rs:
        if meta.get("round_feats") and r["X"].shape[1] == len(FEATURES):
            r = dict(r, X=tracks.with_round_feats(r["X"]))
        pa, ph = predict((m, h), r, cols)
        keys_ = list(zip(r["pos"].tolist(), r["ref"].tolist(), r["alt"].tolist()))
        pa = tracks.apply_source_thr(meta, r["X"], tracks.apply_class_thr(meta, keys_, pa), pa)
        pred = assemble(r, pa, tracks.apply_gt_rule(meta, r["X"], ph), meta["threshold"], meta.get("bed_rule", False))
        t = truth_of(r)
        per[r["name"]] = prf(pred, t)
        for k, (z, syn) in t.items():
            g = ("inserted " if syn else "own ") + cls(k)
            agg[g]["n"] += 1
            agg[g]["found"] += k in pred
            agg[g]["gt_ok"] += pred.get(k, ("", 0))[0] == z
        for k in pred:
            if k not in t:
                agg["extra " + cls(k)]["n"] += 1
    return per, agg


def train(track, rs, name):
    # kg: the read comparison's measurements when every round was prepared with it (and it worked)
    kg_rd = track == "kg" and all(r["info"].get("kg_readdiff") and "error" not in (r["info"].get("readdiff") or {})
                                  for r in rs)
    kg_hg = kg_rd and all(r["info"].get("kg_hapgeno") for r in rs)
    feats = tracks.features_for(track, kg_rd, kg_hg)
    if track == "kg" and ROUND_FEATS:
        if not kg_rd:
            sys.exit("--round-feats needs every round prepared with the read comparison")
        feats = feats + tracks.ROUND_FEATS
        for r in rs:
            r["X"] = tracks.with_round_feats(r["X"])
    cols = [tracks.col_index(f) for f in feats]
    if track == "kg" and KG_OWN_RULE:                # the person's part of each round's truth, by the rule
        from answer_set.k2 import kg
        from answer_set.k2.util import parse_region
        for r in rs:
            ch, a_, b_ = parse_region(r["info"]["region"])
            r["_own"] = set(kg.gatk_person_own(ch, a_, b_, r["info"]["kg"]["person"])[0])
    out = tracks.vdir(track, name)
    out.mkdir(parents=True, exist_ok=True)
    buf = io.StringIO()
    with redirect_stdout(buf):
        print(f"{track} track: {len(rs)} rounds, {len(cols)} measurements")
        names = sorted(r["name"] for r in rs)
        fold_of = {names[i]: n % FOLDS for n, i in enumerate(np.random.default_rng(0).permutation(len(names)))}
        oof = {}
        cost = class_cost(rs)
        if cost:
            print("error cost per class (SNP = 1):", {c: round(v, 1) for c, v in cost.items()}, "alpha", POINT_ALPHA)
        for f in range(FOLDS):
            pair = fit(track, [r for r in rs if fold_of[r["name"]] != f], cols, cost)
            joblib.dump(pair[0], out / f"fold{f}_model.joblib")
            joblib.dump(pair[1], out / f"fold{f}_hom.joblib")
            for r in rs:
                if fold_of[r["name"]] == f:
                    pa, ph = predict(pair, r, cols)
                    oof[r["name"]] = (pa, tracks.apply_gt_rule({"gt_from_gatk": track == "kg"}, r["X"], ph))
        truths = {r["name"]: truth_of(r) for r in rs}
        bed_rule = track == "giab"
        table = {t: float(np.median([prf(assemble(r, *oof[r["name"]], t, bed_rule), truths[r["name"]])[2]
                                     for r in rs])) for t in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7)}
        thr = max(table, key=table.get)
        print("held-out median F1 by threshold", {k: round(v, 5) for k, v in table.items()})
        thr_rd = None
        if track == "giab" and "only_rd" in tracks.features_for(track):
            def f1_at(t2):
                mt = {"threshold": thr, "thr_only_rd": t2}
                return float(np.median([prf(assemble(r, tracks.apply_source_thr(mt, r["X"], oof[r["name"]][0]),
                                                     oof[r["name"]][1], thr, bed_rule), truths[r["name"]])[2]
                                        for r in rs]))
            t2tab = {t2: f1_at(t2) for t2 in (thr, 0.6, 0.7, 0.8, 0.9, 0.95)}
            thr_rd = max(t2tab, key=t2tab.get)
            print("held-out median F1 by cut-off for read-comparison-only candidates",
                  {k: round(v, 5) for k, v in t2tab.items()})
        thr_class = {}
        print("held-out points by class cut-off (weight x F1, mean per round):")
        for c in CLASS_W:
            tab = {t: class_points(rs, oof, truths, c, t, track) for t in THR_GRID}
            thr_class[c] = max(tab, key=tab.get)
            print(f"   {c:7s} at {thr}: {tab[min(THR_GRID, key=lambda x: abs(x - thr))]:.6f}  best {thr_class[c]}: "
                  f"{tab[thr_class[c]]:.6f}  (+{70 * (tab[thr_class[c]] - tab[min(THR_GRID, key=lambda x: abs(x - thr))]):.4f} pts)")
        agg = defaultdict(Counter)
        print(f"\n{'region':28s} {'person':8s} {'precision':>9s} {'recall':>7s} {'F1':>7s}")
        mt = {"threshold": thr, "thr_only_rd": thr_rd, "thr_class": thr_class}
        for r in sorted(rs, key=lambda r: r["name"]):
            keys_ = list(zip(r["pos"].tolist(), r["ref"].tolist(), r["alt"].tolist()))
            pc = tracks.apply_class_thr(mt, keys_, oof[r["name"]][0])
            pred = assemble(r, tracks.apply_source_thr(mt, r["X"], pc, oof[r["name"]][0]), oof[r["name"]][1], thr, bed_rule)
            t = truths[r["name"]]
            p, rc, f1 = prf(pred, t)
            person = r["info"]["mode_donor"] if track == "giab" else r["info"]["kg"]["person"]
            print(f"{r['info']['region']:28s} {person:8s} {p:9.4f} {rc:7.4f} {f1:7.4f}")
            for k, (z, syn) in t.items():
                g = ("inserted " if syn else "own ") + cls(k)
                agg[g]["n"] += 1
                agg[g]["found"] += k in pred
            for k in pred:
                if k not in t:
                    agg["extra " + cls(k)]["n"] += 1
        print("\nper class (held-out):")
        for g in sorted(agg):
            c = agg[g]
            print(f"   {g:20s} {c['n']:7d}" if g.startswith("extra") else
                  f"   {g:20s} found {c['found']}/{c['n']} ({100 * c['found'] / c['n']:.2f}%)")
        final = fit(track, rs, cols, cost)
        joblib.dump(final[0], out / "model.joblib")
        joblib.dump(final[1], out / "hom.joblib")
        meta = {"track": track, "created": tracks.now(), "features": feats, "kg_readdiff": int(kg_rd),
                "rounds": names, "fold_of": fold_of, "threshold": thr, "heldout_f1": round(table[thr], 5),
                "bed_rule": bed_rule, "n_targets_median": 580.0, "gt_from_gatk": track == "kg",
                "thr_only_rd": thr_rd, "thr_class": thr_class, "point_alpha": POINT_ALPHA if cost else None}
        if track == "kg" and KG_OWN_RULE:
            meta["kg_own_rule"] = 1
        if track == "kg" and ROUND_FEATS:
            meta["round_feats"] = 1
        if kg_hg:
            meta["kg_hapgeno"] = min(int(r["info"]["kg_hapgeno"]) for r in rs)
        if track == "giab":
            from answer_set.k2 import features as _features
            meta["rd_1bp_factor"] = _features.RD_1BP_FACTOR
            meta["hap_genotyper"] = min(int(r["info"].get("hap_genotyper") or 0) for r in rs)
        (out / "meta.json").write_text(json.dumps(meta, indent=1))
    (out / "report.txt").write_text(buf.getvalue())
    return meta


def main():
    global POINT_ALPHA, KG_OWN_RULE, ROUND_FEATS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--track", required=True, choices=tracks.TRACKS)
    ap.add_argument("--same-rounds-as", help="train on the rounds (of this track) that a track version used")
    ap.add_argument("--promote", action="store_true")
    ap.add_argument("--exclude-rounds", help="file with round dirs to leave out (e.g. a test set)")
    ap.add_argument("--point-alpha", type=float, default=POINT_ALPHA,
                    help="weight rows by class error cost ** alpha (0: unweighted)")
    ap.add_argument("--kg-own-rule", action="store_true",
                    help="kg: the person's part is copied by Minos's exact rule; learn only the rest")
    ap.add_argument("--round-feats", action="store_true",
                    help="kg: add the round-level measurements (tracks.ROUND_FEATS: which insertion pattern the round has)")
    args = ap.parse_args()
    POINT_ALPHA = args.point_alpha
    KG_OWN_RULE = args.kg_own_rule
    ROUND_FEATS = args.round_feats
    want = None
    if args.same_rounds_as:
        n = args.same_rounds_as
        want = set(tracks.meta(args.track, n)["rounds"])
    rs = [r for r in load(want) if round_track(r) == args.track]
    if args.track == "giab" and any(f in tracks.features_for("giab") for f in ("rd_support",)):
        rs = [r for r in rs if not r["info"].get("no_rd_features") and not r["info"].get("no_rdx_features")
              and not r["info"].get("no_sib_features") and not r["info"].get("no_hg_features")]
    if args.track == "kg":                         # the kg track needs the spot measurements
        rs = [r for r in rs if not r["info"].get("no_spot_features") and not r["info"].get("no_mask_features")
              and not r["info"].get("no_annot_features") and not r["info"].get("no_gatkonly_features")]
    if args.exclude_rounds:
        ex = set(open(args.exclude_rounds).read().split())
        rs = [r for r in rs if r["name"] not in ex]
    if want is not None:
        need = {x for x in want if (PREPARED / f"{x}.json").exists()}
        got = {r["name"] for r in rs} | {r["name"] for r in load(want) if round_track(r) != args.track}
        if need - got:
            sys.exit(f"{len(need - got)} requested rounds are not prepared with the current features")
    if len(rs) < FOLDS:
        sys.exit(f"only {len(rs)} {args.track} rounds prepared; need at least {FOLDS}")
    name = tracks.new_name(args.track)
    t0 = time.time()
    meta = train(args.track, rs, name)
    print((tracks.vdir(args.track, name) / "report.txt").read_text())
    print(f"{args.track} version {name} saved ({time.time() - t0:.0f}s); held-out F1 {meta['heldout_f1']}")
    if args.promote:
        tracks.promote(args.track, name, "manual train --promote")
        print("promoted")


if __name__ == "__main__":
    main()
