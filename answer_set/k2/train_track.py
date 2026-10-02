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
    else:
        m = np.ones(len(r["X"]), dtype=bool)
    if for_hom:
        m = m & (r["y_any"] == 1)
    return m


def fit(track, rs, cols):
    X = np.vstack([r["X"][rows(track, r)][:, cols] for r in rs])
    y = np.concatenate([r["y_any"][rows(track, r)] for r in rs])
    Xh = np.vstack([r["X"][rows(track, r, True)][:, cols] for r in rs])
    yh = np.concatenate([r["y_hom"][rows(track, r, True)] for r in rs])
    return model().fit(X, y), model().fit(Xh, yh)


def predict(pair, r, cols):
    m, h = pair
    X = r["X"][:, cols]
    return m.predict_proba(X)[:, 1], h.predict_proba(X)[:, 1]


def evaluate(track, name, rs):
    """{round: (precision, recall, f1)} and per-class counts for a track version on rounds."""
    m, h, meta, _ = tracks.load(track, name)
    cols = [FEATURES.index(f) for f in meta["features"]]
    per, agg = {}, defaultdict(Counter)
    for r in rs:
        pa, ph = predict((m, h), r, cols)
        pa = tracks.apply_source_thr(meta, r["X"], pa)
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
    cols = [FEATURES.index(f) for f in tracks.features_for(track)]
    out = tracks.vdir(track, name)
    out.mkdir(parents=True, exist_ok=True)
    buf = io.StringIO()
    with redirect_stdout(buf):
        print(f"{track} track: {len(rs)} rounds, {len(cols)} measurements")
        names = sorted(r["name"] for r in rs)
        fold_of = {names[i]: n % FOLDS for n, i in enumerate(np.random.default_rng(0).permutation(len(names)))}
        oof = {}
        for f in range(FOLDS):
            pair = fit(track, [r for r in rs if fold_of[r["name"]] != f], cols)
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
        agg = defaultdict(Counter)
        print(f"\n{'region':28s} {'person':8s} {'precision':>9s} {'recall':>7s} {'F1':>7s}")
        mt = {"threshold": thr, "thr_only_rd": thr_rd}
        for r in sorted(rs, key=lambda r: r["name"]):
            pred = assemble(r, tracks.apply_source_thr(mt, r["X"], oof[r["name"]][0]), oof[r["name"]][1], thr, bed_rule)
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
        final = fit(track, rs, cols)
        joblib.dump(final[0], out / "model.joblib")
        joblib.dump(final[1], out / "hom.joblib")
        meta = {"track": track, "created": tracks.now(), "features": tracks.features_for(track),
                "rounds": names, "fold_of": fold_of, "threshold": thr, "heldout_f1": round(table[thr], 5),
                "bed_rule": bed_rule, "n_targets_median": 580.0, "gt_from_gatk": track == "kg",
                "thr_only_rd": thr_rd}
        if track == "giab":
            from answer_set.k2 import features as _features
            meta["rd_1bp_factor"] = _features.RD_1BP_FACTOR
        (out / "meta.json").write_text(json.dumps(meta, indent=1))
    (out / "report.txt").write_text(buf.getvalue())
    return meta


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--track", required=True, choices=tracks.TRACKS)
    ap.add_argument("--same-rounds-as", help="train on the rounds (of this track) that a track version used")
    ap.add_argument("--promote", action="store_true")
    ap.add_argument("--exclude-rounds", help="file with round dirs to leave out (e.g. a test set)")
    args = ap.parse_args()
    want = None
    if args.same_rounds_as:
        n = args.same_rounds_as
        want = set(tracks.meta(args.track, n)["rounds"])
    rs = [r for r in load(want) if round_track(r) == args.track]
    if args.track == "giab" and any(f in tracks.features_for("giab") for f in ("rd_support",)):
        rs = [r for r in rs if not r["info"].get("no_rd_features") and not r["info"].get("no_rdx_features")
              and not r["info"].get("no_sib_features")]
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
