"""Train the window-decision model (hapwin.py) on top of a GIAB track version, as a new version.

  python -m answer_set.k2.train_window --base g010_20261003T0918
The new version = the base version's files + window.joblib (and leave-out fold window models) + meta window_thr.
Per training round of the base: the base's leave-out answer (its fold model never saw the round), the genotyper's
windows saved by prepare (work/hg_windows.pkl.gz), the windows inside the confident region whose called pair differs
from that answer, labelled by the released truth (is the called pair the truth's DNA?). Same folds as the base. The
cut-off maximises points on held-out folds: a window fixed gains, a window broken loses, each by the class cost of
its events (train_track.class_cost: a 1-bp indel ~20 SNP errors).
Owner: answer_set.
"""
import argparse
import json
import shutil
import sys
import time

import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier

from answer_set.k2 import hapwin, tracks
from answer_set.k2.features import FEATURES
from answer_set.k2.train import assemble, cls, load, truth_of
from answer_set.k2.train_track import class_cost
from answer_set.k2.util import CACHE, Bed, donor_bed, parse_region

GRID = [round(0.3 + 0.05 * i, 2) for i in range(14)]


def window_model():
    return HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, l2_regularization=1.0,
                                          random_state=0)


def round_rows(r, base, cost):
    """(features, label, answer right, gain if fixed, loss if broken) for a prepared round's candidate windows."""
    tag = r["name"]
    wpath = CACHE / "prepare" / tag / "hg_windows.pkl.gz"
    if not wpath.exists():
        return None
    m, h, meta, _ = tracks.load("giab", base, exclude_round=tag)
    cols = [FEATURES.index(f) for f in meta["features"]]
    keys = list(zip(r["pos"].tolist(), r["ref"].tolist(), r["alt"].tolist()))
    p_raw = m.predict_proba(r["X"][:, cols])[:, 1]
    ph = h.predict_proba(r["X"][:, cols])[:, 1]
    pa = tracks.apply_source_thr(meta, r["X"], tracks.apply_class_thr(meta, keys, p_raw), p_raw)
    answer = assemble(r, pa, ph, meta["threshold"], meta.get("bed_rule", True))
    bed = Bed(donor_bed(r["info"]["mode_donor"]), parse_region(r["info"]["region"])[0])
    cands = hapwin.candidates(hapwin.load_windows(wpath), answer, dict(zip(keys, p_raw)), r["_donor"], bed)
    if not cands:
        return None
    y, ok = hapwin.label(cands, truth_of(r))
    gain = np.array([max([cost[cls(k)] for k in c[4]] + [cost[cls(k)] for k, _ in c[3]] + [1.0]) for c in cands])
    return np.vstack([c[5] for c in cands]), y, np.array(ok), gain


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True, help="GIAB track version whose answers the windows correct")
    args = ap.parse_args()
    t0 = time.time()
    bmeta = tracks.meta("giab", args.base)
    rs = load(set(bmeta["rounds"]))
    cost = class_cost(rs)
    data = {}
    for r in rs:
        x = round_rows(r, args.base, cost)
        if x is not None:
            data[r["name"]] = x
    print(f"{len(data)}/{len(rs)} rounds with genotyper windows; candidate windows "
          f"{sum(len(v[1]) for v in data.values())}, called pair right {sum(int(v[1].sum()) for v in data.values())}, "
          f"answer right {sum(int(v[2].sum()) for v in data.values())}")
    fold_of = bmeta["fold_of"]
    folds = sorted({fold_of[n] for n in data})
    oof, fold_models = {}, {}
    for f in folds:
        tr = [n for n in data if fold_of[n] != f]
        mdl = window_model().fit(np.vstack([data[n][0] for n in tr]), np.concatenate([data[n][1] for n in tr]))
        fold_models[f] = mdl
        for n in data:
            if fold_of[n] == f:
                oof[n] = mdl.predict_proba(data[n][0])[:, 1]

    def points(t):
        fix = brk = 0.0
        nf = nb = 0
        for n, (X, y, ok, g) in data.items():
            take = oof[n] >= t
            nf += int(np.sum(take & (y == 1) & ~ok)); nb += int(np.sum(take & (y == 0) & ok))
            fix += float(np.sum(g[take & (y == 1) & ~ok])); brk += float(np.sum(g[take & (y == 0) & ok]))
        return fix - brk, nf, nb
    table = {t: points(t) for t in GRID}
    thr = max(table, key=lambda t: table[t][0])
    for t in GRID:
        v, nf, nb = table[t]
        print(f"  cut-off {t:.2f}: windows fixed {nf:5d}, broken {nb:5d}, net (SNP-error units) {v:9.1f}{'  <-' if t == thr else ''}")
    final = window_model().fit(np.vstack([d[0] for d in data.values()]), np.concatenate([d[1] for d in data.values()]))
    name = tracks.new_name("giab")
    src, dst = tracks.vdir("giab", args.base), tracks.dev_dir("giab", name)
    dst.mkdir(parents=True, exist_ok=True)
    for f in list(tracks.RELEASE_FILES) + [f"fold{i}_{k}.joblib" for i in range(10) for k in ("model", "hom")]:
        for d in (src, tracks.dev_dir("giab", args.base)):
            if (d / f).exists():
                shutil.copy2(d / f, dst / f)
                break
    joblib.dump(final, dst / "window.joblib")
    for f, mdl in fold_models.items():
        joblib.dump(mdl, dst / f"fold{f}_window.joblib")
    meta = dict(bmeta, created=tracks.now(), window_thr=thr, window_features=hapwin.FEATS, window_base=args.base,
                window_heldout={"fixed": table[thr][1], "broken": table[thr][2], "net_snp_units": round(table[thr][0], 1)})
    (dst / "meta.json").write_text(json.dumps(meta, indent=1))
    tracks.log("giab", "train_window", version=name, base=args.base, threshold=thr)
    print(f"giab version {name} = {args.base} + window decisions (cut-off {thr}) saved ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    sys.exit(main())
