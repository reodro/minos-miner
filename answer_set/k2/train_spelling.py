"""Train the spelling model of replaced windows (hapwin.respell) on a GIAB version that has window decisions, as a
new version (the version's files + spelling.joblib and leave-out fold spelling models).

  python -m answer_set.k2.train_spelling --version g013_20261003T1402
Examples, per training round of the version: windows whose called pair is the truth's DNA - every one the window
decision can fix (the window base's leave-out answer wrong) and a 5% sample of the rest; for each called sequence
every spelling option (hapwin.spelling_options), labelled 1 when all its events are truth records. Same folds as the
version. Reported on the fixable windows, held out: exactly the truth's records with the fewest-events spelling vs
the learned one.
Owner: answer_set.
"""
import argparse
import json
import shutil
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

import joblib
import numpy as np
import pysam
from sklearn.ensemble import HistGradientBoostingClassifier

from answer_set.k2 import hapwin, tracks
from answer_set.k2.features import FEATURES
from answer_set.k2.train import assemble, load, truth_of
from answer_set.k2.util import CACHE, Bed, donor_bed, parse_region, ref_path


def chooser():
    return HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, random_state=0)


def round_rows(args):
    """[(window id, sequence id, features, label, fixable, truth records, events)] for one training round."""
    tag, base = args
    r = load({tag})[0]
    m, h, meta, _ = tracks.load("giab", base, exclude_round=tag)
    cols = [FEATURES.index(f) for f in meta["features"]]
    keys = list(zip(r["pos"].tolist(), r["ref"].tolist(), r["alt"].tolist()))
    p_raw = m.predict_proba(r["X"][:, cols])[:, 1]
    ph = h.predict_proba(r["X"][:, cols])[:, 1]
    pa = tracks.apply_source_thr(meta, r["X"], tracks.apply_class_thr(meta, keys, p_raw), p_raw)
    answer = assemble(r, pa, ph, meta["threshold"], meta.get("bed_rule", True))
    chrom = parse_region(r["info"]["region"])[0]
    p_of = dict(zip(keys, p_raw))
    truth = truth_of(r)
    tk = sorted(truth)
    tp = [k[0] for k in tk]
    fa = pysam.FastaFile(str(ref_path(chrom)))
    wpath = CACHE / "prepare" / tag / "hg_windows.pkl.gz"
    if not wpath.exists():
        return tag, []
    wins = hapwin.load_windows(wpath)
    cands = hapwin.candidates(wins, answer, p_of, r["_donor"], Bed(donor_bed(r["info"]["mode_donor"]), chrom))
    y, ok = hapwin.label(cands, truth)
    fixable = {(c[0], c[1]) for c, yy, o in zip(cands, y, ok) if yy and not o}
    rng = np.random.default_rng(sum(map(ord, tag)))
    rows = []
    for wi, (s, e, w) in enumerate(wins):
        if "pair" not in w:
            continue
        fx = (s, e) in fixable
        if not fx and rng.random() > 0.05:
            continue
        T = hapwin._inside(tk, tp, s, e)
        if not T:
            continue
        ref = w["ref"]
        tps = hapwin.pairs_of(ref, s, [(k, truth[k][0]) for k in T])
        if not tps or tuple(sorted(w["pair"])) not in tps:
            continue
        for hi, hh in enumerate(sorted(set(w["pair"]) - {ref})):
            opts = hapwin.spelling_options(fa, chrom, s, ref, hh)
            if not opts or (len(opts) < 2 and not fx):
                continue
            X = hapwin.spelling_features(opts, p_of)
            for (ev, _), x in zip(opts, X):
                rows.append((wi, hi, x, int(set(ev) <= set(T)), fx, frozenset(T), ev))
    return tag, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", required=True, help="GIAB version with window decisions")
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    t0 = time.time()
    vmeta = tracks.meta("giab", args.version)
    base = vmeta.get("window_base")
    if vmeta.get("window_thr") is None or not base:
        sys.exit(f"{args.version} has no window decisions")
    with ProcessPoolExecutor(args.workers) as ex:
        data = {t: rows for t, rows in ex.map(round_rows, [(t, base) for t in vmeta["rounds"]]) if rows}
    fold_of = vmeta["fold_of"]
    res, fold_models = Counter(), {}
    for f in sorted({fold_of[t] for t in data}):
        tr = [x for t in data if fold_of[t] != f for x in data[t]]
        mdl = chooser().fit(np.array([x[2] for x in tr]), np.array([x[3] for x in tr]))
        fold_models[f] = mdl
        for t in data:
            if fold_of[t] != f:
                continue
            q = mdl.predict_proba(np.array([x[2] for x in data[t]]))[:, 1]
            by = {}
            for x, qq in zip(data[t], q):
                by.setdefault(x[0], {}).setdefault(x[1], []).append((qq, x))
            for haps in by.values():
                first = next(iter(haps.values()))[0][1]
                if not first[4]:
                    continue
                T = set(first[5])
                learned = set().union(*[max(v, key=lambda t_: t_[0])[1][6] for v in haps.values()])
                fewest = set().union(*[min(v, key=lambda t_: (t_[1][2][0], t_[1][2][7]))[1][6] for v in haps.values()])
                res["fixable"] += 1
                res["learned"] += learned == T
                res["fewest"] += fewest == T
    n = max(1, res["fixable"])
    print(f"{len(data)} rounds; fixable windows {res['fixable']}, held out, exactly the truth's records: fewest events "
          f"{res['fewest']} ({100 * res['fewest'] / n:.1f}%), learned {res['learned']} ({100 * res['learned'] / n:.1f}%)")
    final = chooser().fit(np.array([x[2] for t in data for x in data[t]]), np.array([x[3] for t in data for x in data[t]]))
    name = tracks.new_name("giab")
    src_rel, src_dev, dst = tracks.vdir("giab", args.version), tracks.dev_dir("giab", args.version), tracks.dev_dir("giab", name)
    dst.mkdir(parents=True, exist_ok=True)
    for d in (src_dev, src_rel):                          # the version's own files (dev copy has the fold models)
        if d.exists():
            for f in d.iterdir():
                if f.suffix == ".joblib" and not (dst / f.name).exists():
                    shutil.copy2(f, dst / f.name)
    joblib.dump(final, dst / "spelling.joblib")
    for f, mdl in fold_models.items():
        joblib.dump(mdl, dst / f"fold{f}_spelling.joblib")
    meta = dict(vmeta, created=tracks.now(), spelling=True, spelling_of=args.version,
                spelling_features=hapwin.SPELL_FEATS,
                spelling_heldout={"fixable": res["fixable"], "fewest": res["fewest"], "learned": res["learned"]})
    (dst / "meta.json").write_text(json.dumps(meta, indent=1))
    tracks.log("giab", "train_spelling", version=name, of=args.version)
    print(f"giab version {name} = {args.version} + spelling model saved ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    sys.exit(main())
