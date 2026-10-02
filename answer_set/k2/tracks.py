"""Two independent model tracks, one per round format.

  giab : rounds whose person is one of the 7 GIAB people (background copied from GIAB's list;
         the model only decides which other candidates are inserted mutations)
  kg   : rounds whose person is from 1000 Genomes (the model decides every variant, using the
         person's DRAGEN calls, the 1000 Genomes panel and the GATK callset as evidence)

Each track has its own models (decision + het/hom), its own measurement list, its own threshold,
its own versions, active version, history and rollback. Released models are part of the code (in git, a few MB);
everything else a developer produces stays in this clone's cache:

  answer_set/k2/models/<track>/ACTIVE                       the released version main runs
  answer_set/k2/models/<track>/<name>/                      released: model.joblib, hom.joblib, meta.json
  <cache_dir>/answer_set/k2/models/<track>/versions/<name>/ every trained version: + fold*_{model,hom}.joblib
                                                            (leave-out models for backtests), report.txt
  <cache_dir>/answer_set/k2/models/<track>/history.jsonl    promotions / rollbacks and why

Releasing a version = promote (copies its three files into answer_set/k2/models and sets ACTIVE) + commit.

  python -m answer_set.k2.tracks                       # list both tracks
  python -m answer_set.k2.tracks --track giab --use g007_...
  python -m answer_set.k2.tracks --track kg --rollback
"""
import argparse
import json
import shutil
import time
from pathlib import Path

from answer_set.k2.features import FEATURES
from answer_set.k2.util import CACHE

TRACKS = ("giab", "kg")
RELEASED = Path(__file__).resolve().parent / "models"      # in git
DEV = CACHE / "models"                                       # this clone
RELEASE_FILES = ("model.joblib", "hom.joblib", "meta.json")
PREFIX = {"giab": "g", "kg": "k"}

# measurements that only make sense for one format
GIAB_ONLY = ["in_mode_donor", "in_donor_bed", "donor_dist_log", "donor_indel_near",
             "rd_support", "rd_cover", "rd_frac",
             "hap_alt", "hap_ref", "hap_n", "hap_frac", "hap_rank", "hap_margin", "hap_cluster",
             "only_rd", "rd_orig", "rd_orig_frac", "rd_sib_max", "rd_sib_ratio"]
KG_ONLY = ["kg_mode", "in_dragen", "dragen_qual", "dragen_gq", "dragen_hom", "in_panel", "only_dragen",
           "dragen_and_gatk", "gatk_person_checked", "gatk_person_called", "gatk_person_same_gt",
           "spot_gt_same", "disagree10", "disagree30", "disagree100", "gatk_alleles_near",
           "in_strict_mask", "in_pilot_mask", "in_giab_easy", "in_giab_repeat",
           "gatk_vqslod", "gatk_an_frac", "gatk_hwe_log", "gatk_exchet_log", "gatk_mq", "gatk_mqranksum",
           "gatk_fs", "gatk_sor", "gatk_site_dp_log", "gatk_inbreeding", "panel_site",
           "gatk_rec_nalt", "gatk_allele_idx", "gatk_rec_overlap", "only_gatk", "gatk_person_hom"]
DROP = ["mode_giab"]                       # constant within a track


def apply_gt_rule(meta, X, p_hom):
    """kg track: where GATK called the variant for this person, Minos's key uses GATK's het/hom
    (100% of 576,785 background variants in 43 rounds, incl. all 1,385 GATK/DRAGEN disagreements)."""
    import numpy as np
    from answer_set.k2.features import FEATURES
    if not meta.get("gt_from_gatk") or "gatk_person_hom" not in FEATURES or X.shape[1] != len(FEATURES):
        return p_hom
    v = X[:, FEATURES.index("gatk_person_hom")]
    return np.where(np.isnan(v), p_hom, v)


def apply_source_thr(meta, X, p_any):
    """giab track: candidates found ONLY by the read comparison must pass their own (stricter)
    cut-off meta['thr_only_rd']; expressed by moving their score across the main cut-off."""
    import numpy as np
    from answer_set.k2.features import FEATURES
    t2 = meta.get("thr_only_rd")
    if t2 is None or "only_rd" not in FEATURES or X.shape[1] != len(FEATURES):
        return p_any
    thr = meta["threshold"]
    o = X[:, FEATURES.index("only_rd")] == 1
    out = p_any.copy()
    out[o & (p_any >= t2)] = np.maximum(p_any[o & (p_any >= t2)], thr)
    out[o & (p_any < t2)] = np.minimum(p_any[o & (p_any < t2)], thr - 1e-6)
    return out


def features_for(track):
    other = KG_ONLY if track == "giab" else GIAB_ONLY
    return [f for f in FEATURES if f not in other and f not in DROP]


def track_of_info(info):
    return "giab" if info["mode_giab"] else "kg"


def now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def dev_dir(track, name):
    return DEV / track / "versions" / name


def released_dir(track, name):
    return RELEASED / track / name


def vdir(track, name):
    """Where a version's files are: the released copy if there is one, else this clone's training output."""
    r = released_dir(track, name)
    return r if (r / "meta.json").exists() else dev_dir(track, name)


def all_versions(track):
    names = {p.parent.name for p in (DEV / track / "versions").glob("*/meta.json")}
    names |= {p.parent.name for p in (RELEASED / track).glob("*/meta.json")}
    return sorted(names)


def new_name(track):
    return f"{PREFIX[track]}{len(all_versions(track)) + 1:03d}_{time.strftime('%Y%m%dT%H%M', time.gmtime())}"


def active_name(track):
    p = RELEASED / track / "ACTIVE"
    return p.read_text().strip() if p.exists() else None


def meta(track, name):
    return json.loads((vdir(track, name) / "meta.json").read_text())


def log(track, event, **kw):
    (DEV / track).mkdir(parents=True, exist_ok=True)
    with open(DEV / track / "history.jsonl", "a") as fh:
        fh.write(json.dumps({"time": now(), "event": event, **kw}) + "\n")


def release(track, name):
    """Copy a trained version's model.joblib, hom.joblib and meta.json into answer_set/k2/models (to be committed)."""
    src, dst = dev_dir(track, name), released_dir(track, name)
    if (dst / "meta.json").exists():
        return dst
    dst.mkdir(parents=True, exist_ok=True)
    for f in RELEASE_FILES:
        shutil.copy2(src / f, dst / f)
    return dst


def _set_active(track, name):
    release(track, name)
    (RELEASED / track / "ACTIVE").write_text(name + "\n")


def promote(track, name, reason):
    prev = active_name(track)
    _set_active(track, name)
    log(track, "promote", version=name, previous=prev, reason=reason)


def previous_active(track):
    h = DEV / track / "history.jsonl"
    if not h.exists():
        return None
    cur = active_name(track)
    for line in reversed(h.read_text().splitlines()):
        e = json.loads(line)
        if e["event"] in ("promote", "rollback") and e.get("version") == cur:
            return e.get("previous")
    return None


def rollback(track, reason):
    prev = previous_active(track)
    if not prev:
        return None
    cur = active_name(track)
    _set_active(track, prev)
    log(track, "rollback", version=prev, previous=cur, reason=reason)
    return prev


def load(track, name=None, exclude_round=None):
    """(decision model, het/hom model, meta, description) for a track version. With exclude_round (backtests of a
    released round): the leave-out fold model that never saw it - those live in the cache of the clone that trained
    the version; without them the full model is used, and the description says so."""
    import joblib
    name = name or active_name(track)
    if not name:
        raise FileNotFoundError(f"no active version for track {track}")
    d = vdir(track, name)
    m = meta(track, name)
    if exclude_round and exclude_round in m.get("fold_of", {}):
        f = m["fold_of"][exclude_round]
        fd = dev_dir(track, name)
        if (fd / f"fold{f}_model.joblib").exists():
            return (joblib.load(fd / f"fold{f}_model.joblib"), joblib.load(fd / f"fold{f}_hom.joblib"), m,
                    f"{track} track {name}: fold{f} (never saw this round)")
        return (joblib.load(d / "model.joblib"), joblib.load(d / "hom.joblib"), m,
                f"{track} track {name}: all its rounds (fold model not on this machine - it saw this round)")
    return (joblib.load(d / "model.joblib"), joblib.load(d / "hom.joblib"), m, f"{track} track {name}: all its rounds")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--track", choices=TRACKS)
    ap.add_argument("--use")
    ap.add_argument("--rollback", action="store_true")
    args = ap.parse_args()
    if args.use or args.rollback:
        if not args.track:
            raise SystemExit("--track is required with --use / --rollback")
        if args.use:
            if args.use not in all_versions(args.track):
                raise SystemExit(f"unknown {args.track} version {args.use}")
            promote(args.track, args.use, "manual")
        if args.rollback:
            print("rolled back to", rollback(args.track, "manual") or "(no previous version)")
    for t in ([args.track] if args.track else TRACKS):
        act = active_name(t)
        print(f"[{t} track] active: {act}")
        for v in all_versions(t):
            m = meta(t, v)
            print(f"  {'*' if v == act else ' '} {v}  trained on {len(m['rounds'])} {t} rounds  "
                  f"held-out F1 {m.get('heldout_f1')}  threshold {m.get('threshold')}  "
                  f"{len(m['features'])} measurements")


if __name__ == "__main__":
    main()
