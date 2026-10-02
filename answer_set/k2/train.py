"""Shared training helpers: prepared rounds, the answer set a model would output for them, scoring. Owner: answer_set.

Prepared rounds live in this clone's <cache_dir>/answer_set/k2/prepared/ (prepare.py; rebuildable): <tag>.npz (features + labels per candidate),
<tag>.json (round info) and <tag>.truth.norm.vcf.gz (the released truth, normalised like the candidates).
"""
import json

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier

from answer_set.k2.features import FEATURES
from answer_set.k2.util import CACHE, parse_region, read_alleles

PREPARED = CACHE / "prepared"
IN_DONOR = FEATURES.index("in_mode_donor")

N_BEFORE_SPOT = FEATURES.index("spot_gt_same")   # measurements appended later are padded for old rounds
N_BEFORE_MASK = FEATURES.index("in_strict_mask")
N_BEFORE_ANNOT = FEATURES.index("gatk_vqslod")
N_BEFORE_RD = FEATURES.index("rd_support")
N_BEFORE_GATKONLY = FEATURES.index("only_gatk")
N_BEFORE_RDX = FEATURES.index("only_rd")
N_BEFORE_SIB = FEATURES.index("rd_sib_max")


def load(names=None):
    """Prepared rounds (all, or those in `names`) as dicts: name, info, X, y_any, y_syn, y_tgt, y_hom, pos, ref, alt."""
    rounds = []
    for f in sorted(PREPARED.glob("*.npz")):
        if names is not None and f.stem not in names:
            continue
        info = json.loads(f.with_suffix(".json").read_text())
        z = np.load(f, allow_pickle=False)
        d = {k: z[k] for k in z.files}
        n = d["X"].shape[1]
        if n > len(FEATURES) or (n < len(FEATURES) and n not in (N_BEFORE_SPOT, N_BEFORE_MASK, N_BEFORE_ANNOT, N_BEFORE_RD,
                                                                   N_BEFORE_GATKONLY, N_BEFORE_RDX, N_BEFORE_SIB)):
            continue                               # prepared with an older feature set (re-run prepare)
        if n < len(FEATURES):                      # prepared before later measurements: unknown (NaN)
            d["X"] = np.hstack([d["X"], np.full((len(d["X"]), len(FEATURES) - n), np.nan)])
            info["no_sib_features"] = True
            info["no_rdx_features"] = n < N_BEFORE_SIB
            info["no_gatkonly_features"] = n < N_BEFORE_RDX
            info["no_rd_features"] = n < N_BEFORE_GATKONLY
            info["no_annot_features"] = n < N_BEFORE_RD
            info["no_mask_features"] = n < N_BEFORE_ANNOT
            info["no_spot_features"] = n == N_BEFORE_SPOT
        rounds.append(dict(name=f.stem, info=info, **d))
    return rounds


def model():
    return HistGradientBoostingClassifier(max_iter=400, learning_rate=0.06, max_leaf_nodes=63,
                                          l2_regularization=1.0, random_state=0)


def cls(k):
    L = abs(len(k[2]) - len(k[1]))
    return "snp" if L == 0 else "ind1" if L == 1 else "ind2_3" if L <= 3 else "ind4_7" if L <= 7 else "ind8+"


def truth_of(r):
    """{key: (zygosity, synthetic)} of the round's released truth, normalised like the candidates."""
    t = {}
    for k, z, info in read_alleles(PREPARED / f"{r['name']}.truth.norm.vcf.gz", keep_info=True):
        t[k] = (z, "SYNTHETIC" in info)
    return t


def assemble(r, p_any, p_hom, thr, bed_rule=True):
    """The answer set this pipeline would output for round r (key -> (zyg, synthetic))."""
    from answer_set.k2.assemble import assemble_records
    from answer_set.k2.util import Bed, donor_bed, load_donor_alleles
    keys = list(zip(r["pos"].tolist(), r["ref"].tolist(), r["alt"].tolist()))
    info = r["info"]
    donor = {}
    if info["mode_giab"]:
        if "_donor" not in r:
            chrom, a, b = parse_region(info["region"])
            r["_donor"] = load_donor_alleles(chrom, a, b)[info["mode_donor"]]
        donor = r["_donor"]
    allowed = None
    if info["mode_giab"] and bed_rule:
        if "_bed" not in r:
            r["_bed"] = Bed(donor_bed(info["mode_donor"]), parse_region(info["region"])[0])
        bed = r["_bed"]
        allowed = lambda k: bed.contains(k[0]) and bed.contains(k[0] + len(k[1]) - 1)  # noqa: E731
    recs = assemble_records(keys, p_any, p_hom, thr, info["mode_giab"], donor, allowed)
    return {(x[0], x[1], x[2]): (x[3], x[4]) for x in recs}


def prf(pred, truth):
    tp = sum(1 for k in pred if k in truth)
    p, rc = tp / max(1, len(pred)), tp / max(1, len(truth))
    return p, rc, 2 * p * rc / max(1e-9, p + rc)
