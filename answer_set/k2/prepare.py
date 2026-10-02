"""Prepare released rounds for training: candidates + features from the BAM, labels from the released truth.
Owner: answer_set.

Output, per round, in <cache_dir>/answer_set/k2/prepared/ (see train.py): <tag>.npz, <tag>.json, <tag>.truth.norm.vcf.gz.
Rounds come from SETTINGS.rounds_dir (python -m backtest.data --last N); work files go to the cache.

  python -m answer_set.k2.prepare                          # every stored round not prepared yet
  python -m answer_set.k2.prepare --tag 20261001T164100Z --redo
"""
import argparse
import json
import os
import sys
import time
from collections import Counter

import numpy as np

from answer_set.k2.features import FEATURES, extract
from answer_set.k2.train import PREPARED, cls
from answer_set.k2.util import CACHE, normalize_vcf, read_alleles, ref_path
from common.settings import SETTINGS
from common.types import Round


def label_round(rnd: Round, work):
    """Released truth and targets of the round, normalised like the candidates."""
    PREPARED.mkdir(parents=True, exist_ok=True)
    tn = normalize_vcf(os.path.join(rnd.released_dir, "truth.vcf.gz"), PREPARED / f"{rnd.tag}.truth.norm.vcf.gz",
                       ref_path(rnd.chrom), rnd.region)
    mn = normalize_vcf(os.path.join(rnd.released_dir, "mutations.vcf.gz"), work / "mut.norm.vcf.gz",
                       ref_path(rnd.chrom), rnd.region)
    truth = {k: (z, "SYNTHETIC" in info) for k, z, info in read_alleles(tn, keep_info=True)}
    targets = {k for k, _, _ in read_alleles(mn)}
    return truth, targets


def prepare(rnd: Round, threads=8):
    t0 = time.time()
    work = CACHE / "prepare" / rnd.tag
    work.mkdir(parents=True, exist_ok=True)
    truth, targets = label_round(rnd, work)
    keys, X, info, _, _ = extract(rnd.bam, rnd.region, work, threads, reuse=True)
    kset = set(keys)
    y_any = np.array([k in truth for k in keys], dtype=np.int8)
    y_syn = np.array([bool(truth.get(k, (0, False))[1]) for k in keys], dtype=np.int8)
    y_tgt = np.array([k in targets for k in keys], dtype=np.int8)
    y_hom = np.array([truth.get(k, ("", 0))[0] == "hom" for k in keys], dtype=np.int8)
    np.savez_compressed(PREPARED / f"{rnd.tag}.npz", X=X, y_any=y_any, y_syn=y_syn, y_tgt=y_tgt,
                        y_hom=y_hom, pos=np.array([k[0] for k in keys]),
                        ref=np.array([k[1] for k in keys]), alt=np.array([k[2] for k in keys]))
    # recall ceiling: which truth alleles never appear as candidates
    miss, tot = Counter(), Counter()
    for k, (z, syn) in truth.items():
        c = ("syn " if syn else "bg ") + cls(k)
        tot[c] += 1
        miss[c] += k not in kset
    info.update(round_dir=rnd.tag, region=rnd.region, n_features=len(FEATURES), n_truth=len(truth),
                n_targets=len(targets), n_syn=sum(1 for v in truth.values() if v[1]),
                ceiling={c: round(1 - miss[c] / tot[c], 4) for c in sorted(tot)},
                truth_counts=dict(tot), seconds=round(time.time() - t0, 1))
    (PREPARED / f"{rnd.tag}.json").write_text(json.dumps(info, indent=1))
    return info


def main():
    from backtest.data import load
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", action="append", help="round tag (repeatable); default: every stored round")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--redo", action="store_true")
    args = ap.parse_args()
    for tag in sorted(args.tag or os.listdir(SETTINGS.rounds_dir)):
        if not args.redo and (PREPARED / f"{tag}.npz").exists():
            continue
        try:
            rnd = load(tag)
            if not os.path.exists(os.path.join(rnd.released_dir, "truth.vcf.gz")):
                print(f"{tag}: truth not released, skipped", flush=True)
                continue
            info = prepare(rnd, args.threads)
        except Exception as e:  # noqa: BLE001  (one bad round must not stop the rest)
            print(f"{tag}: FAILED {e}", flush=True)
            continue
        print(f"{tag} {info['region']:27s} mode={'GIAB:' + info['mode_donor'] if info['mode_giab'] else 'other'} "
              f"cands={info['n_candidates']} truth={info['n_truth']} syn={info['n_syn']} {info['seconds']}s", flush=True)


if __name__ == "__main__":
    sys.exit(main())
