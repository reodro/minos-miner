"""Prepare released rounds for training: candidates + features from the BAM, labels from the released truth.
Owner: answer_set.

Output, per round, in <cache_dir>/answer_set/k2/prepared/ (see train.py): <tag>.npz, <tag>.json, <tag>.truth.norm.vcf.gz.
Rounds come from SETTINGS.rounds_dir (python -m backtest.data --last N); work files go to the cache.

  python -m answer_set.k2.prepare                          # every stored round not prepared yet
  python -m answer_set.k2.prepare --tag 20261001T164100Z --redo
"""
import argparse
import itertools
import json
import os
import sys
import time
from collections import Counter

import numpy as np
import pysam

from answer_set.k2 import features
from answer_set.k2.features import FEATURES, extract
from answer_set.k2.train import PREPARED, cls
from answer_set.k2.util import CACHE, normalize_vcf, read_alleles, ref_path
from common.settings import SETTINGS
from common.types import Round

# Teach the events of genotyper windows that give the truth's DNA as right (window_labels; models from g038 on).
# Off by default so tables prepared for earlier versions stay as they were; set with --window-labels and
# MINOS_PREPARED to keep the tables apart.
WINDOW_LABELS = False


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


def _apply(seq, w0, alleles):
    out, i = [], 0
    for pos, ref, alt in sorted(alleles):
        j = pos - w0
        if j < i:
            return None
        out.append(seq[i:j])
        out.append(alt)
        i = j + len(ref)
    return "".join(out) + seq[i:]


def equivalent_labels(keys, truth, fasta, chrom, near=30):
    """{candidate index: zygosity} for candidates that are not in the truth as spelled but give the SAME DNA as one or
    two inserted truth records within `near` bp - alone or together with one neighbouring candidate (Minos often inserts a
    deletion and an insertion up to ~20 bp apart that together give the DNA of one SNP or one shifted indel; hap.py, which
    scores the answer set, compares DNA, not spellings). Only where the truth's own spelling is not a candidate, so two
    spellings of one event are never both taught as right."""
    fa = pysam.FastaFile(str(fasta))
    syn = sorted(k for k, v in truth.items() if v[1])
    kset = set(keys)
    by_pos = sorted((k[0], i) for i, k in enumerate(keys) if k not in truth)
    pos_only = [p for p, _ in by_pos]
    out = {}
    import bisect
    for i, c in enumerate(keys):
        if c in truth:
            continue
        j0, j1 = bisect.bisect_left(syn, (c[0] - near,)), bisect.bisect_right(syn, (c[0] + len(c[1]) + near,))
        tn = [t for t in syn[j0:j1] if t not in kset]
        if not tn:
            continue
        # truth records nearby that ARE candidates (donor or inserted): the same context on both sides
        ctx = [t for t in truth if t in kset and abs(t[0] - c[0]) <= near + 10]
        k0, k1 = bisect.bisect_left(pos_only, c[0] - near), bisect.bisect_right(pos_only, c[0] + near)
        partners = [keys[x] for _, x in by_pos[k0:k1] if keys[x] != c]
        span = [c] + tn + partners + ctx
        lo = min(x[0] for x in span) - 5
        hi = max(x[0] + len(x[1]) for x in span) + 5
        seq = fa.fetch(chrom, lo - 1, hi).upper()
        targets = {}
        for n in (1, 2):
            for S in itertools.combinations(tn, n):
                h = _apply(seq, lo, list(S) + ctx)
                if h is not None:
                    targets[h] = "hom" if all(truth[t][0] == "hom" for t in S) else "het"
        h = _apply(seq, lo, [c] + ctx)
        if h in targets:
            out[i] = targets[h]
            continue
        for p in partners:
            h = _apply(seq, lo, [c, p] + ctx)
            if h in targets:
                out[i] = targets[h]
                break
    return out


def window_labels(keys, truth, windows_path):
    """{candidate index: zygosity} for candidates that are events of a genotyper window (hapgeno, hg_windows.pkl.gz)
    whose called pair of local DNA sequences is the truth's pair (some phasing of the truth's het records), while the
    candidate is not in the truth as spelled. hap.py compares DNA, so writing such a window's events scores like the
    truth's own records; equivalent_labels covers one candidate (or two) against one or two truth records, this covers
    any number. 28 development rounds (2026-10-10): at 59 of 112 clustered inserted 1-bp indels g026 missed, the
    genotyper's called pair was the truth's, and the model rejected most of its events (taught as wrong). A window is
    skipped when a truth record inside it is a candidate the window did not call (two spellings never both taught as
    right) or a record crosses its edges."""
    import gzip
    import pickle
    try:
        wins = pickle.load(gzip.open(windows_path))
    except (OSError, EOFError, pickle.UnpicklingError):
        return {}
    index = {k: i for i, k in enumerate(keys)}
    tkeys = sorted(truth)
    tpos = [k[0] for k in tkeys]
    import bisect
    out = {}
    for w in wins:
        if not (isinstance(w, (tuple, list)) and len(w) == 3 and isinstance(w[2], dict)):
            continue
        s, e, v = w
        pair, called, ref = v.get("pair"), v.get("called") or [], v.get("ref")
        if not pair or not called or not ref:
            continue
        j0, j1 = bisect.bisect_left(tpos, s - 60), bisect.bisect_right(tpos, e)
        near = [k for k in tkeys[j0:j1] if k[0] <= e and k[0] + len(k[1]) - 1 >= s]
        if any(k[0] < s or k[0] + len(k[1]) - 1 > e for k in near):
            continue                                        # a truth record crosses the window's edge
        evs = {ev: c.get("hom") for c in called for ev in c["events"]}
        if any(k in index and k not in evs for k in near):
            continue                                        # the truth's own spelling is a candidate here
        if any(ev[0] < s or ev[0] + len(ev[1]) - 1 > e for ev in evs):
            continue
        hom = [k for k in near if truth[k][0] == "hom"]
        het = [k for k in near if truth[k][0] != "hom"]
        if len(het) > 8:
            continue
        pairs = set()
        for mask in range(2 ** max(0, len(het) - 1) if het else 1):
            h1 = _apply(ref, s, hom + [k for i, k in enumerate(het) if mask >> i & 1])
            h2 = _apply(ref, s, hom + [k for i, k in enumerate(het) if not mask >> i & 1])
            if h1 is not None and h2 is not None:
                pairs.add(tuple(sorted((h1, h2))))
        if tuple(sorted(pair)) not in pairs:
            continue
        for ev, h in evs.items():
            if ev in index and ev not in truth:
                out[index[ev]] = "hom" if h else "het"
    return out


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
    eq = equivalent_labels(keys, truth, ref_path(rnd.chrom), rnd.chrom)
    for i, z in eq.items():                     # same DNA as inserted truth records, other spelling
        y_any[i], y_syn[i], y_hom[i] = 1, 1, int(z == "hom")
    weq = window_labels(keys, truth, work / "hg_windows.pkl.gz") if WINDOW_LABELS else {}
    for i, z in weq.items():                    # events of a genotyper window that gives the truth's DNA
        if not y_any[i]:
            y_any[i], y_syn[i], y_hom[i] = 1, 1, int(z == "hom")
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
                hap_genotyper=int(features.HAP_GENOTYPER), equivalent_labels=len(eq),
                window_labels=len(weq) if WINDOW_LABELS else None,
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
    ap.add_argument("--no-hap", action="store_true", help="without the haplotype genotyper (models trained before g008)")
    ap.add_argument("--no-kg-rd", action="store_true", help="1000G rounds without the read comparison (models before k008)")
    ap.add_argument("--kg-hg", action="store_true", help="1000G rounds with the haplotype genotyper too")
    ap.add_argument("--kg-site", type=int, default=0, help="1000G rounds with the per-site candidate rule version "
                                                           "(features.KG_SITE_RULE: 1 = k026, 2 = k027); set MINOS_PREPARED "
                                                           "to keep these apart")
    ap.add_argument("--hc", action="store_true", help="GIAB rounds with the GATK evidence (hcrun.py); set "
                                                      "MINOS_PREPARED to keep these tables apart")
    ap.add_argument("--kg-hc", action="store_true", help="1000G rounds with the GATK evidence (hcrun.py); set "
                                                         "MINOS_PREPARED to keep these tables apart")
    ap.add_argument("--hapgeno-version", type=int, default=0, help="genotyper version (hapgeno.VERSION if 0; 2 = "
                                                                   "repeat-aware windows); set MINOS_PREPARED")
    ap.add_argument("--window-labels", action="store_true", help="GIAB: also teach events of genotyper windows that "
                                                                 "give the truth's DNA as right (window_labels)")
    args = ap.parse_args()
    global WINDOW_LABELS
    WINDOW_LABELS = args.window_labels
    from answer_set.k2 import hapgeno
    features.HAP_GENOTYPER = 0 if args.no_hap else (args.hapgeno_version or hapgeno.VERSION)
    features.KG_READDIFF = 0 if args.no_kg_rd else 1
    features.KG_HAPGENO = hapgeno.VERSION if args.kg_hg else 0
    features.KG_SITE_RULE = args.kg_site
    from answer_set.k2 import hcrun
    features.HC_EVIDENCE = hcrun.VERSION if args.hc else 0
    features.KG_HC = hcrun.VERSION if args.kg_hc else 0
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
