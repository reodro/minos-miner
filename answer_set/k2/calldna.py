"""DNA-aware call decision for 1000G indel candidates. Owner: answer_set.

The kg call model sees alignment and read-comparison counts. For an indel candidate that is not one of the person's
GATK calls, a second call model (call_dna.joblib) replaces p_any: it sees the version's own measurements plus the
repeat shape at the event and homrep's reads judged by their DNA (spanning reads carrying the event, showing the
reference, or other; CIGAR stutter).
Trained per fold on k033's 78 new-pattern rounds (test/answer_set/kg_call_feats.py, kg_call_eval.py): with k033's
class cut-offs, missed/false indels 916/183 -> 904/162. On only for a version whose meta sets call_dna; build.py
applies it right after the call model, before the class cut-offs.
"""
import numpy as np

from answer_set.k2.homrep import read_feats, repeat_at


def feats(bam, fa, chrom, k):
    """repeat shape (unit, span, copies, event in units, whole units) + homrep.read_feats for candidate k"""
    u, span, cp = repeat_at(fa, chrom, k)
    L = abs(len(k[1]) - len(k[2]))
    return [u, span, cp, L / u if u else 0.0, float(L % u == 0) if u else 0.0] + read_feats(bam, fa, chrom, k, span)


def apply(p_any, keys, X, person_called, bam, fa, chrom, model):
    """p_any with the DNA call model's value for every indel candidate that is not a person's GATK call
    -> (p_any, n replaced)"""
    out = p_any.copy()
    idx = [i for i, k in enumerate(keys) if len(k[1]) != len(k[2]) and np.nan_to_num(person_called[i]) != 1]
    if not idx:
        return out, 0
    N = np.array([feats(bam, fa, chrom, keys[i]) for i in idx], dtype=float)
    out[idx] = model.predict_proba(np.hstack([X[idx], N]))[:, 1]
    return out, len(idx)
