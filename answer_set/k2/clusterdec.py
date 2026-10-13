"""Edit-based cluster decoder (GIAB). Owner: answer_set.

Minos puts each inserted mutation on one copy at random, and every read it edits for a copy carries all of that copy's
edits (training rounds 2026-10-13: pairs of inserted records within 20 bp on the same copy - median 7 reads carry both,
0 carry only one). So a cluster's two local DNA copies are read off the reads exactly; what is hard is writing them as
the events Minos wrote. Per cluster (>= 2 candidates within 20 bp, >= 1 indel):
  1. the genotyper's best DNA pair (hapgeno.genotype_window, edit-aware with the donor's original reads);
  2. every minimal-cost spelling of each copy (affine costs, all optimal alignment paths, left-aligned events);
  3. a learned chooser picks the spelling Minos used (held-out training folds: 99.7% of copies vs 97.0% for "fewest
     events");
  4. only copies carried by EDITED reads (differ from / not found among the originals) are written; the donor's
     records and the donor copy's DNA are never rewritten;
  5. where the decode and the per-variant answer disagree on the inserted records, a learned arbiter decides
     (held-out training folds, cut-off 0.6: 265 clusters fixed, 38 broken in 108 rounds).
On only for a model whose meta names it (cluster_decoder); build.py runs it after assemble_records.
"""
import collections
import sys

import joblib
import numpy as np
import pysam

from answer_set.k2 import hapcand, hapgeno
from answer_set.k2.util import trim

MM, GO, GE = 4, 3, 1
GAP = 20                  # candidates this close form a cluster
MAX_SPAN = 70
EDIT_MIN, ORIG_MAX_SHARE = 3, 0.2


def spellings(ref, hap, s, fa, chrom, cap=64):
    """All minimal-cost decompositions of hap vs ref (Gotoh DP, 3 states), as sets of left-aligned events."""
    n, m = len(ref), len(hap); INF = 10 ** 9
    M = [[INF] * (m + 1) for _ in range(n + 1)]; D = [[INF] * (m + 1) for _ in range(n + 1)]; I = [[INF] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0
    for i in range(n + 1):
        for j in range(m + 1):
            if i and j: M[i][j] = min(M[i-1][j-1], D[i-1][j-1], I[i-1][j-1]) + (0 if ref[i-1] == hap[j-1] else MM)
            if i: D[i][j] = min(M[i-1][j] + GO, D[i-1][j] + GE, I[i-1][j] + GO)
            if j: I[i][j] = min(M[i][j-1] + GO, I[i][j-1] + GE, D[i][j-1] + GO)
    best = min(M[n][m], D[n][m], I[n][m])
    out = set()
    # enumerate optimal paths backward; ops list of (op, i) in ref coords
    sys.setrecursionlimit(10000)
    def back(st, i, j, ops):
        if len(out) >= cap: return
        if i == 0 and j == 0:
            out.add(tuple(ops[::-1])); return
        if st == "M":
            v = M[i][j]; c = 0 if ref[i-1] == hap[j-1] else MM
            for p, T in (("M", M), ("D", D), ("I", I)):
                if T[i-1][j-1] + c == v: back(p, i - 1, j - 1, ops + ([("X", i - 1, hap[j-1])] if c else []))
        elif st == "D":
            v = D[i][j]
            for p, T, g in (("M", M, GO), ("D", D, GE), ("I", I, GO)):
                if T[i-1][j] + g == v: back(p, i - 1, j, ops + [("D", i - 1, None)])
        else:
            v = I[i][j]
            for p, T, g in (("M", M, GO), ("I", I, GE), ("D", D, GO)):
                if T[i][j-1] + g == v: back(p, i, j - 1, ops + [("I", i, hap[j-1])])
    for st, T in (("M", M), ("D", D), ("I", I)):
        if T[n][m] == best: back(st, n, m, [])
    res = set()
    for ops in out:
        ev = []; k = 0
        while k < len(ops):
            op, i, b = ops[k]
            if op == "X": ev.append((s + i, ref[i], b)); k += 1; continue
            run = [ops[k]]
            while k + 1 < len(ops) and ops[k+1][0] == op and (ops[k+1][1] == run[-1][1] + (1 if op == "D" else 0)): k += 1; run.append(ops[k])
            k += 1
            if op == "D":
                i0 = run[0][1]
                if i0 == 0: break
                ev.append((s + i0 - 1, ref[i0 - 1:i0 + len(run)], ref[i0 - 1]))
            else:
                i0 = run[0][1]
                if i0 == 0: break
                ev.append((s + i0 - 1, ref[i0 - 1], ref[i0 - 1] + "".join(x[2] for x in run)))
        else:
            res.add(tuple(sorted({hapcand.left_align(fa, chrom, *e) for e in ev})))
    return res, best



def _hp(fa, chrom, p):
    seq = fa.fetch(chrom, p - 12, p + 12).upper(); i = 11
    a = b = i
    while a > 0 and seq[a - 1] == seq[i]: a -= 1
    while b + 1 < len(seq) and seq[b + 1] == seq[i]: b += 1
    return b - a + 1


def spell_feats(ev, fa, chrom, rd, bc):
    n = len(ev); snp = [e for e in ev if len(e[1]) == len(e[2])]; ind = [e for e in ev if len(e[1]) != len(e[2])]
    L = [abs(len(e[1]) - len(e[2])) for e in ind]
    pos = sorted(e[0] for e in ev); gaps = [b - a for a, b in zip(pos, pos[1:])]
    hps = [_hp(fa, chrom, e[0] + 1) for e in ind]
    ext = sum(1 for e in ind if len(e[2]) > len(e[1]) and len(set(e[2][1:])) == 1 and e[2][1] == fa.fetch(chrom, e[0], e[0] + 1).upper())
    sup = [rd.get(e, (0, 0, 0))[0] for e in ev]; ori = [rd.get(e, (0, 0, 0))[2] for e in ev]
    return [n, len(snp), len(ind), sum(1 for x in L if x == 1), sum(L), max(L or [0]), min(gaps or [99]), max(gaps or [0]),
            max(hps or [0]), sum(hps) / len(hps) if hps else 0, ext, sum(1 for e in ev if e in bc), sum(1 for e in ev if e in rd),
            min(sup or [0]), sum(sup), min(ori or [0]), sum(1 for e in ind if len(e[2]) > len(e[1])), sum(1 for e in ind if len(e[1]) > len(e[2]))]


def is_edited(c):
    m = c["edit"] + c["unk"]
    return m >= EDIT_MIN and c["orig"] <= ORIG_MAX_SHARE * (m + c["orig"])


def decode(r, ch):
    """{event: zygosity} the edited copies' chosen spellings give, minus donor events; chooser confidence; reads"""
    donor_ev = set(r["donor_in"])
    for c, sp, p in ch:
        if not is_edited(c): donor_ev |= set(sp)
    out, conf, rmin = {}, 1.0, 99
    for c, sp, p in ch:
        if not is_edited(c): continue
        conf = min(conf, p); rmin = min(rmin, c["reads"])
        for ev in sp:
            if ev in donor_ev: continue
            out[ev] = "hom" if r["hom"] or out.get(ev) == "het" else "het"
    return out, conf, rmin


def arb_feats(r, ch, d, conf, rmin):
    qw, P = r["ours"], r.get("p", {})
    only_d = [k for k in d if k not in qw]; only_q = [k for k in qw if k not in d]
    zyg = sum(1 for k in d if k in qw and d[k] != qw[k])
    pd = [P.get(k, (0.0, 0.0))[0] for k in only_d]; pq = [P.get(k, (0.0, 0.0))[0] for k in only_q]
    cand_d = sum(1 for k in only_d if k in P)
    ind = lambda ks: sum(1 for k in ks if len(k[1]) != len(k[2]))
    eds = [(c["edit"], c["orig"], c["unk"], c["n_ori"], c["reads"]) for c, sp, p in ch] + [(0, 0, 0, 0, 0)] * 2
    nsp = [len(c["spellings"]) for c, sp, p in ch] + [0, 0]
    return [len(only_d), len(only_q), zyg, max(pd or [0]), min(pd or [1]), max(pq or [0]), min(pq or [1]), cand_d,
            ind(only_d), ind(only_q), conf, rmin, r["post"], r["anchored"], r["ref_reads"], int(r["hom"]), len(r["donor_in"]),
            eds[0][0], eds[0][1], eds[0][2], eds[0][3], eds[1][0], eds[1][1], eds[1][2], eds[1][3], nsp[0], nsp[1],
            sum(1 for k in qw if k in d)]


def clusters(cands):
    ks = sorted(cands); used = set(); out = []
    for i, k in enumerate(ks):
        if k in used: continue
        cl = [k]
        for o in ks[i + 1:]:
            if o[0] - (cl[-1][0] + len(cl[-1][1])) <= GAP: cl.append(o)
            else: break
        if len(cl) < 2 or not any(len(x[1]) != len(x[2]) for x in cl): continue
        used.update(cl)
        s, e = cl[0][0] - 6, max(x[0] + len(x[1]) for x in cl) + 6
        if e - s <= MAX_SPAN: out.append((s, e))
    return out


def cluster_row(bam, fa, chrom, ob, s, e, donor, rd, bc):
    """the decoder's view of one cluster (as in training): copies with spellings + features, or None"""
    reads, haps, post, L = hapgeno.genotype_window(bam, fa, chrom, s, e, ob, 2)
    if not post: return None
    p0, a, b, _ = post[0]; ref = haps[0]
    assigned = collections.Counter(); st_of = collections.defaultdict(collections.Counter); ori_of = collections.defaultdict(set)
    for row, rr in zip(L, reads):
        mx = max(row); hi = row.index(mx); assigned[hi] += 1; st_of[hi][rr[2]] += 1
        if rr[2] == "edit": ori_of[hi].add(rr[3])
    copies = []
    for hh in sorted({a, b}):
        if hh == 0: copies.append(None); continue
        sp, _ = spellings(ref, haps[hh], s, fa, chrom)
        sp = sorted(sp)
        copies.append({"reads": assigned[hh], "edit": st_of[hh]["edit"], "orig": st_of[hh]["orig"], "unk": st_of[hh]["?"],
                       "n_ori": len(ori_of[hh]), "spellings": sp, "feats": [spell_feats(ev, fa, chrom, rd, bc) for ev in sp]})
    return {"s": s, "e": e, "post": p0, "anchored": len(reads), "hom": a == b, "ref_reads": assigned[0], "copies": copies,
            "donor_in": [k for k in donor if s <= k[0] and k[0] + len(k[1]) - 1 <= e]}


def read_inputs(workdir):
    """(read-scan candidates, read-comparison edits {key: (sup, cov, ori)}) from a build's work dir"""
    bc = set()
    try:
        for r in pysam.VariantFile(str(workdir / "cand.norm.vcf.gz")):
            for a in r.alts or ():
                if a not in ("<*>", "*"): bc.add(trim(r.pos, r.ref, a))
    except Exception:  # noqa: BLE001
        pass
    rd = {}
    p = workdir / "readdiff.vcf"
    if p.exists():
        for l in open(p):
            if l[0] == "#": continue
            f = l.split("	"); info = dict(x.split("=") for x in f[7].strip().split(";"))
            rd[trim(int(f[1]), f[3], f[4])] = (int(info["SUP"]), int(info["COV"]), int(info.get("ORI", info["SUP"])))
    return bc, rd


def apply(recs, bam_path, fasta, chrom, orig_bam, donor, workdir, keys, p_any, p_hom, models, cut):
    """recs (assemble_records output) with disputed clusters replaced where the arbiter says so -> (recs, n)"""
    chooser, arbiter = models
    bc, rd = read_inputs(workdir)
    cand = set(bc) | {k for k, v in rd.items() if v[0] >= 2}
    bam, fa = pysam.AlignmentFile(str(bam_path)), pysam.FastaFile(str(fasta))
    ob = pysam.AlignmentFile(str(orig_bam)) if orig_bam else None
    P = {k: (float(x), float(y)) for k, x, y in zip(keys, p_any, p_hom)}
    syn = {(r[0], r[1], r[2]): r for r in recs if r[4]}
    done, n = [], 0
    for s, e in clusters(cand):
        r = cluster_row(bam, fa, chrom, ob, s, e, donor, rd, bc)
        if r is None: continue
        if len([c for c in r["copies"] if c and c["spellings"]]) != len([c for c in r["copies"] if c]): continue
        ch = []
        for c in r["copies"]:
            if not c: continue
            if len(c["spellings"]) == 1:
                ch.append((c, c["spellings"][0], 1.0)); continue
            pr = chooser.predict_proba(np.array([f + [len(c["spellings"])] for f in c["feats"]]))[:, 1]
            i = int(np.argmax(pr)); ch.append((c, c["spellings"][i], float(pr[i])))
        r["ours"] = {k: v[3] for k, v in syn.items() if s <= k[0] and k[0] + len(k[1]) - 1 <= e}
        r["p"] = {k: v for k, v in P.items() if s - 5 <= k[0] <= e + 5}
        d, conf, rmin = decode(r, ch)
        if d == r["ours"]: continue
        if arbiter.predict_proba(np.array([arb_feats(r, ch, d, conf, rmin)]))[0, 1] < cut: continue
        done.append((s, e, d)); n += 1
    if not done: return recs, 0
    out = [x for x in recs if not (x[4] and any(s <= x[0] and x[0] + len(x[1]) - 1 <= e for s, e, _ in done))]
    have = {(x[0], x[1], x[2]) for x in out}
    for s, e, d in done:
        out.extend((k[0], k[1], k[2], z, True, 1.0) for k, z in d.items() if k not in have)
    return out, n
