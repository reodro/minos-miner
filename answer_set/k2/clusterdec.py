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
After that rewrite, shorten_repeated takes a homozygous run back out of a nearby heterozygous
allele, and then drop_second_spelling removes a second overlapping indel the forest scored under
the calling threshold. Shortening runs first: the longer allele is the true one-base allele
plus the run, and dropping it first would lose that allele. spell_run_edge then writes an
uncalled 2 bp insertion inside a long homopolymer as the run base at the left edge plus the
other base where the reads put it.

1000G (meta kg_one_spelling): the same biology, without the GIAB chooser. A cluster is rewritten only when the
inserted calls already are a spelling of the edited copies' DNA and another minimal spelling is the one to keep
(most read-comparison support, then fewest events). Copies the model never called, including a 1-of-20 allele,
are left alone. Donor records stay.
"""
import collections
import sys
from pathlib import Path

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


def _hp_at(seq, i):
    if i < 0 or i >= len(seq): return 0
    a = b = i
    while a > 0 and seq[a - 1] == seq[i]: a -= 1
    while b + 1 < len(seq) and seq[b + 1] == seq[i]: b += 1
    return b - a + 1


def shared_indel_rows(r, ch, ref, fa, chrom):
    """Shared-indel pattern (2026-10-14): Minos puts a 1-bp indel k on both copies inside a run and a second edit next to
    it on one copy; that copy's DNA then reads as if k were absent, so every per-variant and DNA method calls k het.
    For each 1-bp indel the answer calls het where one decoded copy carries k and the other equals (ref + k) plus one
    small event: [(k, features)] for the resolver (training rounds, held out: AUC 0.946; cut 0.7 fixes 52, breaks 14
    in 108 rounds)."""
    out = []
    if len(ch) != 2: return out
    s = r["s"]; P = r.get("p", {}); qw = r["ours"]
    for k, z in qw.items():
        if abs(len(k[1]) - len(k[2])) != 1 or z != "het": continue
        carr = [hapgeno._carries(ref, s, c["hap"], hapgeno._ed(ref, c["hap"]), k) for c, sp, p in ch]
        if sum(carr) != 1: continue
        (A, spA, pA) = [x for x, cr in zip(ch, carr) if cr][0]
        (B, spB, pB) = [x for x, cr in zip(ch, carr) if not cr][0]
        rk = hapgeno._apply_key(ref, s, k)
        if rk is None: continue
        cb, ck = hapgeno._ed(ref, B["hap"]), hapgeno._ed(rk, B["hap"])
        if ck > cb + 1: continue
        sp, _ = spellings(rk, B["hap"], s, fa, chrom)
        single = [ev for ev in sp if len(ev) == 1 and abs(len(ev[0][1]) - len(ev[0][2])) <= 5]
        if not single: continue
        e = single[0][0]; i = k[0] - s
        bev = spB[0] if len(spB) == 1 else None
        out.append((k, [ck - cb, cb, len(spB), _hp_at(ref, i + 1), abs(e[0] - k[0]), int(len(e[1]) == len(e[2])),
                        len(e[2]) - len(e[1]), int(len(k[2]) > len(k[1])), A["reads"], B["reads"], A["edit"], B["edit"],
                        A["orig"], B["orig"], A["unk"], B["unk"], A["n_ori"], B["n_ori"], r["ref_reads"],
                        P.get(k, (0, 0))[0], P.get(k, (0, 0))[1], P.get(bev, (0, 0))[0] if bev else -1,
                        int(bev in qw) if bev else -1, r["post"], len(qw)]))
    return out


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
        copies.append({"hap": haps[hh], "reads": assigned[hh], "edit": st_of[hh]["edit"], "orig": st_of[hh]["orig"], "unk": st_of[hh]["?"],
                       "n_ori": len(ori_of[hh]), "spellings": sp, "feats": [spell_feats(ev, fa, chrom, rd, bc) for ev in sp]})
    return {"s": s, "e": e, "post": p0, "anchored": len(reads), "hom": a == b, "ref_reads": assigned[0], "copies": copies,
            "donor_in": [k for k in donor if s <= k[0] and k[0] + len(k[1]) - 1 <= e]}


def pick_spelling(spellings, rd):
    """The one spelling of a copy: most read-comparison support, then the fewest events."""
    def key(ev):
        return (sum(rd.get(e, (0, 0, 0))[0] for e in ev), -len(ev))
    return max(spellings, key=key)


def _apply_events(ref, s, events):
    h = ref
    for k in sorted(events, key=lambda k: -k[0]):
        h2 = hapgeno._apply_key(h, s, k)
        if h2 is None:
            return None
        h = h2
    return h


def _carries_copy(ref, s, copy, k):
    return hapgeno._carries(ref, s, copy["hap"], hapgeno._ed(ref, copy["hap"]), k)


def _reconstructs(ref, s, hap, spellings, carried):
    """True when `carried` is a spelling of hap: applying it gives hap, or it is one listed spelling, or the union
    of every listed spelling (both ways of writing the same DNA were emitted)."""
    if _apply_events(ref, s, carried) == hap:
        return True
    have = set(carried)
    sets = [set(sp) for sp in spellings if sp]
    if have in sets:
        return True
    return bool(sets) and len(sets) >= 2 and have == set().union(*sets)


def spelling_replacement(r, rd, ref):
    """{event: zyg} to write instead of the cluster's inserted calls, or None.

    None when the inserted calls are already that one spelling, when no copy is edited, or when the inserted calls
    are not themselves a spelling of the edited DNA (a different variant, not a different way to write the same one).
    """
    present = [c for c in r["copies"] if c]
    if not present or any(not c["spellings"] for c in present):
        return None
    ch = [(c, pick_spelling(c["spellings"], rd), 1.0) for c in present]
    d, _, _ = decode(r, ch)
    if d == r.get("ours", {}):
        return None
    edited = [c for c, sp, p in ch if is_edited(c)]
    if not edited:
        return None
    ours = list(r.get("ours", {}))
    s = r["s"]
    for k in ours:
        if not any(_carries_copy(ref, s, c, k) for c in edited):
            return None
    for c in edited:
        carried = [k for k in ours if _carries_copy(ref, s, c, k)]
        if not _reconstructs(ref, s, c["hap"], c["spellings"], carried):
            return None
    return d


def apply_one_spelling(recs, bam_path, fasta, chrom, orig_bam, donor, workdir):
    """1000G: one spelling per edited haplotype. Rewrites a cluster only when its inserted calls are a different
    minimal spelling of the same edited DNA (spelling_replacement). Donor records stay. -> (recs, n replaced)"""
    bc, rd = read_inputs(workdir)
    cand = set(bc) | {k for k, v in rd.items() if v[0] >= 2}
    ins = [(r[0], r[1], r[2]) for r in recs if r[4] and len(r[1]) != len(r[2])]
    if not ins:
        return recs, 0
    bam, fa = pysam.AlignmentFile(str(bam_path)), pysam.FastaFile(str(fasta))
    ob = pysam.AlignmentFile(str(orig_bam)) if orig_bam else None
    try:
        syn = {(r[0], r[1], r[2]): r for r in recs if r[4]}
        done = []
        for s, e in clusters(cand):
            if not any(s <= k[0] and k[0] + len(k[1]) - 1 <= e for k in ins):
                continue
            r = cluster_row(bam, fa, chrom, ob, s, e, donor, rd, bc)
            if r is None:
                continue
            r["ours"] = {k: v[3] for k, v in syn.items() if s <= k[0] and k[0] + len(k[1]) - 1 <= e}
            ref = fa.fetch(chrom, s - 1, e).upper()
            d = spelling_replacement(r, rd, ref)
            if d:
                done.append((s, e, d))
    finally:
        bam.close()
        fa.close()
        if ob is not None:
            ob.close()
    if not done:
        return recs, 0
    out = [x for x in recs if not (x[4] and any(s <= x[0] and x[0] + len(x[1]) - 1 <= e for s, e, _ in done))]
    have = {(x[0], x[1], x[2]) for x in out}
    for s, e, d in done:
        out.extend((k[0], k[1], k[2], z, True, 1.0) for k, z in d.items() if k not in have)
    return out, len(done)


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


def _net_len(r):
    return len(r[2]) - len(r[1])


def _similar_indel(a, b):
    """Same-direction indels whose spans touch and whose length change differs by one base.
    A wider gap merges a real nearby indel into the spelling (chr16:19946638 on the 15:30 round)."""
    na, nb = _net_len(a), _net_len(b)
    if na == 0 or nb == 0 or (na > 0) != (nb > 0):
        return False
    if abs(abs(na) - abs(nb)) != 1:
        return False
    a1, b1 = a[0] + len(a[1]) - 1, b[0] + len(b[1]) - 1
    return not (a1 < b[0] - 1 or b1 < a[0] - 1)


def drop_second_spelling(recs, p_of, thr):
    """Drop a synthetic indel the forest scored under thr when a similar indel scores higher.

    Two heterozygous copies of that edit become one homozygous record. A donor record is never
    removed and its genotype stays as published. A lone below-threshold indel stays: the cluster
    decoder is allowed to add one spelling the per-variant cutoff missed.
    """
    def score(r):
        if not r[4]:
            return 1.0
        return float(p_of.get((r[0], r[1], r[2]), 0.0))

    drop, promote = set(), set()
    for i, a in enumerate(recs):
        if not a[4] or _net_len(a) == 0 or score(a) >= thr:
            continue
        best, best_p = None, score(a)
        for j, b in enumerate(recs):
            if i == j or j in drop or not _similar_indel(a, b):
                continue
            pb = score(b)
            if pb > best_p:
                best, best_p = j, pb
        if best is None:
            continue
        drop.add(i)
        b = recs[best]
        if b[4] and a[3] == "het" and b[3] == "het":
            promote.add(best)
    out = []
    for i, r in enumerate(recs):
        if i in drop:
            continue
        if i in promote:
            r = r[:3] + ("hom",) + tuple(r[4:])
        out.append(r)
    return out


def _payload(ref, alt):
    """('ins' or 'del', the bases that differ) after the shared ends are trimmed, or (None, None)."""
    a, b = ref, alt
    while a and b and a[0] == b[0]:
        a, b = a[1:], b[1:]
    while a and b and a[-1] == b[-1]:
        a, b = a[:-1], b[:-1]
    if a == "" and b:
        return "ins", b
    if b == "" and a:
        return "del", a
    return None, None


def _window_seq(fa, chrom, events, s, e):
    seq = fa.fetch(chrom, s - 1, e).upper()
    for pos, ref, alt in sorted(events, reverse=True):
        i = pos - s
        if i < 0 or seq[i:i + len(ref)] != ref:
            return None
        seq = seq[:i] + alt + seq[i + len(ref):]
    return seq


def shorten_repeated(recs, fa, chrom):
    """A homozygous insertion of one repeated base is already on both copies. A synthetic
    heterozygous insertion within 40 bp that starts or ends with those bases writes them again.
    Shorten the heterozygous allele so the two records together equal that allele on its own.
    A donor record is never rewritten. An ambiguous pair is left as called."""
    homs = []
    for r in recs:
        if r[3] != "hom" or len(r[1]) == len(r[2]):
            continue
        kind, pay = _payload(r[1], r[2])
        if pay and len(set(pay)) == 1:
            homs.append((r, kind, pay))
    if not homs:
        return recs
    replace = {}
    for i, e in enumerate(recs):
        if not e[4] or e[3] == "hom" or len(e[1]) == len(e[2]):
            continue
        kind_e, pay_e = _payload(e[1], e[2])
        if not pay_e:
            continue
        found = None
        for h, kind_h, pay_h in homs:
            if h[0] == e[0] or abs(h[0] - e[0]) > 40 or kind_h != kind_e or pay_e == pay_h:
                continue
            if pay_e.startswith(pay_h):
                rest = pay_e[len(pay_h):]
            elif pay_e.endswith(pay_h):
                rest = pay_e[:-len(pay_h)]
            else:
                continue
            if not rest or pay_e not in (e[2] if kind_e == "ins" else e[1]):
                continue
            if kind_e == "ins":
                new = (e[0], e[1], e[2].replace(pay_e, rest, 1))
            else:
                new_ref = e[1].replace(pay_h, "", 1)
                if not new_ref:
                    continue
                new = (e[0], new_ref, e[2])
            if len(new[1]) == len(new[2]):
                continue
            s = min(h[0], e[0]) - 3
            end = max(h[0] + len(h[1]), e[0] + len(e[1])) + len(pay_e) + 5
            alone = _window_seq(fa, chrom, [(e[0], e[1], e[2])], s, end)
            both = _window_seq(fa, chrom, [(h[0], h[1], h[2]), new], s, end)
            if alone is None or alone != both:
                continue
            if found is not None:
                found = None
                break
            found = new
        if found is not None and found not in {(r[0], r[1], r[2]) for r in recs}:
            replace[i] = found
    if not replace:
        return recs
    out = []
    for i, r in enumerate(recs):
        if i in replace:
            r = replace[i] + tuple(r[3:])
        out.append(r)
    return out


RUN_EDGE_MIN = 12          # homopolymer this long, or the 2 bp insertion is its own allele
RUN_EDGE_READS = 6
RUN_EDGE_TOP = 5
RUN_EDGE_ORI = 4           # original reads edited with the foreign base; fewer is alignment noise


def _run_of(fa, chrom, pos, base):
    """(start, end) 1-based inclusive of the homopolymer of `base` that covers pos, or None."""
    seq = fa.fetch(chrom, max(0, pos - 80), pos + 80).upper()
    origin = max(0, pos - 80)
    idx = pos - 1 - origin
    if idx < 0 or idx >= len(seq) or seq[idx] != base:
        return None
    a = b = idx
    while a > 0 and seq[a - 1] == base:
        a -= 1
    while b + 1 < len(seq) and seq[b + 1] == base:
        b += 1
    if b - a + 1 < RUN_EDGE_MIN:
        return None
    return origin + a + 1, origin + b + 1


def _split_slip(fa, chrom, event, run_base, ta):
    """One inserted run base plus a different base, written as the left-edge run insertion and the
    other base alone. None unless those two records rebuild the same sequence."""
    pos, ref, alt = event
    kind, pay = _payload(ref, alt)
    if kind != "ins" or not pay or len(pay) != 2 or pay[0] != run_base or pay[1] == run_base:
        return None
    foreign = pay[1]
    prev = fa.fetch(chrom, ta - 2, ta - 1).upper()
    if not prev:
        return None
    edge = trim(ta - 1, prev, prev + run_base)
    rem = trim(pos, ref, alt.replace(pay, foreign, 1))
    s = min(edge[0], rem[0], pos) - 3
    e = max(edge[0] + len(edge[1]), rem[0] + len(rem[1]), pos + len(ref)) + 8
    alone = _window_seq(fa, chrom, [(pos, ref, alt)], s, e)
    both = _window_seq(fa, chrom, [edge, rem], s, e)
    if not alone or alone != both:
        return None
    return edge, rem


def spell_run_edge(recs, bam_path, fa, chrom, rd, cand_vcf=None):
    """Write an uncalled 2 bp insertion inside a long homopolymer as two 1 bp records.

    The reads show one run base plus a foreign base as a single insertion. The score wants the
    run base at the left edge of the run and the foreign base where it sits. That split is used
    only when several of the donor's original reads were edited with the foreign base; the same
    pileup pattern without that edit is a gap that slid along the run. A donor record is never
    added or removed, and a position that already carries a different allele is left alone.
    Returns (records, number of alleles added).
    """
    if not cand_vcf or not Path(str(cand_vcf)).exists():
        return recs, 0
    seeds = set()
    for rec in pysam.VariantFile(str(cand_vcf)).fetch():
        for alt in rec.alts or ():
            if alt in ("*", "<*>"):
                continue
            kind, pay = _payload(rec.ref, alt)
            if kind == "ins" and pay and len(pay) == 2 and pay[0] != pay[1]:
                span = _run_of(fa, chrom, rec.pos, pay[0])
                if span and span[0] + 2 <= rec.pos <= span[1] - 2:
                    seeds.add((span[0], span[1], pay[0]))
    if not seeds:
        return recs, 0
    bam = pysam.AlignmentFile(str(bam_path))
    have = {(r[0], r[1], r[2]) for r in recs}
    at_pos = {}
    for r in recs:
        at_pos.setdefault(r[0], set()).add((r[1], r[2]))
    add = []
    try:
        for ta, tb, base in seeds:
            lo = max(1, ta - 4)
            hi = min(tb + 4, fa.get_reference_length(chrom))
            ref = fa.fetch(chrom, lo - 1, hi).upper()
            seqs = collections.Counter()
            for read in bam.fetch(chrom, lo - 1, hi):
                if (read.mapping_quality < 20 or read.is_secondary or read.is_supplementary
                        or not read.query_sequence or not read.cigartuples):
                    continue
                if read.reference_start > lo - 1 or (read.reference_end or 0) < hi:
                    continue
                q, qi, rp = [], 0, read.reference_start
                for op, length in read.cigartuples:
                    if op in (0, 7, 8):
                        for _ in range(length):
                            if lo - 1 <= rp < hi:
                                q.append(read.query_sequence[qi])
                            rp += 1
                            qi += 1
                    elif op == 1:
                        if lo - 1 <= rp <= hi:
                            q.append(read.query_sequence[qi:qi + length])
                        qi += length
                    elif op == 2:
                        rp += length
                    elif op == 4:
                        qi += length
                    elif op == 3:
                        rp += length
                    elif op == 5:
                        qi += length
                seqs["".join(q)] += 1
            total = sum(seqs.values())
            if total < RUN_EDGE_READS or seqs.get(ref, 0) / total > 0.20:
                continue
            hap, n = seqs.most_common(1)[0]
            if hap == ref or n < RUN_EDGE_TOP or n / total < 0.50:
                continue
            options, _ = spellings(ref, hap, lo, fa, chrom)
            if len(options) != 1 or len(next(iter(options))) != 1:
                continue
            event = next(iter(options))[0]
            pieces = _split_slip(fa, chrom, event, base, ta)
            if pieces is None or not (ta + 2 <= event[0] <= tb - 2):
                continue
            if event in have:
                continue
            foreign = _payload(event[1], event[2])[1][1]
            ori = max((v[2] for k, v in rd.items() if k[0] == event[0] and foreign in k[2]), default=0)
            if ori < RUN_EDGE_ORI:
                continue
            fresh = []
            blocked = False
            for key in pieces:
                if key in have:
                    continue
                others = at_pos.get(key[0], set())
                if any(a != (key[1], key[2]) for a in others):
                    blocked = True
                    break
                fresh.append(key)
            if blocked:
                continue
            for key in fresh:
                add.append(key + ("hom", True, 0.99))
                have.add(key)
                at_pos.setdefault(key[0], set()).add((key[1], key[2]))
    finally:
        bam.close()
    if not add:
        return recs, 0
    return recs + add, len(add)


def apply(recs, bam_path, fasta, chrom, orig_bam, donor, workdir, keys, p_any, p_hom, models, cut, resolver=None,
          resolver_cut=0.7):
    """recs (assemble_records output) with disputed clusters replaced where the arbiter says so, and (with a resolver)
    shared-indel 1-bp calls set hom -> (recs, n replaced)"""
    chooser, arbiter = models
    to_hom = set()
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
        if resolver is not None:
            ref_ = fa.fetch(chrom, s - 1, e).upper()
            for k_, f_ in shared_indel_rows(r, ch, ref_, fa, chrom):
                if resolver.predict_proba(np.array([f_]))[0, 1] >= resolver_cut: to_hom.add(k_)
        if d == r["ours"]: continue
        if arbiter.predict_proba(np.array([arb_feats(r, ch, d, conf, rmin)]))[0, 1] < cut: continue
        done.append((s, e, d)); n += 1
    if to_hom:
        recs = [(x[0], x[1], x[2], "hom") + tuple(x[4:]) if x[4] and (x[0], x[1], x[2]) in to_hom else x for x in recs]
    if not done: return recs, 0
    out = [x for x in recs if not (x[4] and any(s <= x[0] and x[0] + len(x[1]) - 1 <= e for s, e, _ in done))]
    have = {(x[0], x[1], x[2]) for x in out}
    for s, e, d in done:
        out.extend((k[0], k[1], k[2], z, True, 1.0) for k, z in d.items() if k not in have)
    return out, n
