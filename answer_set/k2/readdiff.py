"""Recover inserted mutations by comparing each round read with its ORIGINAL read.

Minos builds a GIAB round BAM from a public sequencing run of the person (giab_origin.py): reads
are trimmed to 128 bp, renamed and re-qualified, and some are EDITED - to carry the inserted
mutations, and (since ~2026-09-12) also with scattered random letter changes. Every unedited
round read is an exact piece of an original read, so:

  round read == piece of an original read   -> unedited
  round read ~= piece of an original read   -> the differences are Minos's edits

Sequencing errors of the original run are in the original too and never show up. Random edits
land in one read each; an inserted mutation shows up in many reads at the same place. So for each
edit we report how many reads carry it and how many matched reads cover it.
"""
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import edlib
import parasail
import pysam

K = 16
END_SKIP = 8          # edits this close to a round read's end are often cut-off pieces of an indel
_MATRIX = parasail.matrix_create("ACGTN", 1, -4)


def _edits_affine(q, s, parent):
    """Edits turning original read s into round read q, aligned with affine gaps (opening a gap
    costs more than extending it, so one indel is not split into pieces), skipping edits within
    END_SKIP letters of the round read's ends. (ref_pos_1based, ref_allele, alt_allele)."""
    tb = parasail.sg_dx_trace_scan_16(q, s, 6, 1, _MATRIX).get_traceback("|", ".", " ")
    qa, ta = tb.query, tb.ref
    qmap = dict(parent.get_aligned_pairs(matches_only=True))
    out, ti, qi, i, n = [], 0, 0, 0, len(qa)
    while i < n and qa[i] == "-":                    # original read before the round piece
        ti += 1
        i += 1
    last = n
    while last > i and qa[last - 1] == "-":
        last -= 1
    while i < last:
        if qa[i] == ta[i]:
            ti += 1
            qi += 1
            i += 1
            continue
        j = i
        while j < last and qa[j] != ta[j]:
            j += 1
        o_part, r_part = ta[i:j].replace("-", ""), qa[i:j].replace("-", "")
        q_end = qi + len(r_part)
        if qi >= END_SKIP and q_end <= len(q) - END_SKIP:
            if len(o_part) == len(r_part):
                for k in range(len(o_part)):
                    rp = qmap.get(ti + k)
                    if rp is not None and o_part[k] != r_part[k]:
                        out.append((rp + 1, o_part[k], r_part[k]))
            elif (ti - 1) in qmap:
                base = s[ti - 1]
                out.append((qmap[ti - 1] + 1, base + o_part, base + r_part))
        ti += len(o_part)
        qi = q_end
        i = j
    return out


def _edits(q, s, loc, parent):
    """Edits turning original read s (aligned piece s[loc0:loc1+1]) into round read q, as
    (ref_pos_1based, ref_allele, alt_allele) using the ORIGINAL read's alignment."""
    s0, s1 = loc
    seg = s[s0:s1 + 1]
    al = edlib.align(q, seg, mode="NW", task="path")
    nice = edlib.getNiceAlignment(al, q, seg)
    qa, ta = nice["query_aligned"], nice["target_aligned"]
    qmap = dict(parent.get_aligned_pairs(matches_only=True))
    out, ti, i = [], s0, 0
    while i < len(qa):
        if qa[i] == ta[i]:
            ti += 1
            i += 1
            continue
        j = i
        while j < len(qa) and qa[j] != ta[j]:
            j += 1
        o_part, r_part = ta[i:j].replace("-", ""), qa[i:j].replace("-", "")
        if len(o_part) == len(r_part):
            for k in range(len(o_part)):
                rp = qmap.get(ti + k)
                if rp is not None:
                    out.append((rp + 1, o_part[k], r_part[k]))
        elif (ti - 1) in qmap:
            base = s[ti - 1]
            out.append((qmap[ti - 1] + 1, base + o_part, base + r_part))
        ti += len(o_part)
        i = j
    return out


def _window(args):
    round_bam, orig_bam, chrom, w0, w1 = args
    rb = pysam.AlignmentFile(str(round_bam))
    ob = pysam.AlignmentFile(str(orig_bam))
    by = defaultdict(list)                                   # (strand, start//50) -> originals
    for p in ob.fetch(chrom, max(0, w0 - 400), w1 + 400):
        if p.is_unmapped or p.is_secondary or p.is_supplementary or not p.query_sequence:
            continue                              # unmapped: placed at the mate, no alignment to map edits through
        by[(p.is_reverse, p.reference_start // 50)].append(p)
    support, spans, stats = defaultdict(int), [], defaultdict(int)
    origins = defaultdict(set)                                # edit -> original reads it came from
    for r in rb.fetch(chrom, max(0, w0 - 1), w1):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or not r.query_sequence:
            continue                              # unmapped (placed at the mate, reference_end None): ~40-200 per round
        if not (w0 - 1 <= r.reference_start < w1):
            continue
        q = r.query_sequence
        b = r.reference_start // 50
        cands = by.get((r.is_reverse, b - 1), []) + by.get((r.is_reverse, b), []) + by.get((r.is_reverse, b + 1), [])
        if any(q in p.query_sequence for p in cands):
            stats["unedited"] += 1
            spans.append((r.reference_start, r.reference_end))
            continue
        kmers = {q[i:i + K] for i in range(0, len(q) - K + 1, 8)}
        best = None
        for p in cands:
            s = p.query_sequence
            if sum(1 for k in kmers if k in s) < 2:
                continue
            al = edlib.align(q, s, mode="HW", task="locations", k=12)
            d = al["editDistance"]
            if d < 0:
                continue
            if best is None or d < best[0]:
                best = (d, p, al["locations"][0])
        if best is None:
            stats["no original found"] += 1
            continue
        stats["edited"] += 1
        spans.append((r.reference_start, r.reference_end))
        for e in _edits_affine(q, best[1].query_sequence, best[1]):
            support[e] += 1
            origins[e].add(best[1].query_name)
    return dict(support), spans, dict(stats), {k: len(v) for k, v in origins.items()}


def recover(round_bam, orig_bam, chrom, start, end, workers=8, step=250_000):
    """{(pos, ref, alt): (reads with this edit, matched reads covering it, distinct original reads
    those reads came from)} and read stats. Alleles are as seen in the reads (not yet normalised).
    A random change Minos made once in one original read can reach several round reads cut from
    it; a real inserted mutation comes from many different originals."""
    jobs = [(str(round_bam), str(orig_bam), chrom, w, min(end, w + step - 1)) for w in range(start, end + 1, step)]
    support, stats = defaultdict(int), defaultdict(int)
    starts, ends = [], []
    orig_n = {}
    with ProcessPoolExecutor(workers) as ex:
        for sup, spans, st, on in ex.map(_window, jobs):
            orig_n.update(on)
            for k, v in sup.items():
                support[k] += v
            for k, v in st.items():
                stats[k] += v
            for s, e in spans:
                starts.append(s)
                ends.append(e)
    import bisect
    starts.sort()
    ends.sort()
    out = {}
    for k, n in support.items():
        p0 = k[0] - 1
        # reads covering p0 = started at/before p0 minus ended at/before p0
        cov = bisect.bisect_right(starts, p0 - 1) - bisect.bisect_right(ends, p0 + len(k[1]))
        out[k] = (n, max(n, cov), orig_n.get(k, n))
    return out, dict(stats)


def normalised(found, chrom, workdir):
    """Left-align/trim the edits the same way as the candidates (bcftools norm, REF set from the
    reference). {key: (reads with the edit, matched reads covering it, distinct originals)}."""
    from answer_set.k2.util import bcftools, ref_path, trim
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    raw = workdir / "readdiff.vcf"
    with open(raw, "w") as fh:
        fh.write(f"##fileformat=VCFv4.2\n##contig=<ID={chrom}>\n"
                 '##INFO=<ID=SUP,Number=1,Type=Integer,Description="reads with the edit">\n'
                 '##INFO=<ID=COV,Number=1,Type=Integer,Description="matched reads covering it">\n'
                 '##INFO=<ID=ORI,Number=1,Type=Integer,Description="distinct original reads">\n'
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for (p, r, a), v in sorted(found.items()):
            s, c, o = (tuple(v) + (v[0],))[:3]
            if r != a and set(r + a) <= set("ACGT"):
                fh.write(f"{chrom}\t{p}\t.\t{r}\t{a}\t.\t.\tSUP={s};COV={c};ORI={o}\n")
    out = workdir / "readdiff.norm.vcf.gz"
    bcftools(f"bcftools norm -c s -m -both -f /r/{chrom}.fa /o/{raw.name} -Oz -o /o/{out.name} 2>/dev/null",
             [(str(workdir), "/o"), (str(ref_path(chrom).parent), "/r")])
    res = {}
    for rec in pysam.VariantFile(str(out)):
        if not rec.alts or rec.alts[0] == rec.ref:
            continue
        k = trim(rec.pos, rec.ref, rec.alts[0])
        s, c = rec.info["SUP"], rec.info["COV"]
        o = rec.info.get("ORI", s)
        if k in res:
            s, c, o = s + res[k][0], max(c, res[k][1]), o + res[k][2]
        res[k] = (s, c, o)
    return res
