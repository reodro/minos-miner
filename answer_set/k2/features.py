"""Turn a round BAM into candidate alleles with features (no variant calling involved). Owner: answer_set."""
import bisect
import os
import shutil
import tempfile
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pysam

from answer_set.k2 import masks
from answer_set.k2.util import (DONORS, Bed, bcftools, donor_bed, load_donor_alleles, normalize_vcf,
                                 parse_region, ref_path, retry, trim)

# GIAB read-comparison-only candidates: a 1-bp indel needs this share of its strongest neighbour's support
# (1.0 = must be the strongest, as for g005). Each GIAB model stores its value (meta rd_1bp_factor); build.py
# sets it from the model before extracting. Env MINOS_RD1BP overrides the default used by prepare.py.
RD_1BP_FACTOR = float(os.environ.get("MINOS_RD1BP", "0.5"))
# GIAB read comparison (Google originals + read-vs-original diff): tries before the build goes on without it.
# Without it the model sees NaN for every rd_* measurement - still scorable, but weaker (2026-10-01 08:17 on the
# k2/107 machine: inserted 1-bp indels 94.9% found instead of 96.9%).
RD_ATTEMPTS, RD_WAIT_S = 3, 10
RD_BUDGET_S = 900               # the whole read comparison, every try: past it the build goes on without it
# GIAB local haplotype genotyper (hapgeno.py): the two local DNA sequences the person most likely carries at each
# spot, as events (new candidates) and per-candidate posteriors. On only for models trained with it (meta
# hap_genotyper; build.py sets it from the model before extracting), so a model trained without it is never fed
# candidates it has not seen. prepare.py turns it on (--no-hap: off).
HAP_GENOTYPER = int(os.environ.get("MINOS_HAPGENO", "0"))     # 0 off, else the hapgeno.VERSION it was trained with
# 1000G rounds: compare each round read with the person's public original read (kg_origin.py), as for GIAB. Like
# HAP_GENOTYPER it changes the candidates, so it follows the 1000G model (meta kg_readdiff; build.py sets it).
KG_READDIFF = int(os.environ.get("MINOS_KG_RD", "0"))
# 1000G rounds: the local haplotype genotyper too (hapgeno.py, against the person's public originals); follows the
# 1000G model (meta kg_hapgeno = the hapgeno.VERSION it was trained with; build.py sets it). Needs KG_READDIFF.
KG_HAPGENO = int(os.environ.get("MINOS_KG_HG", "0"))
HG_CALL = 0.5                   # genotyper events at least this likely become candidates,
HG_EDIT_ORI = 2                 # and so do events carried by edited reads from at least this many originals
HG_FEATS = ("hg_post", "hg_hom", "hg_reads", "hg_anchored", "hg_lr", "hg_n_haps", "hg_edit", "hg_edit_ori", "hg_unedit",
            "hg_win_edit_other", "hg_edit_q")

FEATURES = [
    "is_indel", "indel_len", "abs_len", "dp", "ref_n", "alt_n", "alt_frac", "alt_frac_ra",
    "alt_fwd_frac", "strand_min", "ref_fwd_frac", "qs_alt", "alt_bq", "ref_bq", "alt_mq",
    "ref_mq", "alt_tail", "ref_tail", "vdb", "sgb", "rpbz", "mqbz", "mqsbz", "bqbz", "scbz",
    "mq0f", "scr_frac", "idv", "imf", "n_alts_site", "hp_len", "hp_len_after", "str_period2",
    "gc50", "dens50", "dens500", "strong50", "donor_count", "bed_count", "dp_rel",
    "in_mode_donor", "mode_giab",
    # repeat-length profile of reads spanning the repeat tract (indels in repeats only)
    "tr_len", "tr_span", "tr_f_ref", "tr_f_this", "tr_f_other", "tr_f_this_nb", "tr_peak2",
    "tr_this_vs_ref", "tr_n_lengths",
    # 1000 Genomes person (synthetic rounds): public DRAGEN calls and phased panel for that person
    "kg_mode", "in_dragen", "dragen_qual", "dragen_gq", "dragen_hom", "in_panel", "only_dragen",
    # raw 1000 Genomes GATK callset (sites only): Minos keeps only calls GATK also made
    "gatk_site", "gatk_pass", "gatk_ac_log", "dragen_and_gatk",
    # raw GATK call for THIS person, looked up at the person's DRAGEN indels >= 4 bp
    "gatk_person_checked", "gatk_person_called", "gatk_person_same_gt",
    # GIAB-person rounds: Minos inserts mutations only inside the person's confident region
    # (497,453/497,453 in 81 rounds); noise often sits right next to the person's own variants
    "in_donor_bed", "donor_dist_log", "donor_indel_near",
    # 1000 Genomes person, long indels: Minos drops calls where DRAGEN and GATK disagree about the
    # spot (whole genotype, or another call within 10-100 bp) or GATK lists many lengths there
    "spot_gt_same", "disagree10", "disagree30", "disagree100", "gatk_alleles_near",
    # public easy/hard-to-read region lists (masks.py): nearly all of Minos's unexplained drops of
    # DRAGEN+GATK calls lie outside the 1000 Genomes strict mask
    "in_strict_mask", "in_pilot_mask", "in_giab_easy", "in_giab_repeat",
    # population annotations of the GATK site, official panel membership, and how GATK wrote this
    # person's record (Minos drops more calls from poor sites and multi-allele/overlapping records)
    "gatk_vqslod", "gatk_an_frac", "gatk_hwe_log", "gatk_exchet_log", "gatk_mq", "gatk_mqranksum",
    "gatk_fs", "gatk_sor", "gatk_site_dp_log", "gatk_inbreeding", "panel_site",
    "gatk_rec_nalt", "gatk_allele_idx", "gatk_rec_overlap",
    # GIAB rounds: each round read compared with its ORIGINAL read (readdiff.py): reads carrying this
    # exact change vs the original, matched reads covering it; and the local haplotype check
    # (localhap.py): which version of the local DNA each read supports
    "rd_support", "rd_cover", "rd_frac",
    "hap_alt", "hap_ref", "hap_n", "hap_frac", "hap_rank", "hap_margin", "hap_cluster",
    # 1000 Genomes person: Minos's key is built from the person's GATK calls (all 576,785 background
    # variants of 43 rounds are GATK calls; 72% of GATK calls DRAGEN lacks are kept), with GATK's het/hom
    # (100%); GATK calls the read scan misses become candidates too
    "only_gatk", "gatk_person_hom",
    # GIAB rounds: changes seen only by the read comparison become candidates (only_rd); how many
    # DIFFERENT original reads carry the change (one random change copied into several round reads
    # cut from the same original looks like a multi-read variant)
    "only_rd", "rd_orig", "rd_orig_frac",
    # read comparison: strongest competing change within 10 bp (broken/shifted versions of one event)
    "rd_sib_max", "rd_sib_ratio",
    # local haplotype genotyper (hapgeno.py): found only by it (only_hg); P(carried), P(hom), reads best explained by a
    # sequence carrying it, reads anchored around the window, log-likelihood margin, candidate sequences in the window
    "only_hg", "hg_post", "hg_hom", "hg_reads", "hg_anchored", "hg_lr", "hg_n_haps",
    # read comparison anywhere in the candidate's genotyped window, whatever the spelling: strongest edit (reads) and
    # the distinct original reads it came from (an inserted mutation is an edit of many originals; donor variants
    # and noise already in the originals are not)
    "hg_rd_sup", "hg_rd_ori",
    # edit-aware genotyper (each read against its Google original): edited reads best explained by a sequence carrying
    # the event, the distinct originals they came from, unedited reads carrying it (the donor's own DNA), and edited
    # reads in the window carrying other sequences (local Minos noise)
    "hg_edit", "hg_edit_ori", "hg_unedit", "hg_win_edit_other", "hg_edit_q",
]


def _tract(seq, i):
    """Repeat tract (start, end) in seq coordinates around the base after index i:
    a homopolymer or a dinucleotide repeat. Returns None if the tract is short."""
    j = i + 1
    if j >= len(seq):
        return None
    c = seq[j]
    a, b = j, j
    while a > 0 and seq[a - 1] == c:
        a -= 1
    while b + 1 < len(seq) and seq[b + 1] == c:
        b += 1
    if b - a + 1 >= 5:
        return a, b
    unit = seq[j:j + 2]
    if len(unit) == 2 and unit[0] != unit[1]:
        a, b = j, j + 1
        while a - 2 >= 0 and seq[a - 2:a] == unit:
            a -= 2
        while seq[b + 1:b + 3] == unit:
            b += 2
        if b - a + 1 >= 8:
            return a, b
    return None


def tract_profiles(bam_path, chrom, lo, seq, indel_keys):
    """For indel candidates inside repeats: distribution of net length change of the
    tract across reads that span it with 3 bp of flank (MAPQ >= 20)."""
    bam = pysam.AlignmentFile(str(bam_path))
    cache = {}
    out = {}
    for k in indel_keys:
        t = _tract(seq, k[0] - lo)
        if t is None:
            continue
        ta, tb = t[0] + lo, t[1] + lo            # 1-based reference coordinates of the tract
        if (ta, tb) not in cache:
            hist = Counter()
            for r in bam.fetch(chrom, ta - 4, tb + 3):
                if r.mapping_quality < 20 or r.is_secondary or r.is_supplementary or r.cigartuples is None:
                    continue
                if r.reference_start > ta - 4 or r.reference_end < tb + 3:
                    continue
                if r.cigartuples[0][0] == 4 and r.reference_start > ta - 10:
                    continue
                pos, delta = r.reference_start, 0      # 0-based
                for op, L in r.cigartuples:
                    if op == 1 and ta - 2 <= pos <= tb:
                        delta += L
                    elif op == 2 and ta - 2 <= pos + L - 1 and pos <= tb:
                        delta -= L
                    if op in (0, 2, 3, 7, 8):
                        pos += L
                hist[delta] += 1
            cache[(ta, tb)] = hist
        hist = cache[(ta, tb)]
        n = sum(hist.values())
        this = len(k[2]) - len(k[1])
        if n == 0:
            out[k] = {"tr_len": tb - ta + 1, "tr_span": 0}
            continue
        f = {d: c / n for d, c in hist.items()}
        top = sorted(f.values(), reverse=True)
        out[k] = {
            "tr_len": tb - ta + 1, "tr_span": n, "tr_f_ref": f.get(0, 0.0), "tr_f_this": f.get(this, 0.0),
            "tr_f_other": 1 - f.get(0, 0.0) - f.get(this, 0.0),
            "tr_f_this_nb": f.get(this - 1, 0.0) + f.get(this + 1, 0.0) - (f.get(0, 0.0) if abs(this) == 1 else 0.0),
            "tr_peak2": top[1] if len(top) > 1 else 0.0,
            "tr_this_vs_ref": f.get(this, 0.0) / max(1e-9, f.get(this, 0.0) + f.get(0, 0.0)),
            "tr_n_lengths": sum(1 for v in f.values() if v >= 0.05),
        }
    return out


def run_candidates(bam, region, out_vcf, threads=8):
    """Every position where reads show a non-reference allele (>= 2 reads), with per-allele counts, strand
    counts and bcftools bias statistics. No calling. The region is split into `threads` equal pieces run in
    parallel, then concatenated in order."""
    chrom, s, e = parse_region(region)
    bam, ref, out_vcf = Path(bam).resolve(), ref_path(chrom).resolve(), Path(out_vcf).resolve()
    step = (e - s + threads) // threads
    chunks = [(i, s + i * step, min(e, s + i * step + step - 1)) for i in range(threads) if s + i * step <= e]
    tmp = Path(tempfile.mkdtemp(prefix=".cand.", dir=out_vcf.parent))
    mounts = [(str(bam.parent), "/b"), (str(ref.parent), "/r"), (str(tmp), "/o")]
    try:
        def one(c):
            i, x, y = c
            bcftools(f"bcftools mpileup -f /r/{ref.name} -r {chrom}:{x}-{y} -d 10000 -q 0 -Q 0 "
                     f"-a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/QS,INFO/AD,INFO/SCR -Ou /b/{bam.name} "
                     f'| bcftools view -i "INFO/AD[1]>=2" -Oz -o /o/c{i}.vcf.gz', mounts)

        with ThreadPoolExecutor(len(chunks)) as ex:
            list(ex.map(one, chunks))
        (tmp / "list").write_text("".join(f"/o/c{i}.vcf.gz\n" for i, _, _ in chunks))
        bcftools("bcftools concat -f /o/list -Oz -o /o/all.vcf.gz && bcftools index -t /o/all.vcf.gz", mounts)
        shutil.move(str(tmp / "all.vcf.gz"), str(out_vcf))
        shutil.move(str(tmp / "all.vcf.gz.tbi"), str(out_vcf) + ".tbi")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out_vcf


def _hp(seq, i):
    """Length of the homopolymer run that contains index i, and the run just after i."""
    def run(j):
        if j < 0 or j >= len(seq):
            return 0
        c, a, b = seq[j], j, j
        while a > 0 and seq[a - 1] == c:
            a -= 1
        while b + 1 < len(seq) and seq[b + 1] == c:
            b += 1
        return b - a + 1
    return run(i), run(i + 1)


def _str2(seq, i):
    """Number of copies of the dinucleotide starting at i+1 repeated in tandem."""
    unit = seq[i + 1:i + 3]
    if len(unit) < 2 or unit[0] == unit[1]:
        return 0
    n, j = 0, i + 1
    while seq[j:j + 2] == unit:
        n, j = n + 1, j + 2
    return n


def detect_donor(cands, donors):
    """Which GIAB donor (if any) is this round's background?

    For each donor, the share of its alleles in the window that the reads support.
    A pure-donor round supports ~99% of that donor's alleles; others share far less.
    """
    supported = {k for k, f in cands.items() if f["alt_n"] >= 2 and f["alt_frac_ra"] >= 0.1}
    scores = {}
    for d, alle in donors.items():
        snps = [k for k in alle if len(k[1]) == 1 and len(k[2]) == 1]
        scores[d] = (sum(1 for k in snps if k in supported) / len(snps)) if snps else 0.0
    best = max(scores, key=scores.get)
    return best, scores


def extract(bam, region, workdir, threads=8, cand_vcf=None, reuse=False):
    """Return (keys, X, info, donors, beds) for all candidate alleles in the region."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    chrom, a, b = parse_region(region)
    cached = workdir / "cand.raw.vcf.gz"
    if cand_vcf:
        raw = Path(cand_vcf)
    elif reuse and cached.exists() and cached.stat().st_mtime > Path(bam).stat().st_mtime:
        raw = cached                     # only if made from this exact BAM (not an older/truncated copy)
    else:
        raw = run_candidates(bam, region, cached, threads)
    norm = normalize_vcf(raw, workdir / "cand.norm.vcf.gz", ref_path(chrom))

    fa = pysam.FastaFile(str(ref_path(chrom)))
    lo = max(1, a - 1000)
    seq = fa.fetch(chrom, lo - 1, b + 1000).upper()

    cands = {}
    site_alts = Counter()
    for rec in pysam.VariantFile(str(norm)):
        alt = rec.alts[0] if rec.alts else None
        if alt is None or alt.startswith("<") or alt == "*":
            continue
        s = rec.samples[0]
        ad = s.get("AD") or (0, 0)
        if len(ad) < 2 or ad[1] is None or ad[1] < 2:
            continue
        info = rec.info
        i16 = info.get("I16")
        adf, adr = s.get("ADF") or (0, 0), s.get("ADR") or (0, 0)
        qs = s.get("QS") or (0, 0)
        dp = info.get("DP") or 0
        ref_n, alt_n = ad[0] or 0, ad[1] or 0
        key = trim(rec.pos, rec.ref, alt)
        L = len(key[2]) - len(key[1])
        f = {
            "is_indel": float(L != 0), "indel_len": float(L), "abs_len": float(abs(L)),
            "dp": float(dp), "ref_n": float(ref_n), "alt_n": float(alt_n),
            "alt_frac": alt_n / dp if dp else 0.0,
            "alt_frac_ra": alt_n / (ref_n + alt_n) if ref_n + alt_n else 0.0,
            "alt_fwd_frac": (adf[1] or 0) / alt_n if alt_n else 0.0,
            "strand_min": min(adf[1] or 0, adr[1] or 0) / alt_n if alt_n else 0.0,
            "ref_fwd_frac": (adf[0] or 0) / ref_n if ref_n else 0.5,
            "qs_alt": float(qs[1]) / max(1e-9, sum(x or 0 for x in qs)) if qs else 0.0,
            "vdb": info.get("VDB", np.nan), "sgb": info.get("SGB", np.nan),
            "rpbz": info.get("RPBZ", np.nan), "mqbz": info.get("MQBZ", np.nan),
            "mqsbz": info.get("MQSBZ", np.nan), "bqbz": info.get("BQBZ", np.nan),
            "scbz": info.get("SCBZ", np.nan), "mq0f": info.get("MQ0F", np.nan),
            "scr_frac": (info.get("SCR") or 0) / dp if dp else 0.0,
            "idv": float(info.get("IDV") or 0), "imf": float(info.get("IMF") or 0),
        }
        if i16 and len(i16) >= 16:
            nr, na = (i16[0] + i16[1]) or 1, (i16[2] + i16[3]) or 1
            f.update(ref_bq=i16[4] / nr, alt_bq=i16[6] / na, ref_mq=i16[8] / nr,
                     alt_mq=i16[10] / na, ref_tail=i16[12] / nr, alt_tail=i16[14] / na)
        else:
            f.update(ref_bq=np.nan, alt_bq=np.nan, ref_mq=np.nan, alt_mq=np.nan,
                     ref_tail=np.nan, alt_tail=np.nan)
        i = key[0] - lo
        f["hp_len"], f["hp_len_after"] = _hp(seq, i)
        f["str_period2"] = _str2(seq, i)
        w = seq[max(0, i - 50): i + 51]
        f["gc50"] = (w.count("G") + w.count("C")) / max(1, len(w))
        site_alts[key[0]] += 1
        cands[key] = f

    # neighbourhood density
    pos = sorted(k[0] for k in cands)
    strong = sorted(k[0] for k, f in cands.items() if f["alt_frac_ra"] >= 0.2 and f["alt_n"] >= 3)
    dps = [f["dp"] for f in cands.values()]
    med_dp = float(np.median(dps)) if dps else 1.0

    def count(arr, p, w):
        return bisect.bisect_right(arr, p + w) - bisect.bisect_left(arr, p - w) - 1

    donors = load_donor_alleles(chrom, a, b)
    beds = {d: Bed(donor_bed(d), chrom) for d in DONORS}
    mode_donor, donor_scores = detect_donor(cands, donors)
    mode_giab = donor_scores[mode_donor] >= 0.95

    # synthetic round: find the 1000 Genomes person and bring in their DRAGEN + panel calls
    kg_info, dragen, panel_p = {}, {}, {}
    gatk_person, gp_checked, gp_layout, psites = {}, set(), {}, set()
    if not mode_giab:
        from answer_set.k2 import kg
        if kg.index_path(chrom).exists():
            frac = {k: f["alt_frac_ra"] for k, f in cands.items()}
            person, sc, sc2, n_sites = kg.identify(chrom, a, b, frac)
            kg_info = {"person": person, "agree": sc, "next": sc2, "sites": n_sites}
            if sc >= 0.9:
                try:
                    dragen = kg.dragen_alleles(chrom, a, b, person, workdir)
                    panel_p = kg.sample_alleles(chrom, a, b, person, workdir)
                except RuntimeError as e:          # network trouble: say so loudly, fall back
                    kg_info["error"] = str(e)
                    dragen, panel_p = {}, {}
                gp_loci = [k for k in dragen if a <= k[0] <= b]      # every DRAGEN call of the person
                if dragen:
                    try:
                        gatk_person = kg.gatk_person_region(chrom, a, b, person)
                        gp_layout = kg.gatk_person_layout(chrom, a, b, person)
                        gp_checked = set(gp_loci)
                        kg_info["gatk_person_checked"] = len(gp_loci)
                    except RuntimeError as e:
                        kg_info["gatk_person_error"] = str(e)
                for k, (z, q, gq) in dragen.items():
                    if k not in cands and a <= k[0] <= b:     # e.g. long indels the pile-up scan misses
                        L = len(k[2]) - len(k[1])
                        f = {n: np.nan for n in FEATURES}
                        f.update(is_indel=float(L != 0), indel_len=float(L), abs_len=float(abs(L)),
                                 alt_n=0.0, only_dragen=1.0)
                        i = k[0] - lo
                        f["hp_len"], f["hp_len_after"] = _hp(seq, i)
                        f["str_period2"] = _str2(seq, i)
                        w = seq[max(0, i - 50): i + 51]
                        f["gc50"] = (w.count("G") + w.count("C")) / max(1, len(w))
                        cands[k] = f
                for k in gatk_person:
                    if k not in cands and a <= k[0] <= b:     # GATK calls neither the scan nor DRAGEN has
                        L = len(k[2]) - len(k[1])
                        f = {n: np.nan for n in FEATURES}
                        f.update(is_indel=float(L != 0), indel_len=float(L), abs_len=float(abs(L)),
                                 alt_n=0.0, only_dragen=0.0, only_gatk=1.0)
                        i = k[0] - lo
                        f["hp_len"], f["hp_len_after"] = _hp(seq, i)
                        f["str_period2"] = _str2(seq, i)
                        w = seq[max(0, i - 50): i + 51]
                        f["gc50"] = (w.count("G") + w.count("C")) / max(1, len(w))
                        cands[k] = f

    from answer_set.k2 import kg as _kg
    gsites = _kg.gatk_sites(chrom, a, b, workdir) or {}
    if dragen:
        psites = _kg.panel_sites(chrom, a, b, workdir)
    dbed = beds[mode_donor] if mode_giab else None
    dlist = donors[mode_donor] if mode_giab else {}
    dpos = sorted(k[0] for k in dlist)
    dind = sorted(k[0] for k in dlist if len(k[1]) != len(k[2]))

    def count_range(arr, lo_, hi_):
        return bisect.bisect_right(arr, hi_) - bisect.bisect_left(arr, lo_)

    def nearest(arr, p):
        j = bisect.bisect_left(arr, p)
        return min([abs(arr[x] - p) for x in (j - 1, j) if 0 <= x < len(arr)] or [10 ** 6])
    # DRAGEN vs GATK disagreements for this person (a call in one but not the other, or other het/hom)
    dz = {k: v[0] for k, v in dragen.items()}
    dis = sorted(k[0] for k in set(dz) | set(gatk_person) if dz.get(k) != gatk_person.get(k))
    gpos = sorted(k[0] for k in gsites)

    dz_pos = sorted(dz)
    gp_pos = sorted(gatk_person)

    def near_calls(calls, sorted_keys, p, e):
        i = bisect.bisect_left(sorted_keys, (p - 1,))
        out = {}
        while i < len(sorted_keys) and sorted_keys[i][0] <= e + 1:
            out[sorted_keys[i]] = calls[sorted_keys[i]]
            i += 1
        return set(out.items())

    rdiff, hap, rd_info = {}, {}, {}
    hg, hg_window = {}, None
    rd_on = mode_giab or bool(KG_READDIFF and dragen)          # 1000G: once the person is known
    if rd_on:
        from answer_set.k2 import giab_origin, kg_origin, localhap, readdiff

        rd_end = time.monotonic() + RD_BUDGET_S

        def compare():
            left = max(0.0, rd_end - time.monotonic())
            orig = (giab_origin.fetch(mode_donor, chrom, a, b, budget_s=left) if mode_giab else
                    kg_origin.fetch(kg_info["person"], chrom, a, b, budget_s=left))
            found, stats = readdiff.recover(bam, orig, chrom, a, b, workers=threads)
            return readdiff.normalised(found, chrom, workdir), stats
        try:
            (rdiff, rd_info), failed = retry(compare, RD_ATTEMPTS, RD_WAIT_S, budget_s=RD_BUDGET_S)
            if failed:
                rd_info["failed_tries"] = failed
        except RuntimeError as e:               # no originals after every try: the model sees NaN
            rd_info = {"error": str(e)[:1000]}
        own = donors[mode_donor] if mode_giab else {}
        # changes the read comparison sees but the read scan missed (mostly indels Minos inserted:
        # 74% of never-candidate insertions are visible here) become candidates too
        rd_keys = sorted(rdiff)
        rd_pos = [k[0] for k in rd_keys]

        def sib_max(k):
            """strongest OTHER read-comparison change within 10 bp (a broken version of the same event)"""
            j0, j1 = bisect.bisect_left(rd_pos, k[0] - 10), bisect.bisect_right(rd_pos, k[0] + len(k[1]) + 10)
            return max([rdiff[o][0] for o in rd_keys[j0:j1] if o != k] or [0])
        n_rd_only = 0
        for k, (sup, cov, ori) in rdiff.items():
            # 1-bp indels: half the strongest neighbour is enough (g005's "must be strongest" rule lost real ones)
            need = RD_1BP_FACTOR * sib_max(k) if abs(len(k[1]) - len(k[2])) == 1 else sib_max(k)
            if (k in cands or k in own or not (a <= k[0] <= b) or sup < 2 or sup < 0.25 * cov or sup < need
                    or (mode_giab and not (dbed.contains(k[0]) and dbed.contains(k[0] + len(k[1]) - 1)))):
                continue
            L = len(k[2]) - len(k[1])
            f = {n: np.nan for n in FEATURES}
            f.update(is_indel=float(L != 0), indel_len=float(L), abs_len=float(abs(L)), alt_n=0.0, only_rd=1.0)
            i = k[0] - lo
            f["hp_len"], f["hp_len_after"] = _hp(seq, i)
            f["str_period2"] = _str2(seq, i)
            w = seq[max(0, i - 50): i + 51]
            f["gc50"] = (w.count("G") + w.count("C")) / max(1, len(w))
            cands[k] = f
            n_rd_only += 1
        rd_info["candidates_added"] = n_rd_only
        if (HAP_GENOTYPER and mode_giab) or (KG_HAPGENO and not mode_giab):
            from answer_set.k2 import hapcand, hapgeno
            seeds = [(k[0], k[0] + len(k[1])) for k, f in cands.items()
                     if len(k[1]) != len(k[2]) or (f.get("alt_n", 0) >= 2 and f.get("alt_frac_ra", 0) >= 0.08)]
            seeds += [(k[0], k[0] + len(k[1])) for k, v in rdiff.items() if v[0] >= 2]
            orig_bam = (giab_origin.path(mode_donor, chrom, a, b) if mode_giab else
                        kg_origin.path(kg_info["person"], chrom, a, b))
            hg = hapgeno.genotype(str(bam), str(ref_path(chrom)), chrom, seeds, threads, known=set(cands) | set(rdiff),
                                  orig_bam=orig_bam if "error" not in rd_info and orig_bam.exists() else None)
            from answer_set.k2 import hapwin
            hapwin.save_windows(workdir / "hg_windows.pkl.gz", hg)        # for window decisions (hapwin.py)
            wins = sorted((k[1], k[2], v) for k, v in hg.items() if k[0] == "window")
            w_starts = [w[0] for w in wins]

            def hg_window(pos):
                """the genotyped window around pos (start, end, summary), or None"""
                i = bisect.bisect_right(w_starts, pos) - 1
                return wins[i] if i >= 0 and wins[i][0] <= pos <= wins[i][1] else None
            n_hg_only = 0
            for k, v in hg.items():
                if (k[0] == "window" or (v["hg_post"] < HG_CALL and not v["hg_edit_ori"] >= HG_EDIT_ORI)
                        or k in cands or k in own or not (a <= k[0] <= b)
                        or abs(len(k[1]) - len(k[2])) > hapcand.MAX_INDEL
                        or (mode_giab and not (dbed.contains(k[0]) and dbed.contains(k[0] + len(k[1]) - 1)))):
                    continue
                L = len(k[2]) - len(k[1])
                f = {n: np.nan for n in FEATURES}
                f.update(is_indel=float(L != 0), indel_len=float(L), abs_len=float(abs(L)), alt_n=0.0, only_hg=1.0)
                i = k[0] - lo
                f["hp_len"], f["hp_len_after"] = _hp(seq, i)
                f["str_period2"] = _str2(seq, i)
                w = seq[max(0, i - 50): i + 51]
                f["gc50"] = (w.count("G") + w.count("C")) / max(1, len(w))
                cands[k] = f
                n_hg_only += 1
            rd_info["hg_candidates_added"] = n_hg_only
            rd_info["hg_windows"] = len(wins)
    keys = list(cands)
    prof = tract_profiles(bam, chrom, lo, seq, [k for k in keys if len(k[1]) != len(k[2])])
    if mode_giab:
        kpos = sorted(k[0] for k in set(keys) | set(own))

        def crowded(k):
            return count(kpos, k[0], 10) > 1
        only = [k for k in keys if k not in own and (len(k[1]) != len(k[2]) or crowded(k))]
        hap = localhap.hap_support(bam, ref_path(chrom), chrom, set(keys) | {k for k in own if a <= k[0] <= b},
                                   only=only)
    for k in keys:
        cands[k].update(prof.get(k, {}))
        f = cands[k]
        f["n_alts_site"] = site_alts[k[0]]
        f["dens50"], f["dens500"] = count(pos, k[0], 50), count(pos, k[0], 500)
        f["strong50"] = count(strong, k[0], 50)
        f["donor_count"] = sum(1 for d in DONORS if k in donors[d])
        f["bed_count"] = sum(1 for d in DONORS if beds[d].contains(k[0]))
        f["dp_rel"] = f["dp"] / med_dp if med_dp else 1.0
        f["in_mode_donor"] = float(mode_giab and k in donors[mode_donor])
        f["mode_giab"] = float(mode_giab)
        d = dragen.get(k)
        f["kg_mode"] = float(bool(dragen))
        f["in_dragen"] = float(d is not None)
        f["dragen_qual"] = d[1] if d else 0.0
        f["dragen_gq"] = d[2] if d else 0.0
        f["dragen_hom"] = float(d is not None and d[0] == "hom")
        f["in_panel"] = float(k in panel_p)
        f.setdefault("only_dragen", 0.0)
        gs = gsites.get(k)
        f["gatk_site"] = float(gs is not None)
        f["gatk_pass"] = float(gs is not None and gs[0])
        f["gatk_ac_log"] = float(np.log1p(gs[1])) if gs else 0.0
        f["dragen_and_gatk"] = float(d is not None and gs is not None and gs[0])
        f["gatk_person_checked"] = float(k in gp_checked)
        f["gatk_person_called"] = float(k in gatk_person)
        f["gatk_person_same_gt"] = float(d is not None and gatk_person.get(k) == d[0])
        f.setdefault("only_gatk", 0.0 if dragen else np.nan)
        gz = gatk_person.get(k)
        f["gatk_person_hom"] = (1.0 if gz == "hom" else 0.0) if gz else np.nan
        if dragen:
            e = k[0] + len(k[1]) - 1
            own = int(dz.get(k) != gatk_person.get(k))   # this call's own disagreement is not counted
            f["spot_gt_same"] = float(near_calls(dz, dz_pos, k[0], e) == near_calls(gatk_person, gp_pos, k[0], e))
            for w in (10, 30, 100):
                f[f"disagree{w}"] = float(count_range(dis, k[0] - w, e + w) - own)
            f["gatk_alleles_near"] = float(count_range(gpos, k[0] - 2, e + 2))
        else:
            f["spot_gt_same"] = f["disagree10"] = f["disagree30"] = f["disagree100"] = np.nan
            f["gatk_alleles_near"] = np.nan
        for m_ in masks.MASKS:
            f[m_] = masks.inside(m_, chrom, k[0], k[1]) if dragen else np.nan
        if dragen:
            ann = gs[2] if gs else {}
            g_ = lambda n: ann.get(n, np.nan)
            f["gatk_vqslod"] = g_("VQSLOD")
            f["gatk_an_frac"] = g_("AN") / 6404.0
            f["gatk_hwe_log"] = -np.log10(max(g_("HWE"), 1e-300)) if not np.isnan(g_("HWE")) else np.nan
            f["gatk_exchet_log"] = -np.log10(max(g_("ExcHet"), 1e-300)) if not np.isnan(g_("ExcHet")) else np.nan
            f["gatk_mq"], f["gatk_mqranksum"] = g_("MQ"), g_("MQRankSum")
            f["gatk_fs"], f["gatk_sor"] = g_("FS"), g_("SOR")
            f["gatk_site_dp_log"] = np.log1p(g_("DP")) if not np.isnan(g_("DP")) else np.nan
            f["gatk_inbreeding"] = g_("InbreedingCoeff")
            f["panel_site"] = float(k in psites)
            lay = gp_layout.get(k)
            f["gatk_rec_nalt"], f["gatk_allele_idx"], f["gatk_rec_overlap"] = (map(float, lay) if lay else (0.0, 0.0, 0.0))
        else:
            for n_ in ("gatk_vqslod", "gatk_an_frac", "gatk_hwe_log", "gatk_exchet_log", "gatk_mq", "gatk_mqranksum",
                       "gatk_fs", "gatk_sor", "gatk_site_dp_log", "gatk_inbreeding", "panel_site",
                       "gatk_rec_nalt", "gatk_allele_idx", "gatk_rec_overlap"):
                f[n_] = np.nan
        if mode_giab:
            f["in_donor_bed"] = float(dbed.contains(k[0]) and dbed.contains(k[0] + len(k[1]) - 1))
            f["donor_dist_log"] = float(np.log1p(nearest(dpos, k[0])))
            f["donor_indel_near"] = float(nearest(dind, k[0]) <= 10)
        else:
            f["in_donor_bed"] = f["donor_dist_log"] = f["donor_indel_near"] = np.nan
        if rd_on and "error" not in rd_info:
            sup, cov, ori = rdiff.get(k, (0, 0, 0))
            f["rd_support"], f["rd_cover"] = float(sup), float(cov)
            f["rd_frac"] = sup / cov if cov else 0.0
            f["rd_orig"] = float(ori)
            f["rd_orig_frac"] = ori / sup if sup else np.nan
            f.setdefault("only_rd", 0.0)
            if np.isnan(f["only_rd"]):
                f["only_rd"] = 0.0
            sm = sib_max(k)
            f["rd_sib_max"] = float(sm)
            f["rd_sib_ratio"] = sup / (sup + sm) if sup + sm else np.nan
        for n_ in ("rd_support", "rd_cover", "rd_frac", "rd_orig", "rd_orig_frac", "only_rd", "rd_sib_max",
                   "rd_sib_ratio") \
                if not rd_on or "error" in rd_info else ():
            f[n_] = np.nan
        if hg_window is not None:
            ev_, w_ = hg.get(k), hg_window(k[0])
            if ev_ is not None:                      # an event of a likely pair
                f.update({n_: ev_[n_] for n_ in HG_FEATS})
            elif w_ is not None:                     # genotyped, but no likely pair carries it
                ws_ = w_[2]
                no_orig = np.isnan(ws_["hg_win_edit_other"])
                f.update(hg_post=0.0, hg_hom=0.0, hg_reads=0.0, hg_anchored=ws_["hg_anchored"], hg_lr=np.nan,
                         hg_n_haps=ws_["hg_n_haps"], hg_edit=np.nan if no_orig else 0.0,
                         hg_edit_ori=np.nan if no_orig else 0.0, hg_unedit=np.nan if no_orig else 0.0,
                         hg_win_edit_other=ws_["hg_win_edit_other"], hg_edit_q=np.nan)
            else:
                for n_ in HG_FEATS:
                    f[n_] = np.nan
            if w_ is not None and "error" not in rd_info:
                j0, j1 = bisect.bisect_left(rd_pos, w_[0]), bisect.bisect_right(rd_pos, w_[1])
                f["hg_rd_sup"] = float(max([rdiff[o][0] for o in rd_keys[j0:j1]] or [0]))
                f["hg_rd_ori"] = float(max([rdiff[o][2] for o in rd_keys[j0:j1]] or [0]))
            else:
                f["hg_rd_sup"] = f["hg_rd_ori"] = np.nan
            if np.isnan(f.get("only_hg", np.nan)):
                f["only_hg"] = 0.0
        else:
            for n_ in ("only_hg", "hg_rd_sup", "hg_rd_ori") + HG_FEATS:
                f[n_] = np.nan
        h_ = hap.get(k)
        for n_ in localhap.FEATS if mode_giab else ("hap_alt", "hap_ref", "hap_n", "hap_frac", "hap_rank",
                                                      "hap_margin", "hap_cluster"):
            f[n_] = float(h_[n_]) if h_ else np.nan
    X = np.array([[cands[k].get(n, np.nan) for n in FEATURES] for k in keys], dtype=float)
    info = {"mode_donor": mode_donor, "mode_giab": bool(mode_giab), "donor_scores": donor_scores,
            "n_candidates": len(keys), "median_dp": med_dp, "kg": kg_info, "n_dragen": len(dragen)}
    info["giab_bed"] = mode_donor if mode_giab else None
    info["readdiff"] = rd_info
    info["kg_readdiff"] = int(bool(KG_READDIFF) and not mode_giab)
    info["kg_hapgeno"] = int(KG_HAPGENO) if (not mode_giab and rd_on and hg_window is not None) else 0
    return keys, X, info, donors, beds
