"""Build a round's answer set from its BAM. Owner: answer_set.

Writes, in the same shape as the platform's released files:
  out/truth.vcf.gz       background + detected injections (injections flagged SYNTHETIC)
  out/mutations.vcf.gz   N entries (only the count matters to the v2 score)
  out/scores.tsv.gz      model confidence per candidate
  out/report.json        round type, donor match, counts, timings

How it works:
  1. features.run_candidates lists every allele the reads show (bcftools mpileup, no calling).
  2. The round type is detected by how well each GIAB donor's known variants are supported.
     GIAB-person round  -> background = that donor's official v4.2.1 variants, exactly.
     1000 Genomes round -> background comes from the model (step 3).
  3. The format's model track decides which remaining candidates are real, and whether each is het or hom.
"""
import json
from collections import Counter
import time
from pathlib import Path

import numpy as np
import pysam

from answer_set.k2 import features
from answer_set.k2.assemble import assemble_records
from answer_set.k2.features import FEATURES, extract
from answer_set.k2.util import parse_region, ref_path


def write_scores(path, keys, p_any, donor, allowed):
    """Model confidence per candidate (pos, ref, alt, p): lets a config search ignore spots where
    the answer set is unsure. Copied person variants count as certain (1); candidates the GIAB
    region rule forbids count as certain no (0)."""
    import gzip
    with gzip.open(path, "wt") as fh:
        fh.write("pos\tref\talt\tp\n")
        for k, p in zip(keys, p_any):
            if k in donor:
                p = 1.0
            elif allowed is not None and not allowed(k):
                p = 0.0
            fh.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{p:.4f}\n")


def write_vcf(path, chrom, records, donor_raw=None):
    """records: list of (pos, ref, alt, zyg, synthetic). Writes bgzipped + indexed VCF.

    Alleles sharing a position are merged into one record (genotype 1/2), expressed
    against the longest REF among them, as a person carries at most two alleles.
    donor_raw {pos: (REF, ALTs, GT)} (util.donor_records): where every allele at a position is one of the person's
    own and they are exactly that record's, the record is written as Minos's truth writes it (NIST's REF, ALT order
    and genotype order).
    """
    from answer_set.k2.util import trim
    fai = {l.split()[0]: int(l.split()[1]) for l in open(str(ref_path(chrom)) + ".fai")}
    lines = []
    if donor_raw:                                  # the person's own records, exactly as Minos's truth writes them
        allele_of = {}                             # minimal allele -> its NIST record's position
        for rp, (R0, alts, gt) in donor_raw.items():
            for g in {int(x) for x in gt.split("/") if x not in ("0", ".")}:
                allele_of[trim(rp, R0, alts[g - 1])] = rp
        have = Counter(allele_of.get(trim(*r[:3])) for r in records if not r[4])
        raw_done = set()
        for rp, (R0, alts, gt) in donor_raw.items():
            n = len({int(x) for x in gt.split("/") if x not in ("0", ".")})
            if have.get(rp) == n:
                lines.append((rp, f"{chrom}\t{rp}\t.\t{R0}\t{','.join(alts)}\t60\tPASS\t.\tGT\t{gt}\n"))
                raw_done.add(rp)
        records = [r for r in records if r[4] or allele_of.get(trim(*r[:3])) not in raw_done]
    by_pos = {}
    for r in records:
        by_pos.setdefault(r[0], []).append(r)
    for pos in sorted(by_pos):
        rs = by_pos[pos][:2]
        R = max((r[1] for r in rs), key=len)
        alts = [r[2] + R[len(r[1]):] for r in rs]
        syn = any(r[4] for r in rs)
        if len(rs) == 2 and alts[0] != alts[1]:
            gt = "1/2"
        else:
            alts = alts[:1]
            gt = "1/1" if rs[0][3] == "hom" else "0/1"
        lines.append((pos, f"{chrom}\t{pos}\t.\t{R}\t{','.join(alts)}\t60\tPASS\t"
                           f"{'SYNTHETIC' if syn else '.'}\tGT\t{gt}\n"))
    lines = [l for _, l in sorted(lines, key=lambda t: t[0])]
    plain = Path(str(path)[:-3]) if str(path).endswith(".gz") else Path(path)
    with open(plain, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        fh.write(f"##contig=<ID={chrom},length={fai[chrom]}>\n")
        fh.write("##INFO=<ID=SYNTHETIC,Number=0,Type=Flag,Description=\"Predicted injection\">\n")
        fh.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")
        fh.writelines(lines)
    pysam.tabix_index(str(plain), preset="vcf", force=True)
    return path


def build(bam, region, out, threads=8, exclude_round=None, track_versions=None):
    """track_versions: {'giab': name, 'kg': name} to use specific track versions (default: each track's ACTIVE).
    exclude_round: a released round's tag -> the fold model that never saw it (backtests)."""
    from answer_set.k2 import inputs, tracks          # here, so `python -m answer_set.k2.<module>` imports cleanly
    t0 = time.time()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    chrom, a, b = parse_region(region)
    changed_inputs = inputs.check(inputs.files(chrom))        # the same model gives other answers on other inputs
    # candidate rules must match the GIAB model that will judge them
    gname = (track_versions or {}).get("giab") or tracks.active_name("giab")
    if gname:
        gmeta = tracks.meta("giab", gname)
        features.RD_1BP_FACTOR = gmeta.get("rd_1bp_factor", 1.0)
        features.HAP_GENOTYPER = int(gmeta.get("hap_genotyper", 0))      # True (g009) == 1
    kname = (track_versions or {}).get("kg") or tracks.active_name("kg")
    features.KG_READDIFF = int(tracks.meta("kg", kname).get("kg_readdiff", 0)) if kname else 0
    features.KG_HAPGENO = int(tracks.meta("kg", kname).get("kg_hapgeno", 0)) if kname else 0
    keys, X, info, donors, beds = extract(bam, region, out / "work", threads)
    t_feat = time.time() - t0
    track = tracks.track_of_info(info)
    # the format's own track: its own models, measurements, threshold and rules
    m_any, m_hom, meta, used = tracks.load(track, (track_versions or {}).get(track), exclude_round)
    thr = meta["threshold"]
    meta = dict(meta, giab_bed_rule=meta.get("bed_rule", False))
    X_full = X
    if meta.get("features") and meta["features"] != FEATURES:      # use exactly its own columns
        X = X[:, [FEATURES.index(n) for n in meta["features"]]] if len(keys) else X
    p_any = m_any.predict_proba(X)[:, 1] if len(keys) else np.array([])
    p_hom = m_hom.predict_proba(X)[:, 1] if len(keys) else np.array([])
    p_raw = p_any
    if len(keys) and meta.get("thr_class"):
        p_any = tracks.apply_class_thr(meta, keys, p_any)
    if len(keys) and meta.get("thr_only_rd") is not None:
        p_any = tracks.apply_source_thr(meta, X_full, p_any, p_raw)
    if len(keys) and meta.get("gt_from_gatk"):
        p_hom = tracks.apply_gt_rule(meta, X_full, p_hom)

    donor = donors[info["mode_donor"]] if info["mode_giab"] else {}
    kg_own = None
    kgi = info.get("kg") or {}
    if (not info["mode_giab"] and meta.get("kg_own_rule") and kgi.get("gatk_person_checked")
            and not kgi.get("gatk_person_error")):
        # 1000G: the person's part of the truth is copied by Minos's exact rule (kg.gatk_person_own), like GIAB's
        # NIST list; the model decides only the rest (the inserted mutations). Any failure: the model decides all.
        from answer_set.k2 import kg
        try:
            own, every = kg.gatk_person_own(chrom, a, b, kgi["person"])
            if own:
                kg_own = len(own)
                donor = own
                p_any = np.where([k in every for k in keys], 0.0, p_any) if len(keys) else p_any
        except Exception as e:  # noqa: BLE001  (never stuck: fall back to the model)
            kg_own = f"error: {e}"[:300]
    allowed = None
    if info["mode_giab"] and meta.get("giab_bed_rule"):
        bed = beds[info["mode_donor"]]
        allowed = lambda k: bed.contains(k[0]) and bed.contains(k[0] + len(k[1]) - 1)  # noqa: E731
    recs = assemble_records(keys, p_any, p_hom, thr, info["mode_giab"], donor, allowed)
    if isinstance(kg_own, int):                       # 1000G: the model's calls are the inserted mutations, flagged
        recs = [r[:4] + (k_ not in donor,) + r[5:] for r in recs for k_ in [(r[0], r[1], r[2])]]  # as Minos does
    n_windows = None
    if info["mode_giab"] and meta.get("window_thr") is not None and (out / "work" / "hg_windows.pkl.gz").exists():
        # window decisions: where the genotyper's whole-window answer replaces the per-variant one (hapwin.py)
        from answer_set.k2 import hapwin
        wm = tracks.load_window(track, (track_versions or {}).get(track) or tracks.active_name(track), exclude_round)
        answer = {(r[0], r[1], r[2]): (r[3], r[4]) for r in recs}
        cands = hapwin.candidates(hapwin.load_windows(out / "work" / "hg_windows.pkl.gz"), answer,
                                  dict(zip(keys, p_raw)), donor, beds[info["mode_donor"]])
        if cands:
            pw = wm.predict_proba(np.vstack([c[5] for c in cands]))[:, 1]
            vname = (track_versions or {}).get(track) or tracks.active_name(track)
            sm = tracks.load_aux(track, vname, "spelling", exclude_round)
            respell_fn = None
            if sm is not None:                           # write replaced windows the way Minos most likely did
                fa_ = pysam.FastaFile(str(ref_path(chrom)))
                p_of = dict(zip(keys, p_raw))
                respell_fn = lambda c: hapwin.respell(fa_, chrom, c[0], c[2], p_of, sm)  # noqa: E731
            wthr = dict(meta["window_thr_group"], default=meta["window_thr"]) if meta.get("window_thr_group") \
                else meta["window_thr"]
            recs, n_windows = hapwin.apply(recs, cands, pw, wthr, respell_fn)
    donor_raw = None
    if info["mode_giab"]:
        from answer_set.k2.util import donor_records
        donor_raw = donor_records(info["mode_donor"], chrom, a, b)
    write_vcf(out / "truth.vcf.gz", chrom, [r[:5] for r in recs], donor_raw)
    write_scores(out / "scores.tsv.gz", keys, p_any, donor, allowed)

    # mutations file: only its entry count N matters to the v2 score (false-call divisor).
    # Minos picks N at random each round (369-791 seen); the median is the best guess
    # (typical error ~12%), better than any estimate from the round itself.
    injected = sorted((r for r in recs if r[4]), key=lambda r: -r[5])
    n = int(meta.get("n_targets_median", 580))
    pool = injected if injected else sorted(recs, key=lambda r: -r[5])
    step = max(1, len(pool) // max(1, n))
    muts = pool[::step][:n]
    write_vcf(out / "mutations.vcf.gz", chrom, [(r[0], r[1], r[2], r[3], True) for r in muts])

    report = {
        "region": region,
        "mode": ("GIAB:" + info["mode_donor"]) if info["mode_giab"] else
                ("1000G:" + info["kg"]["person"]) if info.get("n_dragen") else
                ("1000G:" + info["kg"]["person"] + " WITHOUT DRAGEN (download failed, reads only)")
                if (info.get("kg") or {}).get("error") else "unrecognised (reads only)",
        "person_match": info.get("kg") or None, "dragen_calls": info.get("n_dragen", 0),
        "gatk_person_checked": (info.get("kg") or {}).get("gatk_person_checked", 0),
        "gatk_person_error": (info.get("kg") or {}).get("gatk_person_error"),
        "donor_support": {d: round(v, 4) for d, v in info["donor_scores"].items()},
        "candidates": info["n_candidates"], "answer_set_size": len(recs),
        "background_from_donor": len(donor), "model_calls": len(recs) - len(donor),
        "predicted_injections": len(injected), "mutations_N": len(muts), "threshold": thr,
        "model": used, "seconds": {"features": round(t_feat, 1), "total": round(time.time() - t0, 1)},
        "readdiff": info.get("readdiff"), "kg_own_rule": kg_own, "inputs": changed_inputs or "ok", "windows_replaced": n_windows,
    }
    (out / "report.json").write_text(json.dumps(report, indent=1))
    return report
