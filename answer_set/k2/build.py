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


def write_vcf(path, chrom, records):
    """records: list of (pos, ref, alt, zyg, synthetic). Writes bgzipped + indexed VCF.

    Alleles sharing a position are merged into one record (genotype 1/2), expressed
    against the longest REF among them, as a person carries at most two alleles.
    """
    fai = {l.split()[0]: int(l.split()[1]) for l in open(str(ref_path(chrom)) + ".fai")}
    by_pos = {}
    for r in records:
        by_pos.setdefault(r[0], []).append(r)
    lines = []
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
        lines.append(f"{chrom}\t{pos}\t.\t{R}\t{','.join(alts)}\t60\tPASS\t"
                     f"{'SYNTHETIC' if syn else '.'}\tGT\t{gt}\n")
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
        features.RD_1BP_FACTOR = tracks.meta("giab", gname).get("rd_1bp_factor", 1.0)
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
    if len(keys) and meta.get("thr_only_rd") is not None:
        p_any = tracks.apply_source_thr(meta, X_full, p_any)
    if len(keys) and meta.get("gt_from_gatk"):
        p_hom = tracks.apply_gt_rule(meta, X_full, p_hom)

    donor = donors[info["mode_donor"]] if info["mode_giab"] else {}
    allowed = None
    if info["mode_giab"] and meta.get("giab_bed_rule"):
        bed = beds[info["mode_donor"]]
        allowed = lambda k: bed.contains(k[0]) and bed.contains(k[0] + len(k[1]) - 1)  # noqa: E731
    recs = assemble_records(keys, p_any, p_hom, thr, info["mode_giab"], donor, allowed)
    write_vcf(out / "truth.vcf.gz", chrom, [r[:5] for r in recs])
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
        "readdiff": info.get("readdiff"), "inputs": changed_inputs or "ok",
    }
    (out / "report.json").write_text(json.dumps(report, indent=1))
    return report
