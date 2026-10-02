"""1000 Genomes (3,202-person high-coverage panel) support: identify a round's person, extract their variants.

Minos's "synthetic" rounds turned out to be real 1000 Genomes individuals: their variants
match one panel person 100% with identical genotypes. This module finds that person from
the round's reads and returns their variants, the same way GIAB rounds use GIAB lists.

  python -m answer_set.k2.kg --index chr18      # one-time: build the lookup index
  python -m answer_set.k2.kg --index all
"""
import argparse
import subprocess
import time
from pathlib import Path

import numpy as np

from answer_set.k2.util import KG_DIR as KG_ROOT, BCFTOOLS_IMG, CACHE, bcftools, read_alleles, ref_path, trim
from common.settings import SETTINGS

KG_DIR = Path(SETTINGS.kg_dir) / "phased"      # the phased panel (scripts/fetch_reference_data.sh)
INDEX_DIR = KG_ROOT / "index"
CHROMS = [f"chr{i}" for i in range(14, 23)]
THIN_BP = 1000          # at most one index site per 1 kb
MIN_AF, MAX_AF = 0.05, 0.95


def panel(chrom):
    return KG_DIR / f"{chrom}.vcf.gz"


def index_path(chrom):
    return INDEX_DIR / f"index_{chrom}.npz"


def build_index(chrom):
    """Common biallelic SNPs (5-95% frequency), thinned to 1 per kb, with every person's genotype."""
    INDEX_DIR.mkdir(parents=True, exist_ok=True)
    cmd = (f"bcftools query -l /k/{panel(chrom).name} > /x/.samples_{chrom}; "
           f"bcftools view -m2 -M2 -v snps -i 'INFO/AF>{MIN_AF} && INFO/AF<{MAX_AF}' /k/{panel(chrom).name} -Ou "
           f"| bcftools query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n'")
    args = ["docker", "run", "--rm", "-v", f"{KG_DIR}:/k", "-v", f"{INDEX_DIR}:/x", BCFTOOLS_IMG, "sh", "-c", cmd]
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    pos, ref, alt, rows = [], [], [], []
    last = -THIN_BP
    code = {"0|0": 0, "0|1": 1, "1|0": 1, "1|1": 2, "0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2}
    for line in proc.stdout:
        p = line.split("\t", 3)
        x = int(p[0])
        if x - last < THIN_BP:
            continue
        last = x
        gts = line.rstrip("\n").split("\t")[3:]
        pos.append(x); ref.append(p[1]); alt.append(p[2])
        rows.append(np.fromiter((code.get(g, 0) for g in gts), dtype=np.int8, count=len(gts)))
    proc.wait()
    samples = (INDEX_DIR / f".samples_{chrom}").read_text().split()
    np.savez_compressed(index_path(chrom), pos=np.array(pos), ref=np.array(ref), alt=np.array(alt),
                        G=np.vstack(rows), samples=np.array(samples))
    return len(pos), len(samples)


def identify(chrom, a, b, cand_frac):
    """Which panel person is this round?

    cand_frac: {(pos, ref, alt): share of reads showing the allele} for read candidates.
    Returns (sample, score, second_score). Score = genotype agreement at index sites where
    either the reads or the person show a variant.
    """
    z = np.load(index_path(chrom))
    m = (z["pos"] >= a) & (z["pos"] <= b)
    pos, ref, alt, G = z["pos"][m], z["ref"][m], z["alt"][m], z["G"][m]
    obs = np.zeros(len(pos), dtype=np.int8)
    for i, (p, r, x) in enumerate(zip(pos, ref, alt)):
        f = cand_frac.get((int(p), r, x), 0.0)
        obs[i] = 0 if f < 0.15 else 2 if f > 0.8 else 1
    informative = (G > 0) | (obs[:, None] > 0)
    agree = ((G == obs[:, None]) & informative).sum(0) / np.maximum(1, informative.sum(0))
    order = np.argsort(-agree)
    return str(z["samples"][order[0]]), float(agree[order[0]]), float(agree[order[1]]), int(len(pos))


def sample_alleles(chrom, a, b, sample, workdir):
    """{key: 'het'|'hom'} for one panel person in the region (normalised like everything else)."""
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / f"kg_{sample}.vcf.gz"
    if not out.exists():
        bcftools(f"bcftools view -s {sample} -r {chrom}:{a}-{b} --min-ac 1 -Ou /k/{panel(chrom).name} "
                 f"| bcftools view -e 'GT=\"0|0\" || GT=\"0/0\" || ALT~\"<\"' -Ou "
                 f"| bcftools norm -m -both -f /r/{chrom}.fa -Oz -o /o/{out.name} && bcftools index -f -t /o/{out.name}",
                 [(str(KG_DIR), "/k"), (str(ref_path(chrom).parent), "/r"), (str(workdir), "/o")])
    return {k: z for k, z, _ in read_alleles(out)}


DRAGEN_URL = "https://1000genomes-dragen.s3.amazonaws.com/data/dragen-3.7.6/hg38-graph-based"


GATK_DIR = KG_ROOT / "gatk_sites"


def gatk_sites(chrom, a, b, workdir):
    """{key: (pass_filter, allele_count)} from the raw 1000 Genomes GATK callset (sites only,
    local copy). Minos's key for a 1000 Genomes person keeps only calls GATK also made, so a
    DRAGEN call at a site GATK never reported is (almost) never in the key."""
    import pysam
    src = GATK_DIR / f"{chrom}.sites.vcf.gz"
    if not src.exists():
        return None
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / "gatk_sites.vcf.gz"
    if not out.exists() or _count_records(out) == 0:
        bcftools(f"set -o pipefail; bcftools view -G -r {chrom}:{a}-{b} /g/{src.name} -Ou 2>/dev/null "
                 f"| bcftools norm -m -both -f /r/{chrom}.fa -Oz -o /o/{out.name} 2>/dev/null "
                 f"&& bcftools index -f -t /o/{out.name}",
                 [(str(GATK_DIR), "/g"), (str(workdir), "/o"), (str(ref_path(chrom).parent), "/r")])
    res = {}
    for rec in pysam.VariantFile(str(out)):
        ac = rec.info.get("AC")
        ac = ac[0] if isinstance(ac, tuple) else (ac or 0)
        passed = list(rec.filter.keys()) in (["PASS"], [])
        ann = {}
        for f in SITE_INFO:
            try:
                v = rec.info.get(f)
            except (ValueError, KeyError):
                v = None
            v = v[0] if isinstance(v, tuple) else v
            ann[f] = float(v) if v is not None else float("nan")
        for alt in rec.alts or ():
            if alt.startswith("<") or alt == "*":
                continue
            res[trim(rec.pos, rec.ref, alt)] = (passed, int(ac or 0), ann)
    return res


# population-level GATK annotations: Minos drops more DRAGEN+GATK calls where these look poor
# (low VQSLOD, many missing genotypes (AN), Hardy-Weinberg / excess-het failures, low MQ)
SITE_INFO = ("VQSLOD", "AN", "HWE", "ExcHet", "MQ", "MQRankSum", "FS", "SOR", "DP", "InbreedingCoeff")


def panel_sites(chrom, a, b, workdir):
    """Alleles of the official filtered 1000 Genomes panel (any person) in the region. Minos keeps
    99.8% of DRAGEN+GATK SNPs at panel sites but only ~63% elsewhere."""
    import pysam
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / "panel_sites.vcf.gz"
    if not out.exists() or _count_records(out) == 0:
        bcftools(f"set -o pipefail; bcftools view -G -r {chrom}:{a}-{b} /k/{panel(chrom).name} -Ou "
                 f"| bcftools norm -m -both -f /r/{chrom}.fa -Oz -o /o/{out.name} && bcftools index -f -t /o/{out.name}",
                 [(str(panel(chrom).parent), "/k"), (str(workdir), "/o"), (str(ref_path(chrom).parent), "/r")])
    return {trim(rec.pos, rec.ref, alt) for rec in pysam.VariantFile(str(out)) for alt in rec.alts or ()}


GATK_RAW = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
            "20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrom}.recalibrated_variants.vcf.gz")
GATK_PERSON_CACHE = CACHE / "kg_gatk_person"


def gatk_person_at(chrom, a, b, sample, loci, chunks=6, attempts=10):
    """What the raw 1000 Genomes GATK callset says for ONE person, fetched only around `loci`
    (the person's DRAGEN indels of >=4 bp). Minos's key for a 1000 Genomes person keeps only calls
    both DRAGEN and GATK made, so this is the missing piece for indels.

    Returns {key: 'het'|'hom'} of GATK's calls near the loci (cached per round region)."""
    import concurrent.futures as cf
    import pysam
    GATK_PERSON_CACHE.mkdir(parents=True, exist_ok=True)
    final = GATK_PERSON_CACHE / f"{sample}_{chrom}_{a}_{b}.vcf.gz"
    if not final.exists():
        wins = []
        for k in sorted(loci):
            s0, e0 = k[0] - 5, k[0] + len(k[1]) + 5
            if wins and s0 - wins[-1][1] < 2000:
                wins[-1][1] = max(wins[-1][1], e0)
            else:
                wins.append([s0, e0])
        if not wins:
            return {}
        tmp = GATK_PERSON_CACHE / f".tmp_{sample}_{chrom}_{a}_{b}"
        tmp.mkdir(exist_ok=True)
        parts = [wins[i::chunks] for i in range(chunks)]

        def fetch(i):
            if not parts[i]:
                return None
            reg = tmp / f"r{i}.tsv"
            reg.write_text("".join(f"{chrom}\t{x}\t{y}\n" for x, y in sorted(parts[i])))
            out = tmp / f"p{i}.vcf.gz"
            for t in range(attempts):
                try:
                    bcftools(f"set -o pipefail; cd /o && bcftools view -s {sample} -R /o/{reg.name} "
                             f"{GATK_RAW.format(chrom=chrom)} -Ou "
                             f"| bcftools view -e 'GT=\"0/0\" || GT=\"./.\" || ALT~\"<\" || ALT=\"*\"' -Ou "
                             f"| bcftools norm -m -both -f /r/{chrom}.fa -Oz -o /o/{out.name} "
                             f"&& bcftools index -f -t /o/{out.name}",
                             [(str(tmp), "/o"), (str(ref_path(chrom).parent), "/r")])
                    return out
                except RuntimeError as e:
                    err = e
                    out.unlink(missing_ok=True)
                    time.sleep(5 * (t + 1))
            raise RuntimeError(f"GATK lookup for {sample} chunk {i} failed after {attempts} attempts: {err}")

        with cf.ThreadPoolExecutor(chunks) as ex:
            outs = [o for o in ex.map(fetch, range(chunks)) if o]
        names = " ".join(f"/o/{o.name}" for o in outs)
        bcftools(f"set -o pipefail; cd /o && bcftools concat -a {names} -Oz -o /o/merged.vcf.gz "
                 f"&& bcftools index -f -t /o/merged.vcf.gz", [(str(tmp), "/o")])
        (tmp / "merged.vcf.gz").rename(final)
        (tmp / "merged.vcf.gz.tbi").rename(str(final) + ".tbi")
        for f_ in tmp.iterdir():
            f_.unlink()
        tmp.rmdir()
    res = {}
    for rec in pysam.VariantFile(str(final)):
        gt = rec.samples[0]["GT"]
        if not gt or any(g is None for g in gt):
            continue
        for i_ in {g for g in gt if g}:
            alt = rec.alleles[i_]
            if not alt.startswith("<") and alt != "*":
                res[trim(rec.pos, rec.ref, alt)] = "hom" if len(set(gt)) == 1 else "het"
    return res


def _person_path(chrom, a, b, sample):
    # _full2: normalised with --old-rec-tag, so the raw record layout is kept (gatk_person_layout)
    return GATK_PERSON_CACHE / f"{sample}_{chrom}_{a}_{b}_full2.vcf.gz"


def gatk_person_layout(chrom, a, b, sample):
    """{key: (alleles in the raw GATK record, which of them is this allele, other raw records
    overlapping it)}. Minos drops more calls from multi-allele / overlapping GATK records."""
    import bisect
    import pysam
    gatk_person_region(chrom, a, b, sample)
    raw, per_key = {}, {}
    for rec in pysam.VariantFile(str(_person_path(chrom, a, b, sample))):
        gt = rec.samples[0]["GT"]
        if not gt or any(g is None for g in gt):
            continue
        o = rec.info.get("OLD")
        o = (",".join(o) if isinstance(o, tuple) else o) if o else f"{rec.chrom}|{rec.pos}|{rec.ref}|{','.join(rec.alts)}"
        p = o.split("|")
        start, ref, alts = int(p[1]), p[2], p[3].split(",")
        raw[(start, ref, p[3])] = start + len(ref) - 1
        for i_ in {g for g in gt if g}:
            alt = rec.alleles[i_]
            if not alt.startswith("<") and alt != "*":
                idx = int(p[4]) if len(p) > 4 else 1
                per_key[trim(rec.pos, rec.ref, alt)] = (len(alts), idx, (start, start + len(ref) - 1))
    spans = sorted((s_, e_) for (s_, _, _), e_ in raw.items())
    starts = [x[0] for x in spans]
    out = {}
    for k, (n, idx, (s_, e_)) in per_key.items():
        j, ov = bisect.bisect_left(starts, s_ - 300), 0
        while j < len(spans) and spans[j][0] <= e_:
            if spans[j][1] >= s_ and spans[j] != (s_, e_):
                ov += 1
            j += 1
        out[k] = (n, idx, ov)
    return out


def gatk_person_region(chrom, a, b, sample, chunks=8, attempts=10):
    """The person's GATK calls for the WHOLE region (all SNPs and indels), fetched as `chunks`
    parallel pieces of the indexed remote callset; cached. Returns {key: 'het'|'hom'}."""
    import concurrent.futures as cf
    import pysam
    GATK_PERSON_CACHE.mkdir(parents=True, exist_ok=True)
    final = _person_path(chrom, a, b, sample)
    if not final.exists():
        tmp = GATK_PERSON_CACHE / f".tmpfull2_{sample}_{chrom}_{a}_{b}"
        tmp.mkdir(exist_ok=True)
        step = (b - a + chunks) // chunks
        pieces = [(a + i * step, min(b, a + (i + 1) * step - 1)) for i in range(chunks) if a + i * step <= b]

        def fetch(i):
            x, y = pieces[i]
            out = tmp / f"p{i}.vcf.gz"
            for t in range(attempts):
                try:
                    bcftools(f"set -o pipefail; cd /o && bcftools view -s {sample} -r {chrom}:{x}-{y} "
                             f"{GATK_RAW.format(chrom=chrom)} -Ou "
                             f"| bcftools view -e 'GT=\"0/0\" || GT=\"./.\" || ALT~\"<\" || ALT=\"*\"' -Ou "
                             f"| bcftools norm -m -both --old-rec-tag OLD -f /r/{chrom}.fa -Oz -o /o/{out.name} "
                             f"&& bcftools index -f -t /o/{out.name}",
                             [(str(tmp), "/o"), (str(ref_path(chrom).parent), "/r")])
                    if _count_records(out) > 0:
                        return out
                except RuntimeError as e:
                    err = e
                out.unlink(missing_ok=True)
                time.sleep(5 * (t + 1))
            raise RuntimeError(f"GATK region fetch for {sample} {chrom}:{x}-{y} failed after {attempts} attempts")

        with cf.ThreadPoolExecutor(len(pieces)) as ex:
            outs = list(ex.map(fetch, range(len(pieces))))
        names = " ".join(f"/o/{o.name}" for o in outs)
        bcftools(f"set -o pipefail; cd /o && bcftools concat -a {names} -Oz -o /o/merged.vcf.gz "
                 f"&& bcftools index -f -t /o/merged.vcf.gz", [(str(tmp), "/o")])
        (tmp / "merged.vcf.gz").rename(final)
        (tmp / "merged.vcf.gz.tbi").rename(str(final) + ".tbi")
        for f_ in tmp.iterdir():
            f_.unlink()
        tmp.rmdir()
    res = {}
    for rec in pysam.VariantFile(str(final)):
        gt = rec.samples[0]["GT"]
        if not gt or any(g is None for g in gt):
            continue
        for i_ in {g for g in gt if g}:
            alt = rec.alleles[i_]
            if not alt.startswith("<") and alt != "*":
                res[trim(rec.pos, rec.ref, alt)] = "hom" if len(set(gt)) == 1 else "het"
    return res


def _count_records(path):
    import gzip
    try:
        return sum(1 for l in gzip.open(path, "rt") if not l.startswith("#"))
    except (OSError, EOFError):
        return 0


def dragen_alleles(chrom, a, b, sample, workdir):
    """{key: (zyg, qual, gq)} from the person's public DRAGEN 3.7.6 calls (PASS only)."""
    import pysam
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / f"dragen_{sample}.vcf.gz"
    for attempt in range(3):
        if out.exists() and _count_records(out) > 0:
            break
        out.unlink(missing_ok=True)
        try:
            # pipefail: a failed/partial download from S3 must fail, not produce an empty file
            bcftools(f"set -o pipefail; cd /o && bcftools view -f PASS -r {chrom}:{a}-{b} "
                     f"{DRAGEN_URL}/{sample}/{sample}.hard-filtered.vcf.gz -Ou "
                     f"| bcftools view -e 'GT=\"0/0\" || GT=\"./.\" || ALT~\"<\"' -Ou "
                     f"| bcftools norm -m -both -f /r/{chrom}.fa -Oz -o /o/{out.name} && bcftools index -f -t /o/{out.name}",
                     [(str(workdir), "/o"), (str(ref_path(chrom).parent), "/r")])
        except RuntimeError:
            pass
        time.sleep(3 * (attempt + 1))
    if not out.exists() or _count_records(out) == 0:
        out.unlink(missing_ok=True)
        raise RuntimeError(f"could not fetch DRAGEN calls for {sample} {chrom}:{a}-{b} after 3 attempts")
    res = {}
    for rec in pysam.VariantFile(str(out)):
        s = rec.samples[0]
        gt = s.get("GT")
        if not gt or any(g is None for g in gt):
            continue
        zyg = "hom" if len(set(gt)) == 1 else "het"
        for i in {g for g in gt if g}:
            alt = rec.alleles[i]
            if alt.startswith("<") or alt == "*":
                continue
            res[trim(rec.pos, rec.ref, alt)] = (zyg, float(rec.qual or 0), float(s.get("GQ") or 0))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", required=True, help="chromosome (e.g. chr18) or 'all'")
    args = ap.parse_args()
    for c in CHROMS if args.index == "all" else [args.index]:
        if not panel(c).exists():
            print(f"{c}: panel not downloaded, skipped")
            continue
        n, s = build_index(c)
        print(f"{c}: {n} index sites x {s} people", flush=True)


if __name__ == "__main__":
    main()
