# answer_set

Rebuilds the round's truth at round start, so configs can be scored without the validator.

Public API: `builder(name)` / `BUILDER_NAMES` (which builder main runs - `builders.py`), `Provider` (the real builder), `K2Provider` (model-based, `k2/`), `OracleProvider` (released truth - backtests and upper bound),
`AnswerSetProvider` (protocol: `build(round)` once, `update(answer_set, callsets)` during the round).

## What the truth is

`truth.vcf` = one donor's **germline** + ~5K injected variants flagged `SYNTHETIC`. The score's core uses **every**
record, so both halves matter equally. `mutations.vcf` (the targets) only provides `N`.

## How `Provider.build` works

| Step | Code | Method | Measured |
|---|---|---|---|
| Donor | `fingerprint.py`, `giab.py` | 1000G: genotype ~1.5K panel SNPs from the BAM, best-matching sample. GIAB: genotype SNPs where the 7 NIST sets disagree | 1000G: right donor >= 0.991 vs <= 0.776; GIAB 24/24, 0.99+ vs <= 0.72 |
| Germline | `germline.py`, `giab.py` | 1000G: the sample's NYGC raw calls + `KG_RULE` (PASS, no `*`, GQ >= 20, DP >= 10). GIAB: the donor's NIST records. Both **as whole records** (multi-allelic records keep every ALT, like the platform's truth) | record-identical to the released truth: 1000G 36/36 rounds, GIAB 24/24 |
| Parent reads | `parent.py` | the donor's own sequencing over the region: 1000G 30x CRAM (S3), GIAB NIST 100-300x BAM (NCBI), parallel chunks + retries | 1000G ~10 s, GIAB ~30-60 s |
| Synthetic pool | `synthetic.py` | bcftools candidates on the round BAM, germline alleles dropped, kept if the **parent shows reference** | 1000G: SNP ~0.998/0.997, indel ~0.99/0.985 recall/precision; GIAB HG005: SNP 0.995/0.985, indel 0.96/0.95 |
| N | `provider.py` | `0.1069 x` confident synthetic count (unknown live) | median error 13% |
| Output | `provider.py` | `truth.vcf.gz` + `mutations.vcf.gz` in the cache, named in `notes` for the scorers | - |

Why the parent works: round reads are **real sequencing** of the donor (downsampled, renamed). Injected variants
are in the round reads but in no real sequencing of the donor; real variants the truth leaves out are in both.
On 2 Mb: 0/1112 synthetic SNPs show the alt in the parent, 3189/3192 germline SNPs do.

Regime (GIAB vs 1000G) is decided by genotype, in this order: 1000G fingerprint >= 0.95 (and not NA12878) ->
1000G; else GIAB match >= 0.95 -> GIAB; else UNKNOWN. On a 1000G round every GIAB donor scores ~0.52; on a GIAB round
every 1000G sample scores <= ~0.69 (HG002-HG007 are not in 1000G). All 60 stored rounds got the right regime.

Degraded paths: no donor -> `Regime.UNKNOWN` (nothing else is built; main submits a static baseline without
searching); no parent reads -> every candidate goes into the truth unfiltered (real donor variants included -
degraded but scorable). `build()` can still raise on missing data; main then falls back to the static baseline.

**Known risk - NA12878**: it is GIAB's HG001 and is also in 1000G. Every NA12878 round seen so far used the NIST truth,
so it is sent the GIAB way. That is a rule from past rounds, not something checked live: if the platform ever uses
NA12878's 1000G calls instead, we would build the wrong germline for that donor without noticing.

## `K2Provider` (default builder: model-based, `answer_set/k2/`)

Ported from the k2/107 answer-set project (live build only; training and analysis stay in that project). Candidates
are every allele the round's reads show (bcftools mpileup, no calling); each gets ~110 measurements from public
evidence (GIAB NIST lists, Google originals of GIAB reads, 1000G panel, DRAGEN and GATK calls, region masks); a
gradient-boosted model per round format decides P(in truth) and P(hom). Same `AnswerSetProvider` interface as
`Provider`; registered in `builders.py` as `builder("k2")`, the default builder main runs (`Provider` stays available as `builder("provider")`).

| Module | Role |
|---|---|
| `k2/provider.py` | `K2Provider`: build into `<cache_dir>/answer_set/k2/<tag>/v1/`, read it into an `AnswerSet` |
| `k2/build.py` | candidates -> features -> track model -> `truth.vcf.gz`, `mutations.vcf.gz` (N=580), `scores.tsv.gz`, `report.json` |
| `k2/features.py` | read scan, format detection (GIAB / 1000G), every measurement |
| `k2/kg.py`, `k2/masks.py` | 1000G person, DRAGEN / GATK lookups, panel index; region masks |
| `k2/giab_origin.py`, `k2/readdiff.py`, `k2/localhap.py` | GIAB original reads, read-vs-original comparison, local haplotypes |
| `k2/tracks.py`, `k2/assemble.py` | model versions per track (ACTIVE, rollback); scores -> records |
| `k2/util.py` | input and cache paths (`SETTINGS`), bcftools in Docker, VCF helpers |

Released models are code: `k2/models/<track>/ACTIVE` names the version main runs (GIAB g007_20261002T1020, 1000G k006) and
`k2/models/<track>/<version>/` holds its `model.joblib`, `hom.joblib`, `meta.json` (~3 MB per track). Releasing a new
version = `tracks.promote` (copies those three files here) + commit; rollback = `--rollback` or git revert.
`python -m answer_set.k2.tracks` lists versions.
Inputs, under `SETTINGS.data_dir` (listed in `k2/util.py`): the shared `reference/`, `giab/<donor>.vcf.gz|.bed` and
`1kg/{phased,index,gatk_sites,masks}` (the last three: the builder's 1000G inputs). This clone's `cache/answer_set/k2/`: one build per
round, downloads (Google originals, per-person GATK), mask intervals, trained versions with their leave-out models,
prepared training rounds, the input manifest. Training: `k2/prepare.py` -> `k2/train_track.py`.

- Needs Docker (the validator's bcftools 1.20 image), internet (S3 DRAGEN, EBI GATK, Google NovaSeq originals),
  `edlib` and `parasail` (requirements.txt). ~1-1.5 min per round on 16 cores.
- Released rounds (`round.released_dir` set) use the leave-out fold model that never saw the round.
- 1000G rounds: nothing is flagged `SYNTHETIC`, so every record is `germline` (the score uses all records).
- Never silently weaker: the GIAB read comparison (Google originals) is tried 3 times; if it still fails the build
  goes on without it (inserted 1-bp indels ~2 points worse). Inputs are checked against a manifest
  (`cache/answer_set/k2/inputs.json`: sha256, size, mtime; `python -m answer_set.k2.inputs record|verify [--full]`)
  - the same model gave different answer sets when the k2/107 machine replaced its HG002 VCF between two runs.
  Both land in `notes["warnings"]` of the AnswerSet and in the log.
- g007 (released 2026-10-02): g006's features, trained on 108 GIAB rounds (09-23..10-01 13:05). On 12 held-out
  rounds: F1 median 0.99905 (g006 0.99895, better in 12/12), v2 core +0.011 pts/round. Donor variants 100%; the
  remaining errors are mostly inserted 1-bp indels absent from the reads (validated before downsampling).

## KPI

Score the same callsets under our answer set and under the released truth, then compare.

| Round | Type | Spearman | Oracle gap (points) |
|---|---|---|---|
| 092900Z | 1000G | 0.984 | 0 |
| 055300Z | 1000G | 0.980 | 0 |
| 081700Z | 1000G | **0.841** | **0.211** |
| 152900Z | GIAB HG002 | 0.998 | 0 (17 similar callsets) |
| 224100Z | GIAB HG005 | 0.995 | 0 (17 similar callsets) |

GIAB rounds score ~1.2 points lower across the board under ours (constant offset: rankings unaffected).
These KPIs were measured before the 1000G germline was written as whole records (it was split per allele, which
changes how the score classifies ~2000 multi-allelic records per round) and with the true N: rerun them.

## Open work

- **081700Z gap**: find the truth records that differ and which difficulty classes they hit (small classes like
  indel >= 8 bp or SNP hom can move with a handful of records).
- **GIAB indels**: lower recall/precision than 1000G (different sequencing in the parent).
- **p=0.5 candidates** (thin or ambiguous parent, ~80-150 per round) are left out of the truth; decide from evidence.
- **N**: measure how much the 13% estimate error moves rankings; try better predictors.
- **Speed**: GIAB `classify` (pileups at 300x) takes ~60-80 s; `build()` must stay within a few minutes.
- **Setup**: the NIST VCFs/BEDs are not downloaded by a script yet (reference genomes: `scripts/fetch_genome.sh`).

## Ideas (not planned)

- **Read-diff against parent BAMs** (legacy): only worked when round BAMs were edited copies of the parent BAM
  ("normal" rounds, none since 09-12). Build only if such rounds return.
- **Read-style detector**: read length + name format tell in seconds how a round's reads were made - old "normal"
  rounds (148/250 bp, original names: edited copies of the donor's public BAM) vs current ones (128 bp, hashed
  names). Only useful if old-style rounds return, for the read-diff trick above. (Not the GIAB/1000G regime:
  that is decided by genotype, see above.)
- **Re-test "the answer set cannot be purified"** (legacy claim that injected and background variants look identical
  in read statistics) before relying on it - the parent filter already sidesteps most of it.
- **`update()` during a live round** (parked by the user).
