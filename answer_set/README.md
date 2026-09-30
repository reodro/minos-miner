# answer_set

Rebuilds the round's truth at round start, so configs can be scored without the validator.

Public API: `Provider` (the real builder), `OracleProvider` (released truth - backtests and upper bound),
`AnswerSetProvider` (protocol: `build(round)` once, `update(answer_set, callsets)` during the round).

## What the truth is

`truth.vcf` = one donor's **germline** + ~5K injected variants flagged `SYNTHETIC`. The score's core uses **every**
record, so both halves matter equally. `mutations.vcf` (the targets) only provides `N`.

## How `Provider.build` works

| Step | Code | Method | Measured |
|---|---|---|---|
| Donor | `fingerprint.py`, `giab.py` | 1000G: genotype ~1.5K panel SNPs from the BAM, best-matching sample. GIAB: genotype SNPs where the 7 NIST sets disagree | 1000G: right donor >= 0.991 vs <= 0.776; GIAB 24/24, 0.99+ vs <= 0.72 |
| Germline | `germline.py`, `giab.py` | 1000G: the sample's NYGC raw calls + `KG_RULE` (PASS, no `*`, GQ >= 20, DP >= 10). GIAB: the donor's NIST records as written | exact: 1000G 36/36 rounds, GIAB 24/24 |
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

Fallbacks (never raise): no parent -> synthetic at p=0.5 (left out of the truth); no donor -> `Regime.UNKNOWN`,
and main submits a static baseline without searching.

**Known risk - NA12878**: it is GIAB's HG001 and is also in 1000G. Every NA12878 round seen so far used the NIST truth,
so it is sent the GIAB way. That is a rule from past rounds, not something checked live: if the platform ever uses
NA12878's 1000G calls instead, we would build the wrong germline for that donor without noticing.

## KPI

`test/answer_set_kpi.py`: score the same callsets under our answer set and under the released truth.

| Round | Type | Spearman | Oracle gap (points) |
|---|---|---|---|
| 092900Z | 1000G | 0.984 | 0 |
| 055300Z | 1000G | 0.980 | 0 |
| 081700Z | 1000G | **0.841** | **0.211** |
| 152900Z | GIAB HG002 | 0.998 | 0 (17 similar callsets) |
| 224100Z | GIAB HG005 | 0.995 | 0 (17 similar callsets) |

GIAB rounds score ~1.2 points lower across the board under ours (constant offset: rankings unaffected).
These KPIs used the true N; rerun them with the estimate.

## Open work

- **081700Z gap**: find the truth records that differ and which difficulty classes they hit (small classes like
  indel >= 8 bp or SNP hom can move with a handful of records).
- **GIAB indels**: lower recall/precision than 1000G (different sequencing in the parent).
- **p=0.5 candidates** (thin or ambiguous parent, ~80-150 per round) are left out of the truth; decide from evidence.
- **N**: measure how much the 13% estimate error moves rankings; try better predictors.
- **Speed**: GIAB `classify` (pileups at 300x) takes ~60-80 s; `build()` must stay within a few minutes.
- **Setup**: the NIST VCFs/BEDs and reference FASTAs/SDFs are not downloaded by a script yet.

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
