# answer_set

Rebuilds the round's truth at round start, so configs can be scored without the validator.

Public API: `builder(name)` / `BUILDER_NAMES` (which builder main runs - `builders.py`), `K2Provider` (the builder
main runs, `k2/`), `OracleProvider` (released truth - backtests and upper bound), `AnswerSetProvider` (protocol:
`build(round)` once, `update(answer_set, callsets)` during the round).

## What the truth is

`truth.vcf` = one donor's **germline** + ~5K injected variants flagged `SYNTHETIC`. The score's core uses **every**
record, so both halves matter equally. `mutations.vcf` (the targets) only provides `N`.

## `K2Provider` (model-based, `answer_set/k2/`)

Ported from the k2/107 answer-set project (live build only; training and analysis stay in that project). Candidates
are every allele the round's reads show (bcftools mpileup, no calling); each gets ~110 measurements from public
evidence (GIAB NIST lists, Google originals of GIAB reads, 1000G panel, DRAGEN and GATK calls, region masks); a
gradient-boosted model per round format decides P(in truth) and P(hom). Registered in `builders.py` as
`builder("k2")`, the default builder main runs; `builder("oracle")` (released rounds only) replays a round guided by
the released truth.

| Module | Role |
|---|---|
| `k2/provider.py` | `K2Provider`: build into `<cache_dir>/answer_set/k2/<tag>/v1/`, read it into an `AnswerSet` |
| `k2/build.py` | candidates -> features -> track model -> `truth.vcf.gz`, `mutations.vcf.gz` (N=580), `scores.tsv.gz`, `report.json` |
| `k2/features.py` | read scan, format detection (GIAB / 1000G), every measurement |
| `k2/kg.py`, `k2/masks.py` | 1000G person, DRAGEN / GATK lookups, panel index; region masks |
| `k2/kg_origin.py` | 1000G original reads: the person's public 1000 Genomes CRAM (EBI), round region only |
| `k2/giab_origin.py`, `k2/readdiff.py`, `k2/localhap.py` | GIAB original reads, read-vs-original comparison, local haplotypes |
| `k2/hapgeno.py`, `k2/hapcand.py` | edit-aware local haplotype genotyper (GIAB): the person's two local sequences per spot |
| `k2/hapwin.py`, `k2/train_window.py`, `k2/train_spelling.py` | window decisions: where the genotyper's whole-window answer replaces the per-variant one, and how it is spelled |
| `k2/tracks.py`, `k2/assemble.py` | model versions per track (ACTIVE, rollback); scores -> records |
| `k2/util.py` | input and cache paths (`SETTINGS`), bcftools in Docker, VCF helpers |

Released models are code: `k2/models/<track>/ACTIVE` names the version main runs (GIAB g026_20261004T0742, 1000G k012_20261004T1750) and
`k2/models/<track>/<version>/` holds its `model.joblib`, `hom.joblib`, `meta.json` (~3 MB per track). Releasing a new
version = `tracks.promote` (copies those three files here) + commit; rollback = `--rollback` or git revert.
`python -m answer_set.k2.tracks` lists versions.
Inputs, under `SETTINGS.data_dir` (listed in `k2/util.py`): the shared `reference/`, `giab/<donor>.vcf.gz|.bed` and
`1kg/{phased,index,gatk_sites,masks}` (the last three: the builder's 1000G inputs). This clone's `cache/answer_set/k2/`: one build per
round, downloads (Google originals, per-person GATK), mask intervals, trained versions with their leave-out models,
prepared training rounds, the input manifest. Training: `k2/prepare.py` -> `k2/train_track.py` (GIAB: then `k2/train_window.py --base <version>`, then
`k2/train_spelling.py --version <that version>`).

- Needs Docker (the validator's bcftools 1.20 image), internet (S3 DRAGEN, EBI GATK, Google NovaSeq originals),
  `edlib` and `parasail` (requirements.txt). ~1-1.5 min per round on 16 cores.
- Released rounds (`round.released_dir` set) use the leave-out fold model that never saw the round.
- 1000G rounds: nothing is flagged `SYNTHETIC`, so every record is `germline` (the score uses all records).
- Never silently weaker, never stuck: the GIAB read comparison (Google originals) gets up to 3 tries within 900 s
  (a bug is not retried); if it still fails the build goes on without it (inserted 1-bp indels ~2 points worse).
  Every download and Docker call has a time limit, because main builds the answer set before it can submit.
- Inputs are checked against a manifest (`cache/answer_set/k2/inputs.json`: sha256, size, mtime;
  `python -m answer_set.k2.inputs record|verify [--full]`) - the same model gave different answer sets when the
  k2/107 machine replaced its HG002 VCF between two runs. A failed read comparison and changed inputs both land in
  `notes["warnings"]` of the AnswerSet and in the log.
- g007 (released 2026-10-02): g006's features, trained on 108 GIAB rounds (09-23..10-01 13:05). On 12 held-out
  rounds: F1 median 0.99905 (g006 0.99895, better in 12/12), v2 core +0.011 pts/round. Donor variants 100%; the
  remaining errors are mostly inserted 1-bp indels in repeats, which the reads show in another spelling.
- g009 (released 2026-10-03): adds the local haplotype genotyper (`k2/hapgeno.py`). Per spot (~50K windows/round) it
  reads every read's DNA between two fixed anchors (whatever the alignment), compares each read with its Google
  original (edited by Minos or not), scores every read against the candidate local sequences with its base qualities
  and picks the most likely pair; the chosen sequences are written as events the way Minos writes them (indels, often
  a deletion + an insertion a few bp apart). Its events become candidates and 11 features (posterior, P(hom), edited
  reads and distinct originals, unedited reads, edit base quality, local edit noise, read-comparison support) feed the
  model, which decides - as a direct override it broke far more than it fixed. Labels also accept other spellings of
  the same DNA within 30 bp. Same 108 training rounds as g007; on 28 rounds (22 unseen + 6 leave-out), hap.py vs the
  released truth: median score 99.627 -> 99.682, better in 24/28; missed inserted 1-bp indels 18.8 -> 12.1 per round,
  extras 6.5 -> 7.3. ~+20-40 s per GIAB build.
- g013 (released 2026-10-03), two changes on g009:
  - training counts points, not records (`train_track.py`): each row weighs (class weight / class truth count) **
    0.5 - one missed 1-bp indel costs the score ~20 SNP errors - and each class gets its own cut-off chosen on
    held-out folds for weight x F1 (base version g011);
  - window decisions (`hapwin.py`): among genotyped windows inside the confident region whose called pair differs
    from the per-variant answer, a second model (trained on the same rounds and folds: is the called pair the
    truth's DNA?) decides where the called sequences replace the answer's records. The person's GIAB records are
    never removed or re-genotyped. Held out on 108 rounds: 850 of 1,404 fixable windows fixed, 212 broken; ~9
    windows replaced per round.
  On 28 rounds (22 unseen + 6 leave-out), hap.py vs the released truth: median 99.682 -> 99.731, mean 99.667 ->
  99.732 (better than its base g011 in 27/28); inserted 1-bp indels missed 12.1 -> 8.9 per round, extras 7.3 -> 5.5;
  points lost per round 0.231 -> 0.179 (1-bp indels 0.154 -> 0.113; before the genotyper, g007: 0.296 / 0.199).
- g017 (released 2026-10-03) = g013 + the spelling of replaced windows. The validator's scorer classes each truth
  record by its spelling, so where our answer set stands in for the truth (local scoring of configs), a window
  written as one SNP where Minos wrote a deletion + an insertion lands in the wrong class (weight 0.18 instead of
  0.40). Options: the decompositions of 8 alignment schemes and every insertion + deletion pair (Minos often shifts
  a repeat that way, 10-20 bp apart); a model trained on the same rounds picks the one Minos most likely wrote.
  Options whose records overlap after left-alignment (an insertion inside a deletion) are dropped: hap.py cannot
  match them. Held out on 108 rounds, windows the window decision fixes: exactly Minos's records 59.9% -> ~72%
  (reachable 83.3%). 28 rounds: hap.py score unchanged (identical to g013 in 27/28, mean 99.732 -> 99.731);
  inserted 1-bp indel records missed 375 -> 354, wrong het/hom 148 -> 131, extra 127 -> 97. Fresh round 2026-10-03
  16:41 (HG003): 99.483 -> 99.508.
- g024 (released 2026-10-04) = g017's models with per-class cut-offs chosen by the validator's FULL score (core +
  false-positive part), not weight x F1: an extra call costs its class's F1 and ~0.0065 points of the false-positive
  part, a missed one only the F1 - so extra SNPs cost ~14x a missed one. SNP 0.80, 1 bp 0.60, 2-3 bp 0.75, 4-7 bp
  0.60, >= 8 bp 0.30. 28 rounds: +0.016 per round (better 20/28; SNP het missed 7.9 -> 11.0, extra calls in the
  region 13.5 -> 11.3); 6 fresh rounds (2026-10-03 16:41-22:41, unused in any choice): +0.025 (5/6).
- g026 (released 2026-10-04) = g024 + window cut-offs per group of windows, chosen by validator points on held-out
  training rounds (a replaced window gains or loses its records' class costs and the false-positive part's cost of
  extra calls): 1-3 bp indel windows 0.55, SNP-only 0.85 (at 0.40 they lost points), longer indels 0.75. 28 rounds
  +0.006 per round (better in 15, worse in 8); fresh rounds not used for tuning: 2026-10-03 +0.029 (6/6), 2026-10-04
  +0.016 (4/5); ~3 fewer windows replaced and 1-3 fewer extra calls per round.
- Writing the person's records: exactly as Minos's truth writes them - NIST's REF, ALT order and genotype order
  (2|1 -> 2/1, `util.donor_records`). The validator classes a two-allele record by the allele its genotype names
  first; writing 1/2 for NIST's 2/1 moved ~20 records per round out of the >= 8 bp class (found by the search
  agent's top-band ranking test on 2026-10-02 11:53: local scores under-rated configs with better 4-7 / >= 8 bp
  calling). 28 rounds, truth records per class vs the released truth: >= 8 bp shift 19.8 -> 0.2, 4-7 bp 6.0 -> 1.4,
  2-3 bp 5.1 -> 2.2; our own score unchanged.

- k007 (1000G, released 2026-10-04) = k006's rounds, folds and measurements, with per-class cut-offs chosen on
  held-out folds by weight x F1 minus each extra record's false-positive cost (~3.75 / N points; 1000G answers
  carried ~110 extras per round): SNP 0.90, 1 bp 0.65, 2-3 bp 0.65, 4-7 bp 0.75, >= 8 bp 0.55. Point weighting
  (alpha 0.5 / 1) did not help on 1000G. Scored like the validator on all 43 past 1000G rounds: 98.018 vs k006
  97.896, +0.122 per round (better 35/43); on the 16 rounds neither model trained on +0.120 (14/16). The gain is the
  false-positive part (extras ~110 -> ~83 per round); indel F1 slightly lower.
- 1000G truth, measured on the released rounds: the person's public calls that BOTH DRAGEN and GATK made (DRAGEN-only
  calls are never kept), minus a soft, quality-dependent share Minos drops (indels kept 89%, SNPs 97%; by GATK
  VQSLOD), plus ~6,700 inserted mutations per round with no `SYNTHETIC` flag. The reads are the person's real public
  1000G reads (2x150, ~34x): ~82% of pairs kept, each read cut to its first 128 bases, mutations and some noise
  added, names hashed, realigned.
- k010 (1000G, released 2026-10-04) = k007 + the read comparison: each round read against the person's public
  original read (`k2/kg_origin.py`, ~4 s for a region; `readdiff.py` unchanged, ~6 s), as for GIAB - it sees every
  inserted mutation (one round checked: 641/641 inserted 1-bp indels, 6,653 of its 6,740 strong changes in the
  truth). Its measurements (`tracks.RD_FEATS`) and read-comparison-only candidates are on only for 1000G models
  trained with them (meta `kg_readdiff`; `build.py` sets `features.KG_READDIFF`). Same 27 rounds and folds as k007.
  All 43 past 1000G rounds: 98.223 vs k007 98.018, +0.205 per round (better 42/43, worst -0.015); the 16 rounds
  neither trained on +0.189 (16/16); a new round (2026-10-04 10:41) +0.141. Gains: 1-bp indels +0.065 (extras 30 ->
  24), the false-positive part +0.067, longer indels +0.064. If the originals cannot be had, the build goes on
  without them (NaN measurements), as for GIAB.
- k011 (1000G) = k010 trained on all 43 past 1000G rounds: +0.064 per round over k010 (better 33/43).
- **1000G truth, the person's part, exactly** (found 2026-10-04): the person's calls in the raw 1000G high-coverage
  GATK callset whose site FILTER is PASS, with the person's DP >= 10 and GQ >= 20, at GATK's genotype - DRAGEN plays
  no part. 590,556 calls in 43 rounds: every call passing is in the truth, none failing is. The earlier "soft,
  quality-dependent drop" was these per-person fields, which no measurement carried.
- k012 (1000G, released 2026-10-04) = k011's model + that rule (`kg.gatk_person_own`, meta `kg_own_rule`): the
  person's part is copied like GIAB's NIST list and the person's other GATK calls are never written; the model decides
  only the rest (the inserted mutations). If the person's GATK calls cannot be had, the model decides everything, as
  before. All 44 1000G rounds: 99.924 vs k011 98.297, +1.627 per round (better 44/44, worst round 99.700); the new
  round 2026-10-04 10:41 99.865. Left: ~4 missed + ~2 extra inserted 1-bp indels per round (~0.08 points).

## Open work

- **Inserted 1-bp indels** still cause most of the remaining loss (~9 missed + 5.5 extra per round with g013): the
  ~40% of fixable windows the window model leaves alone, sequences noise splits so the genotyper never lists them,
  spots with no window.
- **Spelling**: ~28% of replaced windows are still not spelled exactly like Minos (~17% need spellings no option
  produces yet).
- **1000G**: the person's part is exact (k012); what is left is the inserted mutations (~4 missed + ~2 extra 1-bp
  indels per round) - the model still trains on every candidate; training it only on candidates outside the person's
  list (as GIAB) may help. A new round (2026-10-04 10:41) had fewer extras (62) and more inserted
  1-bp indels (19% of inserted, before ~11%) - watch whether Minos changed its 1000G mix.
- **`update()` during a live round** (parked by the user).
