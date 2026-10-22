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
| `k2/nest.py` | N when the platform gives none: from a raised pool of inserted mutations (>= 680/Mb: q x pool / 8, q ~0.97 from cleanly raised past rounds), else 1/E[1/N] of past same-type draws that fit this pool (N too high costs the pick more than too low). Walk-forward, 474 rounds: error in N 17.4% (median of the last 10 rounds) -> 12.7%; score lost by the pick (field submissions as candidates) 0.0104 -> 0.0069 per round |
| `k2/build.py` | candidates -> features -> track model -> `truth.vcf.gz`, `mutations.vcf.gz` (N is the setting `answer_set.n`, default 548), `scores.tsv.gz`, `report.json` |
| `k2/features.py` | read scan, format detection (GIAB / 1000G), every measurement |
| `k2/kg.py`, `k2/masks.py` | 1000G person, DRAGEN / GATK lookups, panel index; region masks |
| `k2/kg_origin.py` | 1000G original reads: the person's public 1000 Genomes CRAM (EBI), round region only |
| `k2/giab_origin.py`, `k2/readdiff.py`, `k2/localhap.py` | GIAB original reads, read-vs-original comparison, local haplotypes |
| `k2/hapgeno.py`, `k2/hapcand.py` | edit-aware local haplotype genotyper (GIAB): the person's two local sequences per spot |
| `k2/hapwin.py`, `k2/train_window.py`, `k2/train_spelling.py` | window decisions: where the genotyper's whole-window answer replaces the per-variant one, and how it is spelled |
| `k2/clusterdec.py` | edit-based cluster decoder (GIAB): a cluster's two copies from the reads, Minos's spelling, an arbiter against the per-variant answer |
| `k2/tracks.py`, `k2/assemble.py` | model versions per track (ACTIVE, rollback); scores -> records |
| `k2/util.py` | input and cache paths (`SETTINGS`), bcftools in Docker, VCF helpers |

Released models are code: `k2/models/<track>/ACTIVE` names the version main runs (GIAB g086_20261022T0508, 1000G k060_20261021T1303) and
`k2/models/<track>/<version>/` holds its `model.joblib`, `hom.joblib`, `meta.json` (~3 MB per track). Releasing a new
version = `tracks.promote` (copies those files, plus window / spelling / cluster models when the version has them, here) + commit; rollback = `--rollback` or git revert.
`python -m answer_set.k2.tracks` lists versions.
Inputs, under `SETTINGS.data_dir` (listed in `k2/util.py`): the shared `reference/`, `giab/<donor>.vcf.gz|.bed` and
`1kg/{phased,index,gatk_sites,masks}` (the last three: the builder's 1000G inputs). This clone's `cache/answer_set/k2/`: one build per
round, downloads (Google originals, per-person GATK), mask intervals, trained versions with their leave-out models,
prepared training rounds, the input manifest. Training: `k2/prepare.py` -> `k2/train_track.py` (GIAB: then `k2/train_window.py --base <version>`, then
`k2/train_spelling.py --version <that version>`).

- Needs Docker (the validator's bcftools 1.20 image), internet (S3 DRAGEN, EBI GATK, Google NovaSeq originals),
  `edlib` and `parasail` (requirements.txt). ~1-1.5 min per round on 16 cores.
- Released rounds (`round.released_dir` set) use the leave-out fold model that never saw the round.
- 1000G rounds: Minos flags the inserted mutations `SYNTHETIC` here too (all 45 rounds checked 2026-10-04); k012+ builds
  flag the model's calls (everything not copied by the person's rule) the same way. The validator scores with
  mutations.vcf and ignores the flag, so scores do not depend on it; `main/review.py` / `status.py` use the split.
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
- g042 (released 2026-10-11) = genotyper version 2 (`hapgeno.py`: window edges moved out of tandem repeats, so the
  12-bp anchors sit on unique sequence) + the decision and het/hom models retrained on it with the GATK evidence (108
  rounds, point_alpha 0.5) + the window model retrained on the version-2 windows with the GATK evidence (held out:
  908 windows fixed, 60 broken; g026's: 850 / 212) + g026's spelling model + per-class cut-offs re-chosen by
  validator points with the false-positive cost (SNP 0.95, 1 bp 0.70, 2-3 bp 0.80, 4-7 bp 0.80, >= 8 bp 0.45; window
  0.60). Aimed at Minos's clustered events (an inserted indel with another record within 20 bp; a deletion + an
  insertion close together), where most remaining 1-bp errors are: 28 development rounds, clustered inserted 1-bp
  indels g026 missed whose local DNA the genotyper calls 59 -> 74 of 112. 23 development rounds +0.049 per round vs
  g026 (better 18, worse 5); frozen rounds (from 2026-10-04 08:17, never used for a choice), first 19: +0.058 vs g026
  (t 3.99, 16/3), +0.033 vs g033 (15/4). Builds take ~2-3 min longer (the GATK run).
- g075 (released 2026-10-13) = g054 + the edit-based cluster decoder (`clusterdec.py`). g054 = g048 (decision and
  het/hom models on genotyper version 4, window model with the `lr_answer` feature) with window cut-off 0.75 and per-class
  cut-offs SNP 0.90, 1 bp 0.50, 2-3 bp 0.80, 4-7 bp 0.80, >= 8 bp 0.45 (chosen on 23 development rounds). The cluster
  decoder: Minos puts each inserted mutation on one copy at random and every read it edits for that copy carries all of
  the copy's edits, so per cluster (>= 2 candidates within 20 bp, one an indel) the two local DNA copies are read off the
  reads; every minimal-cost spelling of each copy is listed and a learned chooser picks Minos's (99.7% of copies on
  held-out training rounds); only copies carried by edited reads are written, never the donor's records; where this
  decode and the per-variant answer disagree a learned arbiter decides (held out: 265 clusters fixed, 38 broken in 108
  training rounds). 23 development rounds: +0.0123 vs g042, +0.0027 vs g054 (13 better, 5 worse). Frozen rounds (82, from
  2026-10-04 08:17, never used for a choice): +0.0273 per round vs g042, t 4.54, better 60, worse 22; missed / false
  1-bp indels 611 / 301 -> 470 / 295, missed het SNPs 1120 -> 902. Builds take ~1 min longer.
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
  VQSLOD), plus ~6,700 inserted mutations per round (flagged `SYNTHETIC`). The reads are the person's real public
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
- **Minos changed its 1000G indel insertion on 2026-10-05** (from the 05:53 round on, every 1000G round since): homozygous
  inserted indels now sit in fewer reads (33% below 80% of reads, was 6%; het median 0.45 -> 0.39), as in GIAB rounds;
  SNPs unchanged. k012 then called many of them het: 98.4 / 99.3 / 98.8 on the first three such rounds.
- k020 (1000G, released 2026-10-05) = k012's setup (person's part by the rule) + the haplotype genotyper on 1000G
  rounds too (`features.KG_HAPGENO`, meta `kg_hapgeno`; against the person's public originals), trained on every
  1000G round including the new-pattern ones (1000G data only). Leave-one-out on the 3 new-pattern rounds (each scored
  by a model that never saw it): +0.70 / +0.38 / +0.65 over k012 (1-bp errors 23/19 -> 7/4, 16/11 -> 2/0, 25/14 ->
  11/4). The 46 old-pattern rounds (fold models): +0.017 per round (t 2.95, better 30 / worse 15, worst -0.061).
- k022 (1000G, released 2026-10-06) = k020's setup + three round-level columns (`tracks.ROUND_FEATS`, meta
  `round_feats`, `train_track --round-feats`), computed from the round's own candidate table: among indel candidates
  the read comparison shows as Minos's edits, (share at 0.6-0.85 of the reads) / (share >= 0.85), the het mode's
  median, and each candidate's read fraction relative to it. Since 2026-10-05 05:53 Minos inserts 1000G indels at
  lower read fractions (hom below 0.8: 1-7% of inserted hom indels per round before, 30-39% after; SNPs unchanged), so
  k020 called many hom indels het; the ratio was 0.15-0.49 in all 46 earlier rounds and 0.87-1.01 after, so the round
  itself says which pattern it has (no date switch). Trained on 52 rounds (7 new-pattern). Against k020 on the 10
  new-pattern rounds (6 by fold models, 4 no version had seen): +0.139 per round (t 2.9, better 8/10, worst -0.109);
  vs the same training without the columns (k021) +0.046 on the fold rounds (6/6). 8 old-pattern rounds: -0.012
  (t -0.95).
- k023 (1000G, released 2026-10-06) = k022's recipe trained on 59 1000G rounds (14 new-pattern instead of 7). Scored
  like the validator against k022 (fold models where a version trained on the round): new-pattern rounds +0.072 per
  round (t 4.67, better 14/16, worst -0.025) and on the 4 newest rounds neither had seen +0.12 (4/4); points lost per
  round fall in every indel class (>= 8 bp 0.089 -> 0.070, 1 bp 0.069 -> 0.059, 2-3 bp 0.069 -> 0.057); old-pattern
  rounds +0.008 (t 2.10, better 26 / equal 5 / worse 15, worst -0.080). Ranking 18 near-top GATK configs on 4 rounds:
  0/411 close pairs misordered, as k022.
- k025 (1000G, released 2026-10-07) = k023's recipe (person's part by the rule, hapgeno, round_feats) plus three
  read-comparison columns in that round context (`rd_frac_rel`, `rd_site_n`, `rd_site_sup`; k023 used only the first
  three of `tracks.ROUND_FEATS`): edited-read fraction against the round's het level of edited indels, other edited
  alleles at the position, and edited-read cover of the site. 147 of 153 inserted het indels k023 missed while fewer
  than 10% of the alignments showed them were seen by the read comparison — at 10-20% of the reads (the new pattern's
  low-fraction hets) or as two inserted alleles splitting the reads. Trained on 64 rounds (k023's 59 plus 5 later
  new-pattern); 8 newest new-pattern rounds held out. hap.py vs k023 on those 8 (full models, never in training):
  mean +0.071 per round (t 3.75, better 7/8, worst -0.0005); mean score 99.506 -> 99.577; FN down in every class.
  Six old-pattern fold rounds: +0.009 (t 0.73, better 4 / worse 2). Worst frozen round 20261006T201700Z (HG01953)
  99.123 -> 99.185.
- k028 (1000G, 2026-10-14, not active) = k025's forest with class cut-offs re-chosen on the 40 prepared 1000G
  rounds since 2026-10-07, charging each extra record 3.75/N with that round's mutations.vcf count (the two
  2026-10-14 rounds held out). Cut-offs: snp 0.90, 1 bp 0.60 -> 0.80, 2-3 bp 0.80 -> 0.85, 4-7 bp 0.75 -> 0.90,
  >= 8 bp 0.45 -> 0.15. Mean class points on the 40 rounds moved by about +0.004. On held-out
  20261014T093000Z the new cut-offs removed one extra 1 bp indel and three true 4-7 bp indels (correct 205 -> 202,
  missed 6 -> 9); extras of the other classes did not fall. A support-based one-spelling pass (most
  read-comparison support, then fewest events) rewrote 6 clusters and dropped true alleles, so it is off.
  Future training uses each round's real N (`train_track.extra_cost`).
- k029 (1000G, released 2026-10-14) = k025's recipe retrained on all 109 prepared 1000G rounds (26 Sep through 13 Oct),
  cut-offs chosen with each round's real mutations.vcf count. The two 14 Oct rounds were held out. Cut-offs: SNP 0.90,
  1 bp 0.80, 2-3 bp 0.90, 4-7 bp 0.85, >= 8 bp 0.70. On 20261014T093000Z versus k025, allele by allele: extra
  indels 9 -> 7 (the two 2-3 bp extras gone), inserted exact 7325/7356 -> 7327/7356, 4-7 bp misses 6 -> 4, one
  >= 8 bp indel newly missed. On 20261014T104200Z the extra SNP is gone and one 2-3 bp indel that was the wrong
  het/hom is now exact; one SNP that was the wrong het/hom is now missed. The person's records stay exact on both.
- g077 (GIAB, released 2026-10-16) = g075 + the read guard (`k2/readguard.py`, meta `read_guard`): an inserted indel
  the per-variant model did not call (written by the cluster decoder, window, emit or run-edge step) is dropped when
  >= 6 spanning reads are judged by their DNA and < 15% carry it. Those steps wrote 106 false and 125 true indels in
  31 monitored rounds. Rule fixed before measuring; frozen GIAB rounds 2026-10-13 22:42 to 10-15 23:54 (31), on the
  monitor's saved answers: +0.0075 per round (t 2.96, better 12 / worse 5 / same 14), 22 false and 6 true indels dropped.
- k030 (1000G, released 2026-10-16) = k029 + a repeat het/hom model (`k2/homrep.py`, `hom_repeat.joblib`, meta
  `hom_repeat`): called inserted indels of >= 2 bp in a 1-6 bp tandem repeat spanning >= 8 bp get p_hom from a model
  that adds repeat unit / span / copies, reads judged by their DNA (carrying / reference / other) and CIGAR stutter.
  k029's folds, true inserted repeat indels of its 109 rounds: genotype errors >= 2 bp 258 -> 210 (1 bp worse, so
  1 bp keeps k029's). Release check, 12 live 1000G rounds 2026-10-15 03:30 to 10-16 05:54 (k029 never trained on
  them), both built and scored like the validator: +0.0096 per round, t 2.22, better 5 / worse 2 / same 5
  (`test/answer_set/kg_release_check.py`). Without 10-16 05:54, examined before the design: +0.0067.
- k038 (1000G, released 2026-10-17) = four changes on k029's recipe, each measured held out:
  - base k033: the k025 recipe trained only on the 78 rounds of Minos's current insertion pattern (from 2026-10-05
    05:53); as good as all 124 rounds (k031) with 37% less data;
  - the repeat het/hom model trained on ALL true inserted indels (86k rows), applied as before (repeat >= 2 bp):
    fold genotype errors 472 -> 452;
  - DNA call model (`k2/calldna.py`, `call_dna.joblib`, meta `call_dna`): every non-person indel candidate is
    re-scored by a model that also sees homrep's reads judged by their DNA and the repeat shape; folds, 78 rounds:
    missed / false indels 916 / 183 -> 904 / 162;
  - the cluster decoder (`clusterdec.py`) on 1000G, against the person's originals (`kg_origin`), with a chooser and
    arbiter trained on 1000G clusters (`test/answer_set/kg_spell_data.py`, `kg_clusterdec_train.py`): arbiter cut
    0.8 fixes 85 / breaks 8 clusters in 78 held-out rounds. GIAB's own chooser/arbiter does not transfer (k037).
  Release check: 18 live 1000G rounds 2026-10-16 07:06 to 10-17 07:06, held out from every model and fixed before
  training, built and scored like the validator (`test/answer_set/kg_compare.py`): vs k030 +0.0251 per round, t 2.15,
  better 12 / worse 6 (k036 = without the decoder: +0.0198). Gains come from the weak rounds: on the 9 where k030
  scored below 99.70, +0.049 each (9/9 better); on the other 9 about even.
- k040 (1000G, released 2026-10-17) = k038's stack on a base that also sees the GATK evidence:
  - base k039: the k033 recipe (78 new-pattern rounds) with the fixed-config GATK run's columns (`hcrun.py`; meta
    `kg_hc`, `features.KG_HC`, `prepare.py --kg-hc`; tables in `cache/answer_set/k2/prepared_kghc`). Probe first: on
    3 held-out rounds GATK was right on 24 of 58 inserted indels k036 got wrong. Fold genotype errors 530 -> 449 for
    the base alone; cut-offs snp .9 / 1 bp .8 / 2-3 .8 / 4-7 .65 / 8+ .4;
  - repeat het/hom (all-indel training) and DNA call model re-trained on k039 (folds: genotype 449 -> ~411;
    missed/false indels 872/210 -> 851/202); k038's 1000G cluster chooser/arbiter.
  Release check: 19 live 1000G rounds 2026-10-16 07:06 to 10-17 08:18, held out from every model: vs k038 +0.0350
  per round, t 3.95, better 16 / worse 3; vs k030 +0.0646 (t 4.61); false records 228 -> 202. Builds run GATK
  (~45 s more feature time, median build ~4 min).
- k047 (1000G, released 2026-10-18) = k040's recipe on 102 current-pattern rounds (78 + 24 through 10-17 14:18)
  with a looser read-comparison entry for 1000G indels (features.KG_RD_ENTRY 2, meta kg_rd_entry, prepare
  --kg-rd-entry; tables in cache/answer_set/k2/prepared_kghc_rd2). Found on k040's 18 live rounds: 28 of 36
  "not a candidate" misses had Minos's edit in the read comparison (e.g. 5 of 30 reads) but the entry rule (>= 25%
  of the cover, >= the strongest change within 10 bp) dropped them - Minos now puts inserted indels in fewer reads.
  Rule 2: indels with >= 3 edited reads enter at >= 15% of the cover if >= that neighbour, or at >= 25% regardless.
  On the 102 training rounds it adds 3.7 candidates per round, 1.74 of them real inserted mutations.
  Release check: 15 live 1000G rounds after the training cut-off (2026-10-17 15:30 .. 10-18 09:30), built and
  scored like the validator: vs k040 +0.0311 per round, t 2.42, better 10 / worse 3; vs k045 (same recipe, old
  rule) +0.0255, t 3.28 - the gain is the entry rule; k045 vs k040 +0.0056 (more data alone). GATK fallback k038.
- k050 (1000G, released 2026-10-21) = k047 + a haplotype-aware het/hom model (`k2/hapgt.py`, `hom_hap.joblib`, meta
  `hom_hap`). On k047's 27 live rounds 62% of the wrong het/hom calls on inserted indels (13% of the right ones) had
  another variant in the read window: reads carrying the indel and its neighbour did not match "the indel alone", so
  a hom indel next to a het neighbour looked het, and the reverse. Each spanning read is now judged against every
  local haplotype (the indel with or without each of up to 5 nearest neighbours: the person's records and the other
  called candidates); a model on k047's columns + homrep's measurements + these counts replaces p_hom for every called
  indel candidate that is not the person's. Held-out on k047's folds (102 rounds, 112,114 true inserted indels):
  genotype errors 544 -> 476, with neighbours 306 -> 242, >= 4 bp 193 -> 147.
  Release check: 41 held-out 1000G rounds (k047's out/kg_release set): +0.0218/round, t 3.87, better 22 / worse 5,
  false records 325 -> 291 (test/answer_set/gt_study). Rollback: `--rollback` (k047).
- k060 (1000G, released 2026-10-21) = k050 + window decisions on 1000G (hapwin.py, the GIAB step, now also for 1000G:
  `window.joblib`, meta `window_thr` 0.70; build.py runs it when the person's records are known, with util.Everywhere
  for the missing confident region; GIAB's emit_from_called stays GIAB-only). On k050's 44 live 1000G rounds, at
  1.86 inserted indels per round that never became candidates the local genotyper's called pair of sequences was
  exactly the truth's (compound clusters in repeats, spelled differently from Minos's records); 1000G had no step that
  writes whole windows. Window model trained on k050's 102 rounds (test/answer_set/gt_study/kg_train_window.py,
  leave-out folds; cut-off 0.70: 117 windows fixed / 16 broken). Release check, 41 held-out 1000G rounds: +0.0087/round
  vs k050, t 2.63, better 17 / worse 7, false records 291 -> 284 (cut-off 0.45, the folds' best: +0.0020, t 0.49 -
  0.70 was chosen after seeing that, so a fresh-round check runs in kg_holdout_grow from 2026-10-21 14:00).
- g086 (GIAB, released 2026-10-22) = GIAB retrained on its current pattern + the haplotype-aware call model. Since
  2026-10-20 Minos's newer GIAB tasks insert 11% more indels per Mb and cluster 23% more of them (g077: 99.868 ->
  99.834 per round). Base g082: train_track on the 229 GIAB rounds of 10-01 .. 10-21 10:58 (prepared_giab_v4c: GATK
  evidence, genotyper version 4; float32 tables - 337 rounds as float64 ran out of memory); window model g083 on top
  (train_window --hc --lr, cut-off 0.55: 2,028 windows fixed / 120 broken in 229 rounds); g077's cut-offs, spelling
  model, cluster decoder and read guard; + call_hap.joblib (hapgt.apply_call, person = the donor's NIST records).
  Release check on GIAB rounds from 2026-10-21 11:30 (never trained or tuned on; test/answer_set/giab_fresh_grow.sh),
  9 rounds: +0.0224/round vs g077, t 2.06, better 7 / worse 2 (g085, without the call model: +0.0134; g081, g077 +
  the call model: +0.0068). The combination was chosen after seeing those two on 8 of the 9 rounds. Rollback: g077.
- Off by default, tested and not released: haplotype-aware CALL model (hapgt.apply_call, meta call_hap; 1000G k057
  -0.0007, GIAB g080/g081 +0.009 t 1.74 - fresh-round check running), GIAB haplotype het/hom (g079 -0.0019), entry
  rule 3 (features.KG_RD_ENTRY 3, compound-het second alleles; k056 retrained +0.0063 t 0.97).
- Off by default, not in any released version: `gt_unsure_drop` (leave out inserted indels whose het/hom p is within
  m of 0.5; k048 m 0.1 +0.0055/round, t 2.32 - a scoring tactic, held back) and the GIAB DNA call model hook
  (`call_dna` with `in_mode_donor` as the person; g078 +0.0036/round, t 0.67 on 30 frozen rounds - not released).

## Open work

- **Inserted 1-bp indels** still cause most of the remaining loss (~9 missed + 5.5 extra per round with g013): the
  ~40% of fixable windows the window model leaves alone, sequences noise splits so the genotyper never lists them,
  spots with no window.
- **Spelling**: ~28% of replaced windows are still not spelled exactly like Minos (~17% need spellings no option
  produces yet).
- **1000G**: k029 is the active forest (k025's recipe on 109 rounds, real N in the cut-offs). The held-out gain on
  14 Oct is a handful of alleles. A 1000 Genomes spelling chooser, trained the way the GIAB one was, is still
  unbuilt; picking the spelling with the most read support dropped true alleles (k028).
- **`update()` during a live round** (parked by the user).
