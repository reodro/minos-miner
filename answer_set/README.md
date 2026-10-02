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

## Open work

- **Inserted 1-bp indels** cause ~60-70% of the remaining score loss: ~3 in 4 errors never become candidates
  (the event sits in a long homopolymer or short repeat and the aligner spells it differently). A local haplotype
  step is parked.
- **`update()` during a live round** (parked by the user).
