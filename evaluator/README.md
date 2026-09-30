# evaluator

Turns configs into scores exactly like the validator does - as cheaply as possible.

```
EvalJob ──▶ Executor (where/when: local cores, later a fleet) ──▶ Runner (GATK) ──▶ CallSet ──▶ Scorer (hap.py) ──▶ ScoreResult
```

Public API (`evaluator/__init__.py`): `GatkRunner`, `HappyScorer`, `DecisionCacheScorer`, `LocalExecutor`,
`SimExecutor`, `FakeRunner`/`FakeScorer`, and the protocols in `base.py`.

## Status (2026-09-30)

| Piece | State |
|---|---|
| `GatkRunner` | the validator template's exact command + safe annotation exclusions; `run`, `run_batch`, `repair`, `variant_cost` |
| GATK batching patch (`gatk_patch/`) | 17 parameters branchable in one run; bit-identical on every test (see below) |
| Divergence repair | diverged batch variants fixed by a repair run (~34 s) instead of a full re-run (132 s) |
| `LocalExecutor` | priority queue, dedup, result cache, batching of queued families, family pricing (`cost_hints`) |
| `DecisionCacheScorer` | exact hap.py scoring with a per-cluster cache; `score_many` packs several callsets into shared runs |
| `SimExecutor` | virtual clock for testing strategies and scheduling without GATK |
| Fleet (`fleet/`) | `FleetExecutor` on main + node agents on sub-VPSs: GATK **and scoring** on the nodes, main only proposes |

## Making one evaluation cheaper (kept)

| Speed-up | Effect |
|---|---|
| Exclude annotations the score never uses (`SAFE_ANNOTATION_EXCLUSIONS`) | 1.135x, identical calls. QualByDepth, DepthPerAlleleBySample and RMSMappingQuality must stay (removing them changes QUAL via random draws) |
| 1 thread per GATK job, `floor(cores x 0.95 - 1)` jobs per machine | 2.06x throughput vs 4 jobs x 4 threads |
| Decision cache (per-cluster hap.py results) | re-scoring 7.2 s -> ~0.6 s, identical |
| Packed scoring (`score_many`) | 2.3-3.1x, identical |

Rejected after measuring (not exact, or no gain): GVCF re-genotyping, QUAL projection, native Java outside the
image, interval restriction, scatter-gather, JVM heap/GC tuning, uncompressed BAM, no phasing.

## Batching: one GATK run, many results

A GATK run is: read pass -> regions -> assembly -> PairHMM -> genotyping. Configs that differ only in a later stage
share everything before it. The patch (`gatk_patch/src`, shadow classes loaded before the official 4.5.0.0 jar)
runs the extra configs ("variants") inside one run and writes one VCF per variant.

| Variant differs only in | Params | Extra cost per variant (full load) |
|---|---|---|
| Genotyping | stand_call_conf, indel_heterozygosity, max_alternate_alleles, contamination | 1.4 s (~1%) |
| PairHMM | base_quality_score_threshold, pcr_indel_model, gap continuation penalty, mismapping rate | 16 s (~12%) |
| Assembly | min_pruning, min_dangling_branch_length, recover_all_dangling_branches, max haplotypes, adaptive pruning error rate, pruning LOD | 72 s (~55%) |
| Activity | active_probability_threshold, min/max_assembly_region_size | 78 s (~59%) |
| anything else (read filters, padding, kmers, GVCF modes, ...) | - | full run, 132 s |

- Costs live in `gatk.VARIANT_COST`; a variant pays the most expensive stage it changes.
- A batch is one job on one core: many assembly/activity variants in one batch finish late (latency).
- `assembly_region_padding` is not branchable (it also sets read-shard padding); GVCF modes never batch.

### Exactness and randomness (read this before touching the patch)

GATK draws random numbers from **fixed-seed** generators, in order, for two things: QD jitter on calls with
QD >= 35 (genotyping), and random read selection at **deep spots** (positions with more than
`max_reads_per_alignment_start` reads starting there). In a batch the read pass is shared, but each variant's own
run would have taken different numbers at a deep spot if it made a different number of QD jitters before it.

- Each variant gets its own copy of both generators; read-pass draws are **replayed** into it. If a replayed value
  differs from what its own run would draw, the variant is flagged **diverged**.
- A diverged variant is **repaired**: its own run skips every region whose reads equal the batch's (logged per
  region: span, read fingerprint, draws) and recomputes only regions touched by deep spots; the output is spliced.
  A position check at every skip falls back to computing everything if anything is off.
- Measured: rounds with deep spots <=> rounds with divergence (59/60). 19 of 60 stored rounds diverge, rising to
  ~55% of the newest; a diverging round typically has only 1-4 deep positions.

Tests (all bit-identical to separate stock runs): genotyping 54/54, PairHMM 62/62, assembly 86/86, activity
102/102, repair 48/48 complete VCF lines on the 3 worst rounds. Scripts: `test/tip_*branch.py`, `test/tip_repair.py`.

**Traps found the hard way:**
- GATK's two generators (java `Random` for QD jitter / reservoir downsampling, `Well19937c` for contamination).
- `AssemblyBasedCallerArgumentCollection` is not Serializable: copying args with Java serialization silently resets
  fields. We use a reflective deep copy plus a deep-equality self-check.
- PairHMM variants need GATK's `realignReadsToTheirBestHaplotype` step.
- `ActivityProfile.add` mutates the state objects it stores: every activity stream needs its own copies.
- `ReservoirDownsampler` takes the same number of draws whatever their values - that is why repair works.

## Fleet (`evaluator/fleet/`)

```
main:   python -m main replay|live ... --fleet 0.0.0.0:8707 --cert fleet.crt --key fleet.key
node:   python -m evaluator.fleet.node --main https://MAIN:8707 --cacert fleet.crt
both:   MINOS_FLEET_TOKEN=<shared secret>   (env only - never in the repo or a config file)
```

- **`FleetExecutor`** (main) implements the same Executor protocol, so search and main do not change. Jobs wait in a
  queue (dedup + result cache); a node asks for as many groups as it has free workers; a **family is always one
  group** (one GATK run). Every hand-out is a **lease**: a node that stops heartbeating loses it after `lease_s`
  and its unfinished jobs go back to the queue; a late result from an expired lease still counts, duplicates are
  ignored. Main serves each round's truth files (and the BAM when no presigned URL is given).
- **Node agent** (sub-VPS): before joining, a **self-test** runs a fixed config on a 3.6 MB BAM slice and must
  reproduce `fleet/selftest.json`'s digest exactly (~11 s); then it pulls work, runs it with `LocalExecutor`
  (batching, divergence repair, dedup), **scores every callset on the node** (DecisionCacheScorer, 2 threads) and
  returns callsets carrying their score. VCFs stay on the node (~1 KB per result over the network).
- **Security**: every request carries the token; off localhost main refuses to start without TLS
  (`scripts/make_fleet_cert.sh`) and nodes refuse plain http.
- **Setup a node**: `scripts/setup_node.sh` (docker, JDK, venv, images, patch jar, reference genomes).
- Tested (`tests/test_fleet.py`, `test/fleet_local.py`): 2 node processes on one machine, 8 configs incl. a family
  of 6 - 8/8 scores identical to the single-machine path; a node killed while holding work had its 7 jobs re-queued
  and finished by the other; results in 130 s. **Not yet run on real separate VPSs.**
- Trap found: `LocalExecutor.poll` used to return at once when idle; callers poll in loops, and the spin starved the
  scoring threads of the GIL (scoring 20x slower). An idle poll now waits for a submit.

## Open work

- **Real multi-VPS run**: rent 2 small nodes, run `setup_node.sh`, a replay with `--fleet`; watch per-node
  throughput, scoring lag and bandwidth (BAM download per round per node).
- **Randomized exactness test** - written, not run: `test/tip_random_exactness.py` (random base configs, random
  multi-parameter families, deep and non-deep rounds). All exactness tests so far used 3 base configs.
- **Scoring threads per node**: 2 by default; at ~59 workers and families of 5, a node may produce ~2 callsets/s -
  measure whether 2 threads keep up, raise if the scoring queue grows.
- **Family cost cap / splitting**: a batch's time is solo + sum of variant costs; split families that would run
  longer than ~4-5 min and spread them over idle workers.

## Ideas (not planned; decide with data)

- **E4 threaded variant execution**: do the shared prework once, then run each variant's assembly -> genotyping on
  a thread pool in the same JVM. 1 base + 9 assembly variants: ~140 s wall instead of ~780 s, same core time.
  Needs per-thread RNG overrides and reference readers. Worth it only if assembly-heavy families are common.
- **Repair stage 2**: run the repair only on small windows around deep spots (seconds instead of ~34 s); needs
  starting GATK mid-region at the exact generator position.
- **Smarter batch main**: divergence depends on the batch's main config; learn per round which groupings stay exact.
- **S1 additive scorer**: per-cluster score contributions so a callset's score is a sum over clusters.
