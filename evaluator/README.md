# evaluator

Turns configs into scores exactly like the validator does - as cheaply as possible.

```
EvalJob ──▶ Executor (where/when: local cores or the fleet) ──▶ Runner (GATK) ──▶ CallSet ──▶ Scorer (hap.py) ──▶ ScoreResult
```

Public API (`evaluator/__init__.py`): the factories main uses - `local_executor(workers, batching)`, `scorer()`,
`fleet_executor(inventory)`, `callset_dir(tag)` (`factory.py`: runner settings and the scorer choice live there, not in
main) - and `GatkRunner`, `HappyScorer`, `LocalExecutor`, `SimExecutor`, `FakeRunner`/`FakeScorer`, the protocols in
`base.py`.

## Status (2026-09-30)

| Piece | State |
|---|---|
| `GatkRunner` | the validator template's exact command + safe annotation exclusions; `run`, `run_batch`, `repair`, `variant_cost` |
| GATK batching patch (`gatk_patch/`) | 17 parameters branchable in one run; bit-identical on every test (see below) |
| Divergence repair | diverged batch variants fixed by a repair run (~34 s) instead of a full re-run (132 s) |
| `LocalExecutor` | priority queue, dedup, result cache, batching of queued families, family pricing (`cost_hints`) |
| `HappyScorer` | exact: upstream's hap.py path, one run per callset (the truth sliced once per round) |
| `SimExecutor` | virtual clock for testing strategies and scheduling without GATK |
| Fleet (`fleet/`) | `FleetExecutor` on main + node agents on sub-VPSs: GATK **and scoring** on the nodes, main only proposes |

## Making one evaluation cheaper (kept)

| Speed-up | Effect |
|---|---|
| Exclude annotations the score never uses (`SAFE_ANNOTATION_EXCLUSIONS`) | 1.135x, identical calls. QualByDepth, DepthPerAlleleBySample and RMSMappingQuality must stay (removing them changes QUAL via random draws) |
| 1 thread per GATK job, `floor(cores x 0.95)` jobs per machine | 2.06x throughput vs 4 jobs x 4 threads |
| Truth sliced once per (truth, region) (`truth_slices`) | upstream slices the truth on every scoring (2 bcftools containers); the cached slice is the same bcftools output |
| GATK containers with `--network none` | no network namespace and interface per container (container start 3.9 s -> 0.33 s measured in the legacy project) |

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
102/102, repair 48/48 complete VCF lines on the 3 worst rounds.

**Traps found the hard way:**
- GATK's two generators (java `Random` for QD jitter / reservoir downsampling, `Well19937c` for contamination).
- `AssemblyBasedCallerArgumentCollection` is not Serializable: copying args with Java serialization silently resets
  fields. We use a reflective deep copy plus a deep-equality self-check.
- PairHMM variants need GATK's `realignReadsToTheirBestHaplotype` step.
- `ActivityProfile.add` mutates the state objects it stores: every activity stream needs its own copies.
- `ReservoirDownsampler` takes the same number of draws whatever their values - that is why repair works.

## Fleet (`evaluator/fleet/`)

Main controls every sub-VPS over SSH; the sub-VPSs run GATK **and** score; main only builds the answer set, proposes
configs and submits.

```
fleet.json in the repo folder on main (git-ignored; format in fleet.example.json; re-read at every round start):
  [{"host": "1.2.3.4", "user": "root", "port": 22, "key_filename": "~/.ssh/id_ed25519",
    "key_passphrase": "...", "workers": null, "enabled": true},
   {"host": "5.6.7.8"}]                                # inherits user, key and passphrase

python -m evaluator.fleet.control provision            # fresh machines: copy repo, run scripts/setup_node.sh
python -m evaluator.fleet.control status | sync | selftest | reap
python -m main replay|live ... --fleet                 # mine with the fleet (main starts the node agents)
```

- **Inventory**: `workers: null` = the node's own rule floor(cores x 0.95); `port` lets several machines share
  one IP; optional `repo` / `data_dir`; an entry without `user` / `key_filename` inherits them from the first entry
  that has them (a new host can be just `{"host": ..., "port": ...}`). `key_passphrase` (optional, also inherited):
  main unlocks the keys in its own `ssh-agent`. Each host may have its own key and passphrase.
- **`FleetExecutor`** (main) implements the same Executor protocol, so search and main do not change. For every
  enabled host it keeps one SSH session running the node agent (`python -m evaluator.fleet.node --stdio`); jobs go
  down and scored callsets come back as JSON lines on that session. One queue with dedup and a result cache; a
  **family is always one group** (one GATK run); a node gets groups while it has fewer open groups than workers.
- **Joining**: main pushes the self-test bundle and starts the agent; the node is used only if its **self-test**
  (a fixed config on a 3.6 MB BAM slice, ~11 s) reproduces `fleet/selftest.json`'s digest exactly **and** its code
  digest equals main's (`control sync` fixes a mismatch).
- **Rounds**: main pushes each round's truth files with rsync (and the BAM in replays); on live rounds each node
  downloads the BAM from the presigned URL itself.
- **Failure**: when a session ends (node died, network dropped), that node's unfinished jobs go straight back to the
  queue; main reconnects later. When main stops, every node stops and kills the GATK containers it started.
- **End of round**: after the final submission main calls `end_round`: queued jobs of the round are dropped and every
  node holding its work gets "abort" - it drops the round's local queue, kills exactly that round's containers
  (names carry the round tag) and reports nothing more for it, so all workers are free at the start of the next round.
  `LocalExecutor.end_round` does the same on one machine.
- **Node side**: `LocalExecutor` with `floor(cores x 0.95)` GATK workers (batching, divergence repair, dedup) and
  **one scorer process per GATK worker**, each scoring one callset at a time with upstream's hap.py path
  (`HappyScorer`), exactly what a validator runs. Main counts a group as
  busy until its scored results arrive, so a worker's slot covers GATK and then scoring, and every finished GATK
  result is scored at once. Processes, not threads: the Python half of scoring (upstream's parsers over hap.py's
  output) takes seconds per result on a large GIAB region, and one interpreter (one GIL) capped a 180-core node at
  ~10 results/min. Only scores
  travel back (~1 KB per result), VCFs stay on the node. A scorer process exits by itself within ~2 s when its
  parent is gone, even after `kill -9`.
- **Host size**: the legacy project measured a full evaluation at 159 s on 67-core hosts, 182 s on 178 cores and
  250 s on 356 cores (hap.py 24 -> 37.5 s): big hosts have lower clocks and less memory bandwidth per core. It
  rented 64-core machines. dockerd starts ~16 containers/s, and ~240 concurrent containers made it refuse new ones
  (rc=125); on a 180-core node with ~160 GATK containers plus scoring containers the load average passed 300 with
  30-45% of the CPU idle.
- **Setup a node**: `control provision` copies the repo and **main's reference genomes** (rsync), then runs
  `scripts/setup_node.sh` over SSH: docker, JDK, venv from `requirements.txt`, the GATK image, the **hap.py image
  at upstream's pinned digest** and the **bcftools image** scoring uses (all pre-pulled, so the first scoring does
  not stall or fail), the patch jar, and any genome main did not have (from the platform's reference endpoint).
  Pushing genomes from main runs at ~1.3 MB/s per node over a 35 ms link (~12 min for 907 MB).
- **Blocked registries**: some providers block Docker Hub and quay.io (both on AWS). Docker Hub images then come
  through Google's mirror (`/etc/docker/daemon.json`: `{"registry-mirrors": ["https://mirror.gcr.io"]}`; digests
  are identical, pinned digests resolve); quay.io's bcftools image is copied from main
  (`docker save | ssh <node> docker load`, 66 MB).
- Tested: `tests/test_fleet.py` runs the real node agent as subprocesses (fake GATK): scored results, families as
  one group, a killed node's jobs re-queued to the other, nodes with a failing self-test or other code rejected.
  Real GATK on 2 local nodes from an inventory (no ssh hop): 8 configs incl. a family of 6, the node
  holding all 8 jobs killed at 61 s - all re-queued, finished by the other node, **8/8 scores identical** to the
  single-machine path in 188 s, no agent or container left behind.
  **Real fleet** (2 x 16-core VPS, 28 workers, 2026-09-30): provisioned in ~4 min each, self-test identical on both;
  replay 092900Z, 20 min: **587 configs** (single machine: 296), both nodes 14/14 containers and 86-99% CPU, main
  ~2% CPU; final submission true score 87.53, rank 36/60 (single machine: 87.22, rank 44).
  **Scorer processes vs upstream**: 20 stored callsets scored through spawned
  scorer processes exactly as a node runs them, each compared with upstream's `HappyScorer`: 20/20 identical
  (13 full ScoreResults equal; 1 normal and 5 SNP-less callsets fail in both - hap.py produces no result; 1 SNP-less
  callset scores 0.0 in both, gate failed). A failed callset is reported as failed, never as 0.0.
- Trap: `LocalExecutor.poll` waits for a submit when idle. Callers poll in loops; a poll returning at once would
  spin and starve the scoring threads of the GIL (scoring 20x slower).

## Open work

- **Scale**: more nodes (just `fleet.json` + `control provision`). Nodes keep the last 3 rounds (`KEEP_ROUNDS`) and
  delete older rounds' BAMs and callsets.
- **Big-node fleet run**: scorer processes are verified exact (above) but not yet run on a fleet; a 20-minute replay
  on large nodes should show `busy` equal to the running GATK containers (scoring keeping pace).
- **Shared cache across scorer processes**: each process caches only what it scored itself. If repeats across
  processes turn out to be common, a disk-backed cache per node (SQLite or files keyed by cluster) shares it.
- **Blocked registries, automated**: `setup_node.sh` could set the Docker Hub mirror when a pull fails, and
  `control provision` could copy quay.io's bcftools image from main when its pull fails.
- **Node-side genome download**: nodes fetching their reference from UCSC (`scripts/fetch_genome.sh`) instead of
  main's rsync push would cut provisioning from ~12 min to ~1-2 min.
- **Scoring in the self-test**: the self-test checks GATK output only; a node that cannot score (missing image)
  passes it. Add one scored callset to it.
- **Randomized exactness test** (not run yet): random base configs, random multi-parameter families, deep and
  non-deep rounds. All exactness tests so far used 3 base configs.
- **Family cost cap / splitting**: a batch's time is solo + sum of variant costs; split families that would run
  longer than ~4-5 min and spread them over idle workers.
- **End-of-round scheduling**: in the last minutes, do not start a job whose expected runtime (solo ~132 s at load,
  a batch ~ solo + its variants' costs) exceeds the time left before the final submission - it would only be killed
  by `end_round`. Prefer cheap work then (genotyping variants of the best configs, +1.4 s each).

## Ideas (not planned; decide with data)

- **E4 threaded variant execution**: do the shared prework once, then run each variant's assembly -> genotyping on
  a thread pool in the same JVM. 1 base + 9 assembly variants: ~140 s wall instead of ~780 s, same core time.
  Needs per-thread RNG overrides and reference readers. Worth it only if assembly-heavy families are common.
- **Repair stage 2**: run the repair only on small windows around deep spots (seconds instead of ~34 s); needs
  starting GATK mid-region at the exact generator position.
- **Smarter batch main**: divergence depends on the batch's main config; learn per round which groupings stay exact.
- **S1 additive scorer**: per-cluster score contributions so a callset's score is a sum over clusters.
