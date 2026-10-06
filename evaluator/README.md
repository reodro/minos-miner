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
| `DecisionCacheScorer` | exact (identical `ScoreResult`s to `HappyScorer`): hap.py runs only for clusters this machine has not scored yet; the cluster cache is a per-round SQLite store shared by every scorer process on the machine (`evaluator.scorer()` returns it) |
| `HappyScorer` | the reference: upstream's hap.py path, one full run per callset (the truth sliced once per round) |
| `SimExecutor` | virtual clock for testing strategies and scheduling without GATK |
| `ScoringPool` | main's own scoring (no fleet) in one scorer process per GATK worker sharing the decision cache - what a node does for its results; one thread scores ~6 results/min, 14 GATK workers make ~17 |
| Fleet (`fleet/`) | `FleetExecutor` on main + node agents on sub-VPSs: GATK **and scoring** on the nodes, main only proposes |

## Making one evaluation cheaper (kept)

| Speed-up | Effect |
|---|---|
| Exclude annotations the score never uses (`SAFE_ANNOTATION_EXCLUSIONS`) | 1.135x, identical calls. QualByDepth, DepthPerAlleleBySample and RMSMappingQuality must stay (removing them changes QUAL via random draws) |
| 1 thread per GATK job, `floor(cores x 0.9)` jobs per machine | 2.06x throughput vs 4 jobs x 4 threads |
| Machine-wide decision cache (`DecisionCacheScorer`) | a cluster scored by any process on the node is never sent to hap.py again by another; 34/34 identical to plain hap.py (3 rounds incl. deep coverage) |
| One scoring unit per GATK batch (node agent) | a batch's callsets are scored one after another by one scorer (`score_many`: uncached clusters packed into few hap.py runs), so a worker slot never runs more than one scoring - no burst of parallel hap.py runs when a family finishes |
| hap.py CPU limit (`MINOS_SCORE_CPUS`, default 1; 0 = none) | upstream runs hap.py with `--threads 4` and no limit; on a 180-core node scoring bursts drove the load to ~1000. At 1 core a full scoring takes ~10 s instead of ~5 s idle, identical results |
| Identical records scored once (node agent) | a callset whose records digest was already scored this round reuses that score (or failure): no hap.py |
| Truth sliced once per (truth, region) (`truth_slices`) | upstream slices the truth on every scoring (2 bcftools containers); the cached slice is the same bcftools output |
| GATK containers with `--network none` | no network namespace and interface per container (container start 3.9 s -> 0.33 s measured in the legacy project) |
| Provably identical configs run once (`equivalence.py`) | a config is run as its canonical config when the round's BAM cannot tell them apart: `max_reads_per_alignment_start` >= the most reads at one start, `min_base_quality_score` / `base_quality_score_threshold` with no compared quality between them (BAM qualities + 0 and 20 from overlapping mates; MAPQ cap). Result delivered under the asked-for key. Stock GATK digests identical for every class tested (2 datasets) |
| GATK servers: persistent JVMs (`gatk_patch` MinosServer, `GatkRunner` server pool) | a worker's JVM runs job after job, so only its first job pays class loading and JIT warm-up: stock GATK full round on one core 162 s cold vs 119-127 s warm; on the self-test slice a sequence of 8 jobs took 163 s vs 334 s cold. Before/after each job every static the patch keeps is reset (`Utils.minosResetAll`: generators, overrides, activity hook, cap report). Server vs cold: identical outputs (see gatk.py). Every 40th solo server run is also run cold and compared (`MINOS_SERVER_VERIFY`); a failing server is replaced and the job runs cold; a server leaves by itself after 10 min without a job (a finished review's executor would otherwise leave 8 idle JVMs per round); `MINOS_GATK_SERVER=0` = off |
| Learned haplotype cap | the patch reports per output the smallest `max_num_haplotypes_in_population` that reproduces it (`.minos.kbest`); a run no k-best search of which hit its cap is reused for every cap >= that. Exact, but only fires when no region reaches the cap (the slice: 72; round 20261004T092900Z: the default 128 is reached somewhere, so nothing is learned) |
| Activity cache (`-Dminos.activity.dir`, `MINOS_ACTIVITY_CACHE=0` = off) | the per-locus activity states (pileups + isActive, ~50% of a run's Java CPU) depend only on the skeleton (`branch_key`), so the first run of a skeleton saves each shard's states under `<callsets>/<tag>/activity/<skeleton>/` and later runs read them; reads are still collected per region. A shard is saved only if its read pass drew no random numbers (deep spots); a cached shard that would draw fails the run. VCF records identical to plain runs (assembly, activity, PairHMM, genotyping variants and a deep-spot config; `test/activity_cache_check.py`). Warm servers 62 s -> 44 s per repeated skeleton (single configs); batches of ~6 configs gain less (202 s -> 176 s under load); 60-min replay of the run-5 code 16.9 -> 17.5 configs/min. ~22 MB per skeleton per round (chr21 slice) |
| Assembly cache (`-Dminos.assembly.dir` + `-Dminos.assembly.skeleton`, `MINOS_ASSEMBLY_CACHE=0` = off) | a region's assembly depends on the skeleton, the branchable assembly settings, the region and its finalized reads - not on PairHMM or genotyping settings - so `MinosAssemblyCache` saves each region's haplotypes (in order, every field later steps read) under `<callsets>/<tag>/assembly/`, keyed by those plus every read's SAM text, and a later run rebuilds the result on GATK's own reference-only result (`replaceAllHaplotypes`, as the assembler finishes) instead of assembling; main config and batch assembly riders alike. The k-best searches it reported to `MinosKBest` are replayed. VCF records identical to plain runs (PairHMM / genotyping / assembly / activity variants, deep spots, a batch; `test/assembly_cache_check.py`). A config differing from an earlier one only in PairHMM / genotyping: 116 s plain, 103 s activity cache, 68 s both (cold, one core). ~0.15 KB per region, but many files (52k for 9 configs) |

Scheduling on this machine's real speed: every cold run records real seconds / table seconds; the deadline uses the
90th percentile of the last 30 (`stats()["scale"]`, ~2.07 on a 6-core i7 with 10 workers), for a whole family, while
packing stays in table seconds. `stats()["timeline"]` (and the log) records each minute's busy/idle workers, queue,
held jobs, callsets being scored and what search submitted - the first thing to read when workers sit idle.

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
- A batch that dies inside GATK (rc=3: an exception, not a bad argument) has its **head rerun alone** and its riders
  fail. The batch path throws `Graph must have ref source and sink vertices` in the main's assembly for some heads
  on some regions (every crashing head so far had `min_dangling_branch_length` 19-20; padding 0-1 heads crashed in
  plain mode only, others warm-style too) while every such config passes on its own. The cause inside the patch is
  open (`run_batch(..., warm=True)` runs warm-style without keeping the JVM, kept for experiments); the search also
  never proposes padding below 5 (never in the top 200 of 7,067 true-scored configs).

### Batch size: small vs big

Every job says how big its batch may get (`EvalJob.batch`, set by the search through
`client.propose_many(..., batch=...)`):

- **big** (default): packed with other queued configs that share its GATK traversal, up to 32 configs and
  `evaluator.big_batch` = 2.3 solo runs of work in total, run on the node with the largest free share. Wave 0, the scout,
  riders: nothing waits on one of them.
- **small**: a result the search waits for (a climbing move, a merge child, a generation child, a midpoint test). It
  is always the head of its own batch, holds at most `evaluator.small_batch` = 1.3 solo runs of work in total (itself plus
  cheap riders such as its genotyping and PairHMM siblings), never rides in another batch, and goes to the fastest
  free node that has its round, or to the fastest free slow node when no fast one has a free worker. Packing
  structural variants of one centre in one batch runs them one after another on one worker - 30 assembly variants
  take ~16 run-times - which is why a waited-for result must stay small.

Both limits are total work, the head included, in the node's own solo time. They are settings in
`common/tuning.py`, changeable on main during a round; `batch_limit()` and `max_batch_s()` in `evaluator/local.py`
read them for the local executor, the coordinator and the simulator alike. Nodes always use the defaults, which only
matters for a node's own re-packing: the groups it runs are formed on main.

Each node's solo-run time is measured from its finished batches (an EMA of wall time / predicted cost units;
`FleetExecutor.nodes()` shows `solo_s`) and scales that node's batch cap and rider costs; the planning constant
`SOLO_S` (132 s) is only the starting value.

A late small batch is not waited for. Width is cheap and latency is not, so while the fleet has a worker to spare
(nothing queued could use it) a small batch that has run longer than 1.5x its expected time on its node, or that a
free faster node could redo before its node is expected to finish, is duplicated there - same head, same riders.
The first result wins; the other node kills its run (`abort_jobs`: the executor's `cancel_running`, the runner's
`kill_configs`) and reports nothing for it. `FleetExecutor.stats()` counts `speculated` and `speculation_won`.

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
python -m evaluator.fleet.control genomes              # push main's reference genomes (endpoint unreachable there)
python -m evaluator.fleet.control status | sync | selftest | reap
python -m main replay|live ... --fleet                 # mine with the fleet (main starts the node agents)
```

- **Inventory**: `workers: null` = the node's own rule floor(cores x 0.9); `port` lets several machines share
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
- **Rounds**: `publish_round` ships the round itself first - main pushes the BAM with rsync in replays, on live
  rounds each node downloads it from the presigned URL - and a node starts GATK the moment it has it, while main is
  still building the answer set. `publish_answer_set` ships the truth files when the answer set exists (and every
  newer version): a node holds the results it made before that and scores them the moment the truth lands. A
  chromosome a node has not got is fetched from the platform's reference endpoint on demand.
- **Failure**: when a session ends (node died, network dropped), that node's unfinished jobs go straight back to the
  queue; main reconnects later. When main stops, every node stops and kills the GATK containers it started.
- **A node on main's own machine** (`{"host": "local"}` in the inventory; the bench runs this way: `--fleet
  test/engine/fleet_local.json`) runs the real agent without ssh and gets `MINOS_NODE_SHARED_DATA=1`: the rounds it
  forgets are main's stored rounds, so it never deletes a round directory.
- **End of round**: after the final submission main calls `end_round`: queued jobs of the round are dropped and every
  node holding its work gets "abort" - it drops the round's local queue, kills exactly that round's containers
  (names carry the round tag) and reports nothing more for it, so all workers are free at the start of the next round.
  `LocalExecutor.end_round` does the same on one machine.
- **Node side**: `LocalExecutor` with `floor(cores x 0.9)` GATK workers (batching, divergence repair, dedup) and
  **one scorer process per GATK worker**, each scoring one GATK batch's callsets at a time, one after another, with
  `DecisionCacheScorer` (upstream's
  hap.py path for the clusters not yet scored on this node; the per-round cluster store in `cache/dcache/<tag>/` is
  shared by all of them and deleted with the round). Before a callset reaches a scorer, the agent checks its records
  digest: identical records already scored this round reuse that score. Main counts a group as
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
- **Setup a node**: `control provision` copies the repo, then runs `scripts/setup_node.sh` over SSH: docker
  (Ubuntu's, unless the image already ships one - vast.ai VMs do), JDK, venv from `requirements.txt`, the GATK
  image, the **hap.py image at upstream's pinned digest** and the **bcftools image** scoring uses (all pre-pulled,
  so the first scoring does not stall or fail), the patch jar, and the reference genomes **chr14-22** (every round
  so far was on one of them; any other is fetched on demand) from the platform's reference endpoint. `control
  genomes` pushes main's copies instead for a host that cannot reach the endpoint - at ~1.3 MB/s per node over a
  35 ms link that is ~17 min for the nine chromosomes, which is why provisioning does not do it.
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
- **Scoring in the self-test**: the self-test checks GATK output only; a node that cannot score (missing image)
  passes it. Add one scored callset to it.
- **Randomized exactness test** (not run yet): random base configs, random multi-parameter families, deep and
  non-deep rounds. All exactness tests so far used 3 base configs.

## Ideas (not planned; decide with data)

- **E4 threaded variant execution**: do the shared prework once, then run each variant's assembly -> genotyping on
  a thread pool in the same JVM. 1 base + 9 assembly variants: ~140 s wall instead of ~780 s, same core time.
  Needs per-thread RNG overrides and reference readers. Worth it only if assembly-heavy families are common.
- **Repair stage 2**: run the repair only on small windows around deep spots (seconds instead of ~34 s); needs
  starting GATK mid-region at the exact generator position.
- **Smarter batch main**: divergence depends on the batch's main config; learn per round which groupings stay exact.
- **S1 additive scorer**: per-cluster score contributions so a callset's score is a sum over clusters.
