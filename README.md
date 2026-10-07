# minos-miner

Miner for Bittensor subnet 107 (Minos, genomic variant calling). Each round the platform publishes a BAM and a
region; miners submit a **GATK HaplotypeCaller config**, validators run it and score the calls against a hidden
truth. A round lasts **72 minutes** and **only the last submission before the deadline counts**.

Our approach, in one line: rebuild the round's truth ourselves (the **answer set**), evaluate thousands of configs
against it as cheaply as possible, search well, and keep the best config submitted.

## The parts

| Package | Owner | Job | Start here |
|---|---|---|---|
| `common/` | everyone | Contracts (`types.py`), the upstream score formula wrapper, parameter ranges, platform client, result store, settings | `common/types.py` |
| `answer_set/` | answer-set owner | Rebuild the round's truth at round start: donor, germline, synthetic pool | [answer_set/README.md](answer_set/README.md) |
| `evaluator/` | evaluator owner | Run configs (GATK, batched), schedule jobs, score callsets exactly like the validator | [evaluator/README.md](evaluator/README.md) |
| `search/` | search owner | Strategies that decide which configs to try (RF_UCB, random, ...) | [search/README.md](search/README.md) |
| `main/` | main owner | The miner itself: round loop, fallback chain, live platform, submission, replay | [main/README.md](main/README.md) |
| `backtest/` | everyone | Released rounds, real field scores, KPI functions | `backtest/kpi.py` |

Rules for working together are in [CONTRIBUTING.md](CONTRIBUTING.md) - read them before your first PR.

## How a round runs

```
round opens ─▶ main: download BAM, announce round
                 │
                 ▼
          answer_set.build ──▶ AnswerSet (germline exact + synthetic pool + N estimate, truth VCFs)
                 │
                 ▼
  search strategies ──propose/propose_many──▶ EvalClient ──EvalJob──▶ Executor ──▶ GatkRunner (batched)
         ▲                                     │                                      │
         └──────────── Outcome (score) ◀── Scorer (hap.py, background thread) ◀── CallSet
                 │
                 ▼
  main: baseline + search start together; the best config is kept submitted (at most every 5 min, forced once
  before the deadline)
```

- Strategies never know where jobs run; the executor never knows why a job was chosen.
- Running and scoring are separate steps: if the answer set changes, callsets are re-scored, never re-run.
  On one machine scoring runs in a background thread; with a fleet (`--fleet`) main controls the sub-VPSs over SSH,
  they run GATK **and** score (in scorer processes), and main only proposes and submits (see evaluator/README.md).
- After the round, `main.review` judges every part against the released truth (answer-set accuracy, true scores,
  pick gap, rank). A live dashboard on main, behind a login, shows the round and edits the search's settings
  (dashboard/README.md); the settings live in `common/tuning.py`, with overrides in a git-ignored `tuning.json`.
- Related configs (a *family*: same config except cheap parameters) go to the executor together and run as **one**
  GATK job; each still gets its own exact result.

## How the score works (v2, from the validator code)

`score = 100 x (0.70 x core + 0.30 x exp(-fp_per_target / 8))`, zeroed if ti/tv or het/hom are implausible.

- **core** = difficulty-weighted F1 over **every** truth record in the region (germline + all synthetic), per class:
  SNP het 0.18, SNP hom 0.02, indel 1 bp **0.40**, 2-3 bp 0.16, 4-7 bp 0.10, >= 8 bp 0.14.
- **fp_per_target** = false positives anywhere in the region / **N**, where N = number of entries in
  `mutations.vcf` (~370-790). N is the *only* thing the target list is used for, and it is **not known live**.
- Never reimplement the formula: `common/scoring.py` wraps the upstream `minos_subnet/utils/scoring.py`.

## What we know about the rounds (measured 2026-09-30)

- **Truth** = one donor's germline + a fresh pool of ~5K injected variants flagged `SYNTHETIC` (~25-35% of truth).
- **Germline source**: GIAB NIST v4.2.1 calls (HG001-HG007) or one 1000 Genomes 30x sample (NYGC raw callset +
  a fixed filter). Both regimes are reconstructed exactly (1000G 36/36 rounds, GIAB 24/24).
- **Reads**: 128 bp, hashed names, but they are **real sequencing of the donor**: real variants that the truth
  leaves out are visible and count as false positives for any caller.
- **Timing**: one GATK run is ~72 s on an idle core, ~132 s at full load (14 jobs on 16 cores).
- **Live policy** (network config): on-chain config commitment **on**, resubmission fee **off**.

## Setup (one machine)

```bash
# system: Python 3.10+, docker, java 17 (for the GATK patch build), bcftools, curl
git clone https://github.com/minos-protocol/minos_subnet.git /minos/minos_subnet   # upstream: formula, templates
git clone <this repo> /minos/minos-miner && cd /minos/minos-miner
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt

docker pull broadinstitute/gatk:4.5.0.0                       # validator's GATK image
docker pull "$(grep -oP 'HAPPY_DOCKER_IMAGE = "\K[^"]+' ../minos_subnet/utils/scoring.py)"      # hap.py, pinned digest
docker pull "$(grep -oP 'BCFTOOLS_DOCKER_IMAGE = "\K[^"]+' ../minos_subnet/utils/scoring.py)"   # bcftools for scoring

bash evaluator/gatk_patch/build.sh              # GATK batching patch -> evaluator/gatk_patch/build/minos-patch.jar
bash scripts/fetch_genome.sh                    # GRCh38 chr1-22: FASTA, index, dict, SDF (validator's files)
bash scripts/fetch_reference_data.sh            # 1000G phased panel (chr14-22 by default) + GIAB index
.venv/bin/python -m backtest.data --last 20     # released rounds (BAM + truth) for backtests
```

Data never goes in git, and lives in three places (`common/settings.py`):
- `DATA_DIR` (default `/minos/data`): **inputs** - reference, 1000G, panel, GIAB, rounds, truths (released truths of
  rounds without a stored BAM). Shareable between clones.
- `cache/` in each clone: **rebuildable caches** (callsets, results, parent slices, truth slices).
- `status/` in each clone: what this clone's miner writes for the dashboard and the review.

Git-ignored folders in each clone (never committed, never synced to fleet nodes - every clone keeps its own):

| Folder | What | Safe to delete? |
|---|---|---|
| `cache/` | rebuildable caches (callsets, results, parent slices, truth slices, the fleet self-test's BAM slice) | yes - rebuilt on the next run |
| `status/` | this clone's miner output: live round, history, review (the dashboard reads it) | no - it is the run history; back it up |
| `test/` | your experiments, benchmarks, exactness checks and their output | your call - nothing in the miner uses it |

The GIAB NIST VCFs + BEDs (`$DATA_DIR/giab/`) are not scripted yet. Fleet nodes need much less: `scripts/setup_node.sh` sets one up.

## Everyday commands

```bash
.venv/bin/python -m unittest discover -s tests -t .              # unit tests (~75 s); RUN_SLOW=1 adds GATK tests
.venv/bin/python -m main replay 20260928T092900Z --minutes 20     # whole miner on a stored round, dry-run
.venv/bin/python -m main live --wallet W --hotkey H --demo        # live miner (untested; demo mode first)
```

The live dashboard (https://<main's IP>:8443, behind a login) runs on main: see [dashboard/README.md](dashboard/README.md).

`replay` prints what the miner would have submitted and reviews the run: answer-set accuracy, the **true** scores
of the final, baseline and top configs, and the rank among the real submissions.
Reference points (20-minute replays): round 092900Z - one 16-core machine 296 configs, true 87.22, rank 44/60;
fleet of two 16-core sub-VPSs 587 configs, true 87.53, rank 36/60. Round 20261001T044100Z (GIAB HG007, chr18) -
fleet of a 32-core and a 180-core sub-VPS: true 78.47, rank 52/68 (field best 79.77).

## Running a fleet test (runbook)

The steps, in order, as run for the engine v4 test on 2026-10-05 (6 nodes, 823 workers). Long steps run in the
background; nothing here edits `fleet.json`, which is the user's file.

1. **Nodes.** The user rents them and lists them in `fleet.json` in the repo folder (git-ignored): host, port,
   user, key, `key_passphrase`.
2. **Clear stale host keys.** Providers reuse addresses, and an old key makes ssh refuse the new machine:
   ```bash
   ssh-keygen -R <ip>                    # port 22
   ssh-keygen -R "[<ip>]:<port>"         # a host on another port
   ```
3. **Provision** new nodes, in parallel, about 3-6 minutes each. It installs Docker, the GATK image and the
   genomes, copies the repo and builds the GATK patch jar there.
   ```bash
   .venv/bin/python -m evaluator.fleet.control provision <host> ...    # no host: every host in fleet.json
   .venv/bin/python -m evaluator.fleet.control status                  # cores, workers, image, jar, genomes, code
   ```
   A node gets min(cores - cores // 20, available RAM / 1.3 GB) workers. A node with little RAM per core gets few workers
   (98 GB -> 73) and still runs out of memory under load (5x slower in the v4 test): rent at least 4 GB per core.
4. **On main, before the run:** rebuild the patch after a patch change (`bash evaluator/gatk_patch/build.sh`), run
   the unit tests, and set the dashboard password once (`.venv/bin/python -m dashboard passwd`).
5. **Run it** in the background; a 72-minute search plus the truth review takes about 80 minutes:
   ```bash
   scripts/fleet_replay.sh 20261002T092900Z 72     # sync the code, wipe the round's node caches, replay, review
   ```
   To check that everything runs without the fleet (not the score): `.venv/bin/python -m main replay <tag>
   --minutes 15` on main's own workers.
6. **Watch** on the dashboard, `https://<main's IP>:8443`. The Settings tab writes `tuning.json`: live settings
   apply within a second, round settings at the next round, and the run's manifest records every change. The
   run's log is `test/runs/replay_<tag>_<HHMM>.log`; `grep "engine: "` shows the engine's decisions.
7. **Result.** The end of the log and the dashboard's Result panel show the submitted config's true score, the
   pick gap, the local-vs-true rank correlation and the rank. To check more configs against the released truth:
   ```bash
   .venv/bin/python scripts/true_scores.py <run_id> --top 7     # the top 7 distinct local results
   ```
8. **Pull the history, then release the nodes.** The VCFs exist only on the nodes:
   ```bash
   .venv/bin/python scripts/pull_history.py <tag> <tag>_<label> test/runs/replay_<tag>_<HHMM>.log
   ```
   It lands in `$DATA_DIR/history/<label>/` (about 30-40 GB, 5-7 minutes). Then the user turns the nodes off.
9. **Record it:** add the run to the engine versions table in [search/README.md](search/README.md).

If main was killed mid-run, node agents can outlive it: the wipe step prints `AGENT-RUNNING` for such a node, and
`.venv/bin/python -m evaluator.fleet.control reap` stops their GATK containers once no agent is running.
