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
  main: baseline first, then keeps the best config submitted (at most every 5 min, forced once before the deadline)
```

- Strategies never know where jobs run; the executor never knows why a job was chosen.
- Running and scoring are separate steps: if the answer set changes, callsets are re-scored, never re-run.
  Today scoring happens on the main machine (background thread); with a fleet it must move to the nodes
  (see evaluator/README.md, open work).
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
docker pull quay.io/biocontainers/samtools:1.20--h50ea8bc_0   # parent-read slices
# hap.py image: the digest pinned in minos_subnet/utils/scoring.py is pulled on first use

bash evaluator/gatk_patch/build.sh              # GATK batching patch -> evaluator/gatk_patch/build/minos-patch.jar
bash scripts/fetch_reference_data.sh            # 1000G callsets/panels/masks + CRAM indexes into $DATA_DIR
.venv/bin/python -m backtest.data --last 20     # released rounds (BAM + truth) for backtests
```

Data never goes in git; it lives under `DATA_DIR` (default `/minos/data`). The reference FASTAs + SDFs per
chromosome (`$DATA_DIR/reference/chrN/`) come from the platform's reference endpoint, and the GIAB NIST VCFs + BEDs
(`$DATA_DIR/giab/`) from NIST; neither is scripted yet.

## Everyday commands

```bash
.venv/bin/python -m unittest discover -s tests -t .              # unit tests (~7 s); RUN_SLOW=1 adds GATK tests
.venv/bin/python -m main replay 20260928T092900Z --minutes 20     # whole miner on a stored round, dry-run
.venv/bin/python -m main live --wallet W --hotkey H --demo        # live miner (untested; demo mode first)
```

`replay` prints what the miner would have submitted, its **true** score and its rank among the real submissions.
Reference point (2026-09-30, one 16-core machine, 20 min): 296 configs evaluated, true 87.22, rank 44 of 60.
