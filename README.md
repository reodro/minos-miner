# minos-miner

Miner for Bittensor subnet 107 (Minos). Each 72-minute round we rebuild the round's
truth (the **answer set**), search thousands of GATK configs against it, and keep
the best one submitted. The last submission before the deadline is the one scored.

## Layout and ownership

| Package | Owner | Job | KPI (measured in `backtest/`) |
|---|---|---|---|
| `common/` | everyone | Shared contracts (`types.py`), the upstream v2 formula (`scoring.py`), parameter ranges, platform client | — |
| `answer_set/` | A | `build(round)` at round start, `update(answer_set, callsets)` during the round | oracle gap, germline exactness, synthetic precision/recall/calibration |
| `evaluator/` | B | `run(round, config) -> CallSet` (expensive, fleet), `score(callset, answer_set) -> ScoreResult` (milliseconds) | s/eval, exactness vs hap.py (<1e-9), evals/s/core |
| `main/` | C | Round loop, search (`Searcher`), submission | backtested rank / P(win) against the real leaderboards |
| `backtest/` | everyone | Released rounds + real field scores, KPI functions | — |

## Rules

1. **`main` imports only the public API** of `answer_set` and `evaluator` (their `__init__.py`).
   Nobody imports another part's internals.
2. **`common/types.py` is a contract.** Changes go through a PR reviewed by all three owners.
3. **Never reimplement the scoring formula.** `common/scoring.py` wraps the upstream
   `minos_subnet/utils/scoring.py`; evaluators produce counts, the formula comes from there.
4. **Scores are tied to an answer-set version.** `run` and `score` are separate so every
   callset can be re-scored in milliseconds when the answer set improves mid-round.
5. **Answer-set updates use independent evidence** (reads, other callers, priors), never a
   vote of the configs being ranked - that circularity is what the previous project lost to.
6. **Noise floor.** Identical full search runs differed by ~0.08 points in the previous project.
   Do not claim an effect below ~0.1 from a single pair; replicate.
7. Data never goes in git: it lives under `DATA_DIR` (default `/minos/data`).
8. **All experiments run in `test/`** (git-ignored): exploratory scripts, backtest runs, scratch output.
   `tests/` holds only the unit tests that guard the contracts; those stay in git.

## Mocks, so the parts can be built in parallel

- `answer_set.OracleProvider` - released truth as an answer set (upper bound; lets `main`
  and `evaluator` work before the real builder exists).
- `evaluator.FakeEvaluator` - deterministic toy landscape through the real formula
  (lets `main` build and test the search loop without GATK).
- `main.DryRunSubmitter` - records submissions instead of sending them.

## Setup

Python 3.10+, standard library only for now. Requires the upstream repo:

    git clone https://github.com/minos-protocol/minos_subnet.git /minos/minos_subnet

Settings are environment variables (see `common/settings.py`): `DATA_DIR`,
`MINOS_SUBNET_DIR`, `PLATFORM_URL`, `FILES_URL`.

    python3 -m unittest discover -s tests -t .      # run tests
    python3 -m backtest.data --last 20              # fetch 20 released rounds (~80-170 MB each)

## What we know about the rounds (2026-09-28)

- Truth = germline of one donor + a fresh ~5K `SYNTHETIC` private pool (~25-35% of truth).
  Score core is computed over **all** of it; `N` = entries in `mutations.vcf` (434-779).
- Germline regimes: **GIAB** (NIST v4.2.1 of HG001-HG007, matched 1.000) until ~09-26;
  **1000 Genomes** 30x sample since (matched 0.83-0.93 against the phased panel;
  exact source still being pinned).
- Details and evidence: see `docs/` (to be written) and the investigation notes.
