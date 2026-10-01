# search

Decides which configs to evaluate. A strategy is an `async` function `strategy(client)` that **proposes** configs
through `EvalClient` (the evaluator evaluates them) and reads their scores; several run at once. Only main
**submits**, to the platform - strategies never do.

Public API: `EvalClient`, `Outcome`, `Strategy`, `RandomSampler`, `random_stream`, `rf_ucb`.

## The contract with main and the evaluator

```python
fut = client.propose(config)                         # one config; await fut -> Outcome(config, result, error)
futs = client.propose_many(family)                   # several configs in ONE executor call -> can share one GATK run
outs = await client.propose_all(configs, quorum=0.9) # barrier (GA generation); stragglers keep running and report
client.latest[(key, FULL)]                           # every scored result, from all strategies (shared knowledge)
client.variant_cost(main, variant)                   # extra cost of `variant` in `main`'s batch (0..1 of a full run)
client.idle_slots()                                  # free workers right now
client.exhausted                                     # budget used up (full-run units; batch riders are cheap)
```

Failed configs come back as `Outcome(ok=False)`; they never raise.

## Cost model the search should know

A **family** = configs identical except cheap parameters. Proposed together with `propose_many`, they run as one
GATK job:

| Extra config differs in | Costs (fraction of a full run) |
|---|---|
| genotyping params | ~0.01 |
| PairHMM params | ~0.12 |
| assembly / activity params | ~0.55-0.59 |
| anything else | 1.0 (its own run) |

So: explore genotyping/PairHMM settings densely around good configs, move assembly/activity/other settings more
sparingly. On about half of current rounds, variants diverge and need a ~34 s repair; still far cheaper than solo.

**Choose configs by score, never for the evaluator's convenience.** The only exclusions are settings that can never
score: `space.PINNED` pins `emit_ref_confidence` to `NONE` (GVCF / BP_RESOLUTION output is rejected by hap.py on
the validator too).

## What exists (2026-09-30)

- `space.py`: config <-> unit-cube features (log scale for wide ranges), valid random sampling, local moves.
- `base.py`: `RandomSampler`, `random_stream` - the baseline every real strategy must beat.
- `rf_ucb.py`: **RF_UCB v1** - scikit-learn ExtraTrees surrogate, upper-confidence picks (mean + kappa x spread of
  the trees), mostly local moves around the top configs plus some random ones, and a small family of cheap variants
  with every pick. Proposes whenever a worker is free; the configs it keeps pending are capped at
  max(256, (family + 1) x the executor's total workers), so a big fleet is never starved (a family fills one worker).
  - Fake landscape (6 seeds): 200 evals - mean best 93.2 vs random 91.7, RF_UCB better on 6/6 seeds; at 120 evals
    it is still warming up (91.8 vs 91.7).
  - Real replay, 092900Z, 20 min on 16 cores: 296 configs, true score 87.22, rank 44 of 60 (field best 88.96).
  - **Not tuned on real rounds**: `kappa`, family size, candidate mix, refit frequency are all first guesses.

## How to work on a strategy

```bash
.venv/bin/python -m unittest tests.test_parts            # strategy tests on the fake landscape (seconds)
.venv/bin/python -m main replay <tag> --minutes 20       # the real thing on a stored round, true score + rank
```

`SimExecutor` + `FakeRunner`/`FakeScorer` give thousands of evals per second for developing scheduling and search
logic; `replay` is the ground truth. Compare strategies on several rounds and replicate - one replay pair is noise.

## Open work

- **Move RF_UCB onto a standard library.** The surrogate is scikit-learn, but the loop around it (candidate
  generation, UCB scoring, refit timing) is hand-made. Candidates: SMAC3 (random-forest Bayesian optimization,
  ask/tell, mixed spaces; check that `pyrfr` builds), scikit-optimize `forest_minimize` (simple, barely
  maintained), Optuna (TPE/CMA-ES, the base for other strategies). Keep our part: expanding each asked config into
  a family of cheap variants for `propose_many`, and telling every result back.
- Tune RF_UCB on real rounds (replays), then add the other strategies and run them concurrently.
- Scale: at tens of thousands of results, refit on a subsample (recent + best) or less often.
- Warm start: seed each round with good configs from past rounds (a portfolio).

## Ideas (not planned)

- Algorithms reported for winners (only **RF_UCB** confirmed): QGA, Policy Ladder, QLock, restarts,
  InjectDrain / GA crossover. Use standard libraries where they exist (optuna for TPE/CMA-ES, pymoo or DEAP for GA).
- Legacy observations (unverified): coordinate sweeps were ~7x more productive than waves; the best configs sometimes
  need >= 2 parameters moved together (single-axis descent misses them).
