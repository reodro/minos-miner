# search

Decides which configs to evaluate. A strategy is an `async` function `strategy(client)` that **proposes** configs
through `EvalClient` (the evaluator evaluates them) and reads their scores; several run at once. Only main
**submits**, to the platform - strategies never do.

Public API: `strategies`, `STRATEGY_NAMES`, `EvalClient`, `Outcome`, `Strategy`, `RandomSampler`, `random_stream`,
`rf_ucb`, `tpe`, `gradient_regions`.

## The contract with main and the evaluator

Main asks `search.strategies(name, ranges, workers)` for the round's lineup (`search/lineup.py`; names in
`STRATEGY_NAMES`, the first is the default) and runs them with one shared `EvalClient`. Adding, retuning or combining
strategies is a change in `lineup.py` only - never in main.

```python
fut = client.propose(config)                         # one config; await fut -> Outcome(config, result, error)
futs = client.propose_many(family)                   # several configs in ONE executor call -> can share one GATK run
outs = await client.propose_all(configs, quorum=0.9) # barrier (GA generation); stragglers keep running and report
client.latest[(key, FULL)]                           # every scored result, from all strategies (shared knowledge)
client.variant_cost(main, variant)                   # extra cost of `variant` in `main`'s batch (0..1 of a full run)
client.idle_slots()                                  # free workers right now
client.exhausted                                     # budget used up (full-run units; batch riders are cheap)
client.time_left()                                   # seconds until main stops the search (None: no deadline)
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

**Choose configs by score, never for the evaluator's convenience.** `space.PINNED` holds three values out of the
search. `emit_ref_confidence` is `NONE` because GVCF / BP_RESOLUTION output is rejected by hap.py on the validator
too. `adaptive_pruning_initial_error_rate` and `pruning_lod_threshold` stay at the GATK defaults because
HaplotypeCaller does not read them unless `--adaptive-pruning` is on, and that switch is not a search parameter.

## What exists (2026-10-01)

- `space.py`: config <-> unit-cube features (log scale for wide ranges), valid random sampling, local moves, and
  GATK cost tiers (expensive / PairHMM / genotyping). `initial_design` starts at the GATK default, takes one step
  off it on each expensive axis, then places two starts outside that neighborhood. Cheap parameters stay at defaults.
- `base.py`: `RandomSampler`, `random_stream` - the baseline every real strategy must beat.
- `rf_ucb.py`: **tiered RF-UCB**. Two ExtraTrees, on `core` and on `fp_per_target`. Candidates are ranked by the v2
  formula of those predictions (mean + kappa x spread of the trees) divided by variant cost. Kappa decays toward
  exploitation. Each pick is proposed with a designed genotyping / PairHMM family. After ~100 scores, parameters
  the trees barely use freeze at the incumbent. In the last three minutes only genotyping variants of the incumbent
  are proposed. Previous rounds' best expensive keys are seeded with cheap parameters back at defaults; their
  scores are not reused. Gate failures are labelled like a crashed config.
  Pending configs are capped at max(256, (family + 1) x workers), so a big fleet is never starved.
  v1 reference (not this loop): fake landscape 200 evals, mean best 93.2 vs random 91.7; replay 092900Z, 20 min,
  16 cores, true 87.22, rank 44/60.
- `tpe.py`: Optuna multivariate TPE beside the forest. Each suggestion uses the same designed family. The default
  lineup (`strategies("rf_ucb", ...)`, in `lineup.py`) starts both.
- `regions.py`: several trust regions. Inside a region the next step is one or two parameters, ranked by a
  forest on core and false positives, not by a finite-difference slope. A gate failure is discarded. The radius
  shrinks only after a run of local scores with no improvement. Restarts must land outside every live radius.
  Select it with `--strategy regions`. It is not the default lineup.

## How to work on a strategy

```bash
.venv/bin/python -m unittest tests.test_parts            # strategy tests on the fake landscape (seconds)
.venv/bin/python -m main replay <tag> --minutes 20       # the real thing on a stored round, true score + rank
```

`SimExecutor` + `FakeRunner`/`FakeScorer` give thousands of evals per second for developing scheduling and search
logic; `replay` is the ground truth. Compare strategies on several rounds and replicate - one replay pair is noise.

## Open work

- Tune kappa, family size and the freeze cutoff on replicated replays. v1 was not tuned on real rounds either.
- SMAC3 only if tiered RF-UCB plateaus and `pyrfr` builds. scikit-optimize is unmaintained.
- Scale: at tens of thousands of results, refit on a subsample (recent + best) or less often.

## Ideas (not planned)

- Algorithms reported for winners (only **RF_UCB** confirmed): QGA, Policy Ladder, QLock, restarts,
  InjectDrain / GA crossover. Use standard libraries where they exist (optuna for TPE/CMA-ES, pymoo or DEAP for GA).
- Legacy observations (unverified): coordinate sweeps were ~7x more productive than waves; the best configs sometimes
  need >= 2 parameters moved together (single-axis descent misses them).
