# search

Decides which configs to evaluate. Main asks `search.searcher(ranges, workers)` for the round's search (one `async`
function of the round's `EvalClient`), runs it until the deadline and submits the local best; the search only
**proposes** configs and reads their scores. The search is `engine.py`, asynchronous basin climbing. The older
strategy modules (`basin.py`, `rf_ucb.py`, `pattern.py`, `tpe.py`, `regions.py`) stay importable for experiments
and tests; main never runs them.

Public API: `searcher`, `EvalClient`, `Outcome`, `Strategy`.

## The contract with main and the evaluator

Main asks `search.searcher(ranges, workers)` for the round's search (today the basin race, `basin.py`) and runs it
with one shared `EvalClient`, in replays and live alike. Changing the search is a change in `search/` only - never in
main; the other strategy modules stay for experiments and tests.

```python
fut = client.propose(config)                         # one config; await fut -> Outcome(config, result, error)
futs = client.propose_many(family)                   # several configs in ONE executor call -> can share one GATK run
futs = client.propose_many(step, lane="latency")     # a result the search waits for: fastest free node, cheap riders only
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

## The engine: asynchronous basin climbing (`engine.py`)

One constraint shapes it: a round allows few sequential steps and unlimited width. A GATK traversal takes minutes
under load (median 4 min, p90 8 min on a mixed fleet), so a 70-minute search is at most ~15-20 dependent steps;
workers only widen each step. So nothing waits at a barrier, every step is as wide as the free workers allow, and
an improvement is used the moment it arrives. Every count is derived from the executor's capacity, so a 16-worker
replay on main and a 2,000-worker fleet run the same code.

| part | what it does | module |
|---|---|---|
| opening | seeds (GATK default, the official seed, the known family shape, earlier rounds' bests) and a designed wave over the dominant axes (`DOMINANT_LEVELS`: 2-3 levels each, in balanced combination); every skeleton carries a genotyping family | `axes.design_wave`, `families` |
| basins | k regions at least `RADIUS` (0.15, RMS unit distance over structural axes) apart; k from capacity (`MOVES_PER_BASIN` in flight each). A basin keeps a centre and a move pool; a move is dispatched whenever the basin has a free worker, in the **latency lane**; an accept re-centres at once and rebuilds the pool; a single-axis move that beat its own (older) centre becomes a merge candidate on the new one | `basins.Basin` |
| moves | steps (one axis, +-1x/0.5x/2x its step), merges of recent winners (acc2..accN, pairs of the top four), transplants of 1-2 axes from other basins' bests, full-range probes at birth and on a stall, model-ranked pair moves | `basins.Basin.refill` |
| axis book | per-axis evidence shared by every basin: step sizes (grow on a win, shrink on failures), axes whose moves change nothing are parked (`DEAD_AFTER`, retested after `PARK_S`), learned cliffs (a failing value bounds the axis) | `axes.AxisBook` |
| scouts | 25% of the width (15% late): crossover and mutations of the basins' bests, unrun design cells, range edges, a few random configs; ranked by the (core, FP) forest by predicted gain per cost; a scout beating the weakest basin by `BIRTH_MARGIN` births a basin if it is a region of its own; basins that converge merge; a basin `RETIRE_BEHIND` the leader after `RETIRE_PATIENCE_S` retires; a retired region is not reborn unless beaten | `engine.scout_candidates`, `ranker.Ranker` |
| finish | the last `FINISH_S`: dense genotyping and PairHMM families on the best skeletons, half steps only | `families` |
| riders | a genotyping rider costs ~1% of a traversal, so every skeleton carries up to 24 (6 on a small fleet); the skeleton's score is the best of its family; riders go in the throughput lane in the same call as their head so they pack into its batch | `families.genotyping_riders` |

The forest ranks pools larger than the free width (scouts, pair moves). It never accepts a move: on a stepped
landscape it cannot resolve the 0.005-0.05 differences that decide a step; acceptance is always an exact
evaluation (measured on the legacy project: exact sweeps 90.97 -> 91.15 in ~1,500 evaluations, a model + GA
14,805 evaluations with no gain).

Scaling (virtual clock, fake landscape): 16 workers -> 1 basin, ~12 moves in flight; 950 -> 4 basins, ~140
moves each + ~240 scouts; 1,900 -> 7 basins. Utilisation stays at ~1.0 in every case.

The test harness: `tests/test_engine.py` runs whole virtual rounds through main's loop (`SimExecutor` batches
riders like the fleet does); `test/p2/engine_sim.py <capacity> <minutes>` prints the engine's log for one.

## How to work on the search

```bash
.venv/bin/python -m unittest tests.test_parts            # strategy tests on the fake landscape (seconds)
.venv/bin/python -m main replay <tag> --minutes 20       # the real thing on a stored round, true score + rank
```

`SimExecutor` + `FakeRunner`/`FakeScorer` give thousands of evals per second for developing scheduling and search
logic; `replay` is the ground truth. Compare strategies on several rounds and replicate - one replay pair is noise.

## Open work (engine)

- Measured on fleet replays, not yet: `RADIUS`, k, `RETIRE_BEHIND`/`RETIRE_PATIENCE_S`, the scout share,
  step rules, what the forest's ranking adds over plain ordering.
- A richer virtual landscape (plateaus, cliffs, several basins) for tuning the rules before fleet time.

## Open work (older strategies)

- Tune kappa, family size and the freeze cutoff on replicated replays. v1 was not tuned on real rounds either.
- SMAC3 only if tiered RF-UCB plateaus and `pyrfr` builds. scikit-optimize is unmaintained.
- Scale: at tens of thousands of results, refit on a subsample (recent + best) or less often.

## Ideas (not planned)

- Algorithms reported for winners (only **RF_UCB** confirmed): QGA, Policy Ladder, QLock, restarts,
  InjectDrain / GA crossover. Use standard libraries where they exist (optuna for TPE/CMA-ES, pymoo or DEAP for GA).
- Legacy observations (unverified): coordinate sweeps were ~7x more productive than waves; the best configs sometimes
  need >= 2 parameters moved together (single-axis descent misses them).
