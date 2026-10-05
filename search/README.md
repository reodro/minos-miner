# search

Decides which configs to evaluate. Main asks `search.searcher(ranges, workers)` for the round's search (one `async`
function of the round's `EvalClient`), runs it until the deadline and submits the local best, or the config the search
names in `client.pick_key`; the search only **proposes** configs and reads their scores.

Public API: `searcher`, `EvalClient`, `Outcome`, `Strategy`.

## The contract with main and the evaluator

```python
futs = client.propose_many(configs, batch="small")   # a result the search waits for: one per worker, cheap riders only
futs = client.propose_many(family)                    # batch="big" (default): pack with configs sharing a GATK run
futs = client.propose_many(group, batches=[...])      # one value per config: a small head with "big" riders
client.latest[(key, FULL)]                            # every scored result (shared knowledge)
client.idle_slots(); client.capacity()                # free and total workers
client.time_left()                                    # seconds until main stops the search (None: no deadline)
client.pick_key = key                                 # main submits this config instead of the best
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

Packing saves compute but a batch runs on one worker, its configs one after another: 30 assembly variants in one
batch save ~45% of the work but the last result arrives after ~16 run-times. So a result the search waits for goes
out as a small batch (itself plus its cheap siblings), and only work nobody waits on is packed big.

**Choose configs by score, never for the evaluator's convenience.** `space.PINNED` holds three values out of the
search. `emit_ref_confidence` is `NONE` because GVCF / BP_RESOLUTION output is rejected by hap.py on the validator
too. `adaptive_pruning_initial_error_rate` and `pruning_lod_threshold` stay at the GATK defaults because
HaplotypeCaller does not read them unless `--adaptive-pruning` is on, and that switch is not a search parameter.

## The engine (v4)

| module | what it does |
|---|---|
| `hills.py` | the one rule for which hill a config is on (4+ settings more than 10% of range apart = another hill); ownership, adoption, the scout's exclusion zone, merging and founding all ask it |
| `genetic.py` | recombination: niched genetic generations, and the merge operator - every combination of the settings where two centres differ (when 8 or fewer, else 256 random ones) plus mutated children |
| `climber.py` | a domain's climb, iterated local search (ParamILS): coarse then fine sweeps, best of each sweep with the other winners re-applied, sideways moves on plateaus, the perturb-and-climb loop (shake 2/3/5/8 settings, a local round, keep better OR equal); it reports facts and never decides its own lifecycle |
| `engine.py` | the opening and the controller: everything about domains is decided in one pass from facts |
| `axes.py` | per-setting steps, learned cliffs, parked settings, the distance between configs, the Latin hypercube |
| `families.py` | genotyping riders for global points |
| `ranker.py` | the (core, FP) forest the scout ranks candidates with |
| `space.py` | the round's search space: published ranges, pinned settings, floors, BAM caps |
| `client.py`, `base.py` | the EvalClient the search talks to; a random baseline strategy for tests |

The round (times for a 72-minute search; they scale with the round's length):

| minutes | stage |
|---|---|
| 0-12 | wave 0 (seeds + a Latin hypercube over every setting, 4 threshold riders) while the answer set is built, then 3 genetic generations, each launched once 60% of the previous is back or workers idle; stops early when a generation gains < 0.05; hard limit 14 |
| 12 | founding: one domain per hill (registry), up to 6 |
| 12-62 | ParamILS in every domain; winning settings shared between domains every 5 min; the scout keeps 15% of the workers (hypercube points, crossovers between hills, model upper confidence bound); a scout find founds a domain only after a midpoint run shows a valley (until minute 34) |
| any time | merge on evidence: two domains on one hill, or a frozen domain on a live one's hill |
| 34-40 | scheduled merges, closest pairs first: at most 3 domains from minute 40 |
| 56-62 | a three-way merge: one domain |
| 62-72 | the final on one point (halving single steps, all pairs, a last pairing with the runner-up), then the plateau-centre pick: among configs within 0.05 of the best, the one whose neighbours score highest |

Every merge evaluates its children in parallel (small batches: one config per worker, only cheap siblings riding) and
takes the best of parents and children as the new centre, which inherits both parents' tested configs and steps.
Riders: 4 on wave 0, 2 on generations and scout points; a climbing move carries 2 siblings (one genotyping, one
PairHMM step); merges and the final carry none.

## Engine versions

All on round 2026-10-02 09:29 (field best 89.912, 66 rows), 72-minute fleet runs unless noted. TRUE = scored on the
released truth; one run's noise is 0.04-0.08 TRUE.

| version | commit | what it was | configs | local | TRUE | rank |
|---|---|---|---|---|---|---|
| legacy pipeline v2 | (legacy project) | the friend's waves, crossover, trial of 6 heads keep 3; single machine | 87 | | 89.16 | 38 |
| v1 | 2026-10-04 | basin climbing, 24 riders per head, first accept | 245K | 88.701 | 89.645 | 34 |
| v2 | 2026-10-04 | basins with best-of-sweep acceptance, valley tests, scouts adopted | 103K | 88.728 | 89.736 | 24 |
| v3 | 2026-10-05 | isolated domains, ParamILS-style climb, funnel 6-3-1, shakes streamed | 32K | 88.661 | 89.612 | 34 |
| v3.1 | 1602b81 | fine mode, freezes, merges with inheritance, siblings | 59K | 88.778 | 89.714 | 25 |
| v3.2 | 155b972 | hypercube waves, model-guided scout, rolling domains | 116K | 88.771 | 89.770 | 23 |
| v4 | 2026-10-05 | hill registry, controller, genetic opening, merges by recombination, ParamILS perturbation loop, plateau-centre pick (hand-tuned during the run: opening.max_min, climb.decisive) | 47K | 88.699 | 89.711 | 25 |

What each step taught: riders inflate counts without information (v1); best-of-sweep beats first accept and stale
centres (v2); random shakes on spare workers win nothing after the first 20 minutes (v3); fine steps pay late in the
round (v3.1); scout finds all merged back into one hill - this round has one plateau entered from many sides (v3.2);
waited-for moves packed into big batches ran one after another on one worker (found while building v4, hence
the small/big batch setting). The v1-v3.2 code is in git history (last at 155b972); this package holds only v4.

## Settings

Every tunable number of the engine (and the evaluator's batch limits and main's resubmission rules) is a setting in
`common/tuning.py`, read with `tune("group.name")`: the opening, domains, schedule, merges, sharing, scout, riders,
climb and pick. Defaults are in that file; `tuning.json` on main overrides them, and the dashboard's Settings tab
edits it. A live setting takes effect at its next use, within about a second; a round setting (`domains.max`,
`riders.wave0`, main's) waits for the next round. Tests and simulations never read the file. Each run's manifest
records the settings at its start and every change made during it, so a hand-tuned run can be told apart later.

The engine publishes what the dashboard draws in `client.notes["engine"]` every 2 s: phase, scaled schedule, waves,
every domain with its origin and parents, merges, scout, pick and the last 200 decisions.

## How to work on the search

```bash
.venv/bin/python -m unittest tests.test_engine           # engine tests on the fake landscape (about a minute)
.venv/bin/python -m main replay <tag> --minutes 60       # the real thing on a stored round, true score + rank
```

`SimExecutor` + `FakeRunner`/`FakeScorer` give thousands of evals per second for the controller's logic; its virtual
clock moves only when results arrive, so timing checks there are as coarse as one run. `replay` is the ground truth;
one replay or fleet run is noise at the 0.04-0.08 TRUE level, so judge changes on several rounds.
