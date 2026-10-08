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

## The engine (v9)

Design and the evidence behind it: [DESIGN.md](DESIGN.md). Three steps, one after the other:

| module | what it does |
|---|---|
| `engine.py` | runs A -> B -> C on a 72-minute schedule scaled to the round, then climbs again from the best config off every hill climbed so far until the pick (the highest local score, `far.pick_s` before the end); publishes `client.notes["engine"]` |
| `candidates.py` | Step A: wave 0 (defaults, the official seed, up to 40% portfolio seeds, a Latin hypercube; 4 calling-threshold riders each), then 3 genetic generations bred from three groups of the scored configs (both: the Pareto front on core and false positives; core: highest core; fp: fewest false positives within 0.01 of the best core), 30% of children with 2 settings redrawn; then one candidate per 500 workers, far apart by callset distance |
| `climb.py` | Step B: each candidate climbs to its peak. Every combination of 5 settings at -d, 0, +d (odd iterations the read-pass settings, even ones 5 of the 7 assembly/activity settings in turn); the cheap settings ride on every point at +-d, and in pairs on the centre's run; the best point is the new centre; a winning setting's step grows in the next grid (2d, 3d...) while it keeps winning; 20 random jumps at a stall |
| `far.py` | Step C: from every peak, one setting at a time, the rungs step growth jumped over and the rungs beyond the nearest tried value; a better point, or a second high point that wins its short climb, moves the peak and Step B climbs on |
| `steps.py` | each setting's ladder: its values at the agreed resolution (d = one rung) |
| `distance.py` | callset distance (class-weighted 1 - Jaccard from the sketches on every CallSet); setting distance when a sketch is missing |
| `portfolio.py` | the opening's seeds: the legacy portfolio and our reviewed runs, the round's read-length profile first, never the replayed round's own |
| `space.py`, `bamgroups.py` | the round's search space: published ranges, pinned settings, floors, BAM groups and caps |
| `axes.py` | the Latin hypercube |
| `client.py`, `base.py` | the EvalClient the search talks to; a random baseline strategy for tests |
| `v422/` | the previous engine, a frozen backup (not selectable) |

## Engine versions

How v4 got here, version by version, is [v422/DEVELOPMENT.md](v422/DEVELOPMENT.md). The table below is the measured
record up to v4.8; v9 has no fleet run yet.

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
| v4.1 | e2e0df4, round 20261006T010500Z | scan first, valley-tested finds, evidence merges. Field best 83.876 of 67 | 25K | 84.138 | 83.554 | 37 |
| v4.2 | efb9a05, same round | steady opening, no evidence merges, find-merge ahead of the scan, plateau pick. Best reviewed TRUE 83.625 | 34K | 83.846 | 83.583 | 34 |
| v4.3 | (this tree) | raw-best pick, scan held until it finishes, one find-merge behind the scan, no new domain after the opening, scout sampler off. Local only: 14 workers, one domain, TRUE 84.832 on 20261007T044100Z, rank 56 of 70 | 633 | 83.636 | 84.832 | 56 |
| v4.4 | (this tree) | the v4.3 search, plus a ledger of each basin's opening score, first scan, and score at the 50-minute merge. Fleet replay of 20261007T055300Z is in progress on this code | | | | |
| v4.5 | (this tree) | the opening founds up to 12 basins. A basin whose first scan gains under 0.05 is retired once and replaced. Workers are split evenly until those first scans are back. Not fleet-tested | | | | |
| v4.6 | (this tree) | no merge at minute 34 or 50. All 12 opening centres climb at once. Fleet replay of 20261007T070500Z is in progress on this code | | | | |
| v4.7 | (this tree) | six of the twelve centres climb. The other six are held and swapped in, one at a time, for a basin whose first scan gained under 0.05. Not on the fleet | | | | |
| v4.8 | (this tree) | a basin whose best has not moved for 10 minutes is retired and the next candidate takes its workers. Not on the fleet | | | | |

What each step taught: riders inflate counts without information (v1); best-of-sweep beats first accept and stale
centres (v2); random shakes on spare workers win nothing after the first 20 minutes (v3); fine steps pay late in the
round (v3.1); scout finds all merged back into one hill - this round has one plateau entered from many sides (v3.2);
waited-for moves packed into big batches ran one after another on one worker (found while building v4, hence
the small/big batch setting). The v1-v3.2 code is in git history (last at 155b972); v4.22 is in `v422/`.

## Settings

Every tunable number of the engine (and the evaluator's batch limits and main's resubmission rules) is a setting in
`common/tuning.py`, read with `tune("group.name")`: `candidates` (Step A), `climb` (Step B) and `far` (Step C and
the pick). Defaults are in that file; `tuning.json` on main overrides them, and the dashboard's Settings tab edits
it. A live setting takes effect at its next use; a round setting waits for the next round. Tests and simulations
never read the file.

The engine publishes `client.notes["engine"]`: version, phase, the waves of Step A, one basin row per climber
(opened, first, first_key, best, best_key, status, iterations, waves) and the last decisions.

## How to work on the search

```bash
.venv/bin/python -m unittest tests.test_engine           # engine tests on the fake landscape (about 100 s)
.venv/bin/python -m main replay <tag> --minutes 60       # the real thing on a stored round, true score + rank
```

`SimExecutor` + `FakeRunner`/`FakeScorer` give thousands of evals per second for the controller's logic; its virtual
clock moves only when results arrive, so timing checks there are as coarse as one run. `replay` is the ground truth;
one replay or fleet run is noise at the 0.04-0.08 TRUE level, so judge changes on several rounds.
