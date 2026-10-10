# Search engine v10: RF_UCB (2026-10-09)

The kings use RF_UCB. A month ago its sub-algorithms were QGA, QLock, Restart, Reinitialize, InjectDrain and
GA-crossover. v10 is that design, rebuilt on what the legacy project measured when it ran the same list (its
"canonical seven", `/minos/legacy/code/dev/optimizer.py`, deleted on 09-13 at median rank ~21) and on what v9 taught us.
v9 is kept in `search/v9/`; `MINOS_ENGINE=v9` selects it, and v10 is the default.

## Why the legacy RF_UCB plateaued, and what v10 does instead

| Legacy (measured) | v10 |
|---|---|
| QGA had 8 bins per setting: steps of 5-23 units. GA mutated by 10% of the range. Draws moved 12-17 settings. The winning moves were 1-4 units of one setting, "invisible to the seven proposers". | Every proposer works on the settings' ladders (search/v9/steps.py): QGA has up to 32 bins (one per rung on a short ladder), GA mutates by 1-3 rungs, and a Local proposer sends single-setting lines and small moves. |
| kappa 1.0 against a tree spread of 0.058, with moves worth 0.01: "near the top, proposals are ranked by model noise". | kappa is rescaled at every ranking: among the 5% best predicted, the bonus is `rf.kappa` x the spread of their predicted means. |
| (Not measured by legacy, the main risk named before v10's first test, 2026-10-09.) A forest does not know where it has not looked: far from the data each tree carries over its nearest box's value, the trees agree, and the spread stays low exactly where an unseen hill could be. At the start the portfolio's boxes fill the empty space, so picks land near the portfolio (the user's concern). | UCB has a third term: distance to the nearest evaluated config (mean ladder distance per setting), scaled like kappa (`rf.kappa_dist`). A candidate in an empty region is favoured because nothing is known there; the bonus fades as the region fills. |
| The sweep lane (coordinate descent from the model's best, on its own hosts, from mid-round) produced every winning config after the opening; the model produced the heads. | 10% of the picks are the sweep: every setting's line through the actual best, unranked, rebuilt at every refit, from the start of the loop. |
| Each model pick was one GATK run. | Each pick carries 3 cheap riders (`rf.riders`: genotyping or PairHMM settings one rung away) in the same run: the forest gets about twice the results per run. |
| Restart after half a wave (21 restarts in 65 minutes), population cut from 192 to 8. | Restart after 2.5 fleets' worth of results without a 0.005 gain. The forest keeps every result. |
| QLock could freeze an axis inside the first wave (12 observations), at a resolution of 0.03. | At least 60 results, resolution 0.003, tree spread included. |
| From a plateau, model + GA gained nothing in 2 hours; sweeps of a QGA-found head gained +0.17. "The model produces the diverse heads; what was missing is sweeping them." | Local proposals (lines, small moves) around the 4 best distinct configs join every pool, and the last 15 minutes are v9's grid climb with lines on the best two distinct configs. |
| 40% of the opening from the portfolio. | 15% (user's concern, 2026-10-09: a forest fitted on clusters of other rounds' configs predicts high near them). Configs in crowded regions also share their weight in the fit. |
| Waves: "a wave ends with its slowest job; lane hosts measured 21-27% busy". | No waves: every free worker takes the next ranked config at once. |
| Two-head forest (core, FP): built, no measurable effect. | One forest of the score. |
| Refit every workers/17 results, synchronously, full history, 8 threads; "the controller must stay far cheaper than the evaluations it feeds". | The same cadence, but self-throttled and in the background: a refit starts when the previous one has finished and workers/17 new results are in, on every result in by then; the loop never waits for it. 14 threads (main has 16 cores and uses few); each tree grows on a random 30% of ALL configs (no region preferred): 3.4 s at 40,000 configs instead of 14 s. Encodings cached. The controller's own time per result is logged (`loop_ms_per_result`). |
| The top of the UCB ranking evaluated as it comes. | One merged pool ranked by UCB, as the original. But a forest predicts the same score across a leaf's box, so the top of the ranking is often near-duplicates: 20% of the picks are mutations (1-3 settings, 1-3 rungs) of the 50 best-ranked proposals, not asked of the forest. They spread the batch in and around the best box and hedge the forest's error where it is most confident. QGA also gets 10% of the picks: in one ranking its many-setting draws never reached the top (0 picks in two simulated rounds), and they are one of the few ways v10 reaches other hills. (Tried and dropped 2026-10-09: several candidates with equal shares, which split RF_UCB into small ones and overrode the forest's global judgement.) |

| (Measured on v10's first fleet replay, 20261008T175300Z, 2026-10-09.) By minute 10 the loop worked both ends of the core/false-positive trade-off; at minute 16.8 its high-core best (core 0.8712, 1,653 FP) was 0.003 of core from the kings (0.8745, 1,655). But every climbing source follows the overall best, which sat ~0.1 higher on the low-FP side, so the high-core share fell from 42% to 3% and that end stalled unclimbed. | The best score and both ends of the trade-off (highest core and fewest false positives within 0.5 of the best) are always among the centres, and the endgame climbs all three. Not "follow the high-core end": each round's winners may sit at either end. |

## Fixes from the first fleet replay (20261008T175300Z, 2026-10-09: TRUE 76.436, rank 11 of 72; v9 76.106, rank 40)

| Seen | Fix |
|---|---|
| The loop started with no ranked pool: about 770 workers freed by the opening took random Inject configs. | The first refit runs before any pick; the loop starts at 50% of the opening back; with nothing ranked, Inject keeps its 10% and the rest waits. |
| After the restart (minute 35) QGA was uniform and its far draws won the novelty bonus: 52% of picks, nothing to climb. | The novelty bonus is capped at 3 typical bonuses. |
| A Reinitialize at minute 49 sent 972 random configs into the loop's last 5 minutes. | No Restart or Reinitialize within 10 minutes (scaled) of the endgame. |
| The endgame's first grid wave took 15 minutes behind the loop's leftovers. | The endgame's proposals carry priority 10; every executor's queue orders by priority. |
| The forest never saw 6,246 final failures or 7,790 zero scores, so it kept proposing next to crash zones. | Final failures (after the evaluator's solo and cold retries) and zero scores are learned by a separate classifier (see the second replay's fixes). |
| The high-core end (0.003 of core from the kings) was found, then dropped. | The best score, best core and best FP are always centres, and the endgame climbs all three (see above). |

## Fixes from the second fleet replay (same round, 2026-10-09: TRUE 76.447, rank 11 again, 0.057 from the kings)

Both runs ended on the same hill (core 0.864, ~1,554 FP); the kings sit at core 0.8745, 1,655 FP.

| Seen | Fix |
|---|---|
| The best-core config (44086fc8: core 0.8745 like the kings, 1,738 FP, 75.997) was not climbed: the centres' callset-distance rule dropped it as "near" the best score. | The three ends (best score, best core, best FP) always take their slots unless they are the same config; only the extra centres need the distance. |
| The loop converged by minute ~30; after that 82% of picks polished the summit for +0.009. | When the best gains less than `rf.converge_gain` (0.02) for `rf.converge_min` (10) minutes, after `rf.converge_after_min` (20; both scaled), the loop ends and the endgame starts at once (on the second replay: about minute 40 instead of 55). Counted in minutes, not results: riders make a fleet's worth of results about 40 s, and a first try counted that way fired at minute 23 of a simulated round whose earlier runs still gained about 5 after it. Simulated rounds gain 0.02 at most 12 minutes apart, so they lose little. A restart starts a new window (175300Z, third run: a restart at minute 40 was cut off by convergence 1.2 minutes later). `rf.endgame_min` stays the latest start. |
| Teaching failures at a floor pulled every prediction down by 1.4-2.1 and lowered the sweep's Spearman from 0.91 to 0.76. | The score forest is fitted on scored configs only. A separate classifier (forest.fit_bad, from 20 bad configs) gives P(bad); UCB subtracts P(bad) x the pool's UCB range, so a likely crash sinks and a safe config is unchanged. |
| Both runs used seed 0: the same 972 opening configs, so not independent tests. | A seed per run (time-based, logged in the first engine line); `MINOS_SEED=<n>` repeats one. |
| The endgame's waves still queued behind the loop's leftovers, and each grid point carried 16 riders. | At the endgame the loop's picks that no worker has started are dropped (`cancel_queued` on every executor: they leave in_flight at once; running ones finish). Endgame grid points carry `rf.endgame_riders` (4) riders. |

Crashed batches (evaluator/gatk.py, BATCH_CRASH_NO_COLD): a batch (a head with riders) that dies on a server with
"Graph must have ref source and sink vertices" is no longer re-run cold. The plain batch path dies the same way on such
heads every time, and the executor re-runs the head alone next, which keeps its cold re-run, so every config still gets
a stock-equal result. A solo job that crashes is still re-run cold: the evaluator's own notes say the batch path crashes
on configs that pass on their own, so the cold solo run is the guard. The real fix, the patch's batch path not crashing
where stock GATK does not, belongs to the evaluator developer. MINOS_COLD_AFTER_BATCH_CRASH=1 restores the old re-run.

## The round

| Minutes (72-minute round, scaled) | What runs |
|---|---|
| 0 to ~8 | **Opening**: the defaults, the official seed, 15% portfolio, the rest a Latin hypercube over every setting; 4 calling-threshold riders each. The loop starts when 80% is back. |
| ~8 to 55 | **RF_UCB loop** (asynchronous). Background refits as above; each rebuilds one pool (QGA 55% and GA 35% of 6,000 proposals, plus Local's lines and small moves around the 4 best distinct configs) and ranks it by UCB. Each free worker takes: 10% a fresh hypercube config, 20% a mutation of a top proposal, 10% QGA's best remaining proposal, 10% the sweep (every setting's line through the actual best, unranked), 50% the next of the ranking. QLock, Restart and Reinitialize as above. |
| converged (from 20) or 55, to the pick | **Endgame**: the loop's queued picks are dropped; v9's grid climb with lines (4 riders a point) on the best score, the best core and the best FP, then v9's Step C around each peak, then "climbing again" from the next best. |
| 70.5 | **Pick**: the highest local score. |

## Parts

- `forest.py`: the forest (sklearn, 200 trees, min leaf 2), ladder encoding, crowding weights, UCB, QLock's flatness test.
- `proposers.py`: QGA, GA-crossover, InjectDrain's inflow, Local, the core/FP front, distinct configs.
- `engine.py`: the opening, the loop (quotas, refits, QLock, Restart, Reinitialize), the endgame, the pick, the notes
  (`client.notes["engine"]["rf"]`: refits, restarts, reinits, frozen settings, picks and new bests by source).
- Settings: the `rf` group in `common/tuning.py`.

## What the test logs

Every pick's prediction (mean, tree spread, novelty) next to its actual score: `<status_dir>/rf_picks/<round>_<time>.jsonl`,
and per source in `client.notes["engine"]["rf"]["calibration"]` (picks back, bias, Spearman of predicted vs actual,
median novelty). It answers whether the forest predicted well and whether UCB sent workers into new regions.

## The bar, agreed before the test

v10 succeeds on `20261008T175300Z` only if it beats v9's 76.106 TRUE by more than run noise (about +0.08), or reaches
the kings' point (core 0.874 at about 1,650 false positives). Anything else is a failure.

## Not known yet

Whether it finds the kings' hill. The test: fleet replays with the true answer set on `20261008T175300Z` (v9: TRUE
76.106, rank 40 of 72; the kings sit at core 0.874 with 1,650 false positives) and 2-3 other rounds, v10 against v9 on
the same rounds.
