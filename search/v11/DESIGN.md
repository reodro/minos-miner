# Search engine v11: RF_UCB over several hills (2026-10-10)

v11 = v10 (forest, proposers, riders, crash classifier, endgame climber, seed per run) + explicit hills.
Select it with `MINOS_ENGINE=v11` (v10 stays the default until v11 is tested on the fleet).

## Why

v10 was meant to work several hills at once and did not:

| Measured (v10, 20261009T153000Z, final best 86.548) | |
|---|---|
| Share of the best 50 configs on the final hill | 10/50 at minute 10, 48/50 at minute 20, 50/50 from minute 30 |
| Share of all evaluations on the final hill | 2% (0-10 min), 15%, 56%, 63%, 55%, 77%, 81% (60-70 min) |
| Configs above 86.3 not on that hill | 0 of 24,176 |

Local tests around that best found nothing better (test/compound_het/README.md: 300 multi-setting moves, 656
sub-rung genotyping values, climbs of two other PCR models, basin recombination), while the field sat 0.056 higher.
The gap is a better hill, not a better polish.

Why v10 collapsed (implementation, fixed here): every proposer but Inject followed the global best (GA bred from the
global top, QGA rotated toward it, Local and Mutate around it, the sweep through it); its "centres" were the ends of one
ridge; a restart kept the forest and elites, so the next ranking found the same hill; the endgame climbed that ridge.
What RF_UCB cannot do by itself (strategy): a forest cannot rate a hill it has not sampled, and near the top its error
(~0.05) is the size of the differences, so late progress is local climbing - single-hill unless funded per hill.

## What v11 does

| Part | v11 |
|---|---|
| Hills | At every refit the best `hill.points` (2,000) results are grouped by what they call: a config joins the first hill whose best is within `hill.callset_far` (0.03) by v10's callset distance (class-weighted 1 - Jaccard of the per-class bottom-128 MinHash sketches every result carries; settings distance `hill.setting_far` only without sketches), else founds a hill if within `hill.window` (1.5) of the best; at most `hill.track` (24). A founder near an old hill's founder keeps that hill's id and gain clock. Settings that change no call (range edges, BAM equivalence classes) do not split a hill. The sketch is an estimate (exact only for classes under 128 calls): Calibrated on saved callsets (175300Z, 044200Z): variants of one hill 0.000-0.016 from the best, different trade-offs (core 0.865 vs 0.859 at FP 1,221-1,236 vs 1,115; another PCR model) 0.067-0.095; 0.05 merged the best-core end into the best's hill in v10. ~70 us a pair, cached per (founder, config). |
| Valleys | The radius is small on purpose (a hill split in two is still climbed; a real hill merged into another is lost). Duplicates are merged by a hill-valley test (HillVallEA) on what is already known: configs whose calls lie between two hills' bests (closer to each than they are to each other, at least `hill.valley_mid` (25%) of that distance from both) and the cross children taking 2-3 of the 5 stage blocks from each side. One config between them within `hill.valley_drop` (0.1, just above the near-peak roughness) of the lower: one hill, merged into the higher. `hill.valley_min` (5) or more, all lower: a valley, kept apart. Else unknown: that pair's midpoints go first in the cross share. No extra evaluations: both kinds are evaluated anyway; the cost is waiting for midpoints (~6 a minute per pair at 1,056 workers). Pairs among the best `hill.valley_hills` (6). |
| Funding | `hill.count` (4) hills by best score, except one `hill.drop` (0.5) below the leader with no gain of `rf.converge_gain` for `hill.patience_min` (10, scaled): its place goes to the next. |
| Pools | Per funded hill, ranked by the forest: GA from its own members, Local around its own trade-off ends, the sweep through its own best, mutations of its own top proposals. |
| Exploration | `rf.share_inject` (10%) hypercube; `hill.share_qga` (10%) QGA rotated toward every funded hill; `hill.share_cross` (10%) children of two funded hills' bests mixing whole stage blocks (read, activity, assembly, PairHMM, genotyping), ranked by the forest. |
| Picks to hills | The rest: the leader `hill.lead_share` (40%), the other funded hills equal shares. Within a hill: `hill.share_sweep` (20%) its sweep, `hill.share_mutate` (20%) a mutation, else its ranked pool. |
| Off | QLock (flat at one hill's best is not flat at another's) and Reinitialize (exploration is continuous); a stall resets QGA. |
| Convergence | Every funded hill quiet (no gain of `rf.converge_gain` for `rf.converge_min`, each on its own clock; a new hill starts its clock when founded), after `rf.converge_after_min`; else `rf.endgame_min`. |
| Endgame | v9's grid climb with lines and Step C on the best of each funded hill within `hill.endgame_window` (0.3) of the leader (at most `hill.endgame`, 4), then "climbing again". An endgame climb gains ~0.1-0.25 (044200Z: the FP end +0.47 to 0.004 below the best), so hills further down cannot catch up and would only slow the leader's waves. |
| Pick | The highest local score (unchanged). |

## The test bar

1. Simulator, multi-hill landscapes (`HILLS`=3), v10 and v11 on the same seeds: v11 keeps more than one hill worked
   past minute 30 and does not lose on the final score.
2. Fleet, a round with the field's numbers (20261010T044200Z: v10 with fixes is the baseline; the search developer's
   v9.18 82.638): v11 at or above v10, and how many funded hills reached within 0.1 of the best.

## Test log

| Run | Result |
|---|---|
| 044200Z, settings-distance hills (2026-10-10 09:07) | Stopped at minute 29 (best 82.520, on v10's pace). Its 4 funded "hills" were one ridge: 15 of 25 settings identical, core 0.856-0.860, FP 1,114-1,161; they differed in genotyping settings that are far on their ladders but change few calls. Hill ids passed 180: one hill re-founded again and again. Hence callset distance. |
