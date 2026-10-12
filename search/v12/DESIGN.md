# Search engine v12: trust regions with asynchronous Thompson sampling (2026-10-10)

Select it with `MINOS_ENGINE=v12` (v10 stays the default until v12 is tested on the fleet). Settings: the `tr.*` group.

## Why

| Measured | |
|---|---|
| v10 commits to one hill by minute ~15-20 (153000Z: 48 of its best 50 configs on the final hill at minute 20) | the gap to the king is a better hill, not polish (test/compound_het/README.md) |
| v11 grouped evaluated configs into hills; on the fleet (044200Z, 2026-10-10) its merge rule ("one config between two hills within 0.1 of the lower one: same hill") cascaded and rebuilt one hill within 5 minutes | a hand-set radius and merge rule decide hill size badly |

v12 does not group configs at all. A hill is a search region whose size is learned from its own successes and
failures, as in the published trust-region methods.

## Sources and what each contributes

| Method | Used for |
|---|---|
| TuRBO (Eriksson et al., NeurIPS 2019) | several trust regions; grow/shrink on success/failure; restart on collapse; Thompson sampling across regions as an implicit bandit; defaults L 0.8 / min 0.5^7 / max 1.6 |
| CASMOPOLITAN (Wan et al., 2021), Bounce (Papenmeier et al., NeurIPS 2023) | regions in mixed and ordinal spaces as "how many variables may change"; Bounce's budget-aware shrinking (here: a late clock cap) |
| Asynchronous parallel Thompson sampling (Kandasamy et al., AISTATS 2018) | every free worker draws one sample of the model and takes its best point; better than synchronous batches under a time limit when run times vary |
| Random forest in TuRBO (BOCoDe benchmark), SMAC | the forest as the model (a Gaussian process does not scale to 60K results); one tree = one Thompson draw |
| SCBO (constrained TuRBO) | candidates the crash model rates more likely bad than not are dropped |

## How it works

| Part | v12 |
|---|---|
| Opening | v10's: defaults, official seed, `rf.portfolio_share` portfolio, Latin hypercube; cheap riders |
| Regions | `tr.regions` (4) start at the opening's best configs at least `tr.new_far` (0.067 callset distance: the measured different-hill range) apart |
| A region | its centre (best config) and length L (`tr.length_init` 0.8). A candidate changes up to ceil(L x settings / 2) settings, each by up to L/2 x its ladder x its box weight (1 / forest importance, geometric mean 1, clipped to [1/5, 5]: settings that matter move less). `tr.candidates` (2,500) per region per refit |
| Steps | every `tr.step_results` (100) of the region's picks back, or `tr.step_max_min` (3, scaled) without that: success if its best gained `tr.success_gain` (0.005). `tr.success_tol` (3) successes double L (max 1.6); `tr.failure_tol` (2) failures halve it |
| Clock cap | L <= 1.6 x (L_min / 1.6)^(f^3), f the share of the loop elapsed: 0.82 at half the loop, 0.10 at 80%, L_min at the endgame |
| Restarts | a region below `tr.length_min` is finished (centre to a tabu list); a new one starts at the best evaluated config `tr.new_far` from every active and finished centre. Two active centres within `tr.dup_far` (0.016: one hill's variants measured 0.000-0.016): the weaker restarts. No new region within `tr.no_restart_min` (20, scaled) of the endgame: duplicates are still closed and their workers go to the other regions; collapses are not closed late (the clock cap shrinks every region then) |
| Picks | `rf.share_inject` (10%) hypercube; each region's floor `tr.min_share` (8%: its best candidates by the forest's mean); the rest asynchronous Thompson sampling: one random tree, its best candidate of all regions not taken. Each pick carries v10's cheap riders |
| Endgame | the loop's queued picks dropped; v9's grid climb with lines and Step C on the centres within `tr.endgame_window` (0.3) of the best (at most `tr.endgame`, 4), plus the overall best if it is in no region; then "climbing again" |
| Pick | the highest local score |
| Off | QGA, QLock, Restart/Reinitialize of v10, v11's grouping and valley rule, the convergence clock |

## Adaptations to our scale (not in the papers)

- Steps count returned results per region (the loop is asynchronous, 1,000 workers), not batches.
- Thompson sampling alone starved the weaker regions in a simulated round (one region had 28 results back in 50
  minutes, never completed a step, never restarted): hence the floor share and the timed steps.
- A plain geometric clock cap set every region's length (learned sizes never mattered): hence the cube.

## Test log

| Run | Result |
|---|---|
| Simulator, 3-hill landscape 0, seed 0, before the floor/timed steps/soft cap | pick 95.064; r3 28 results in 50 min; no restarts; L followed the cap |
| Same after the fixes | pick 95.499; all regions 450-2,000 results; r1 (started 91.85, below r0's 92.42) overtook r0 at minute 22 and ended best; no restarts (all regions kept improving) |
| Simulator A/B, 3-hill landscapes, 2 landscapes x 3 engine seeds, 300 workers (2026-10-10) | v12 won 5 of 6, mean +0.155: landscape 0 95.499/95.499/95.501 vs v10 95.045/95.063/95.498 (v10 stuck on a lower hill twice); landscape 1 97.378/97.512/97.512 vs 97.495/97.378/97.496. The loss came with 3 region closures at minutes 38-40 (restart churn): hence no new regions in the last 20 minutes before the endgame |

## Fixes after the first fleet test (044200Z, 2026-10-10: TRUE 82.715, rank 14; v10 82.655, rank 24)

What the test showed: the winner (82.715) came from the endgame's line scans on region r2, which never led in the loop
(82.47). Its gain needed base_quality_score_threshold 16 -> 33 (+17 rungs) and the gap penalty 28 -> 21 (-7): the 394
loop configs near r2's centre all had BQST 15-18. Restarts never started (every good config within 0.067 of a centre),
the best (82.653) sat in no region for 30 minutes, and Inject took 18,861 results (23%, best 81.40) mostly in minutes
40-60 when the two last regions were too small to have new candidates.

| Fix | |
|---|---|
| Lines per region | `tr.share_lines` (25%) of the picks walk a region's lines: every setting's untried rungs along its whole ladder through the centre (cheap settings included: v10's sweep skipped them), coarse to fine (both ends and the middle first), interleaved across settings; unranked, not limited by L, rebuilt at every refit around the current centre |
| Equal box weights | `tr.importance_weights` 0: importance is not a lengthscale; the setting that needed +17 rungs was kept on the shortest leash |
| Early endgame | in the last `tr.no_restart_min`, when no region has gained for `tr.spent_min` (6, scaled), the endgame starts (044200Z: ~minute 40 instead of 55) |
| Fresh restarts | when no evaluated config is `tr.new_far` from every centre, a new region starts at a fresh hypercube config far (by settings) from all of them, its centre evaluated first (TuRBO's restart from a fresh design) |
| Reopen the best | the best config in no live region gets a small region (`tr.reopen_length` 0.1), late too, once per config |

## Finding the right area: successive halving (2026-10-10)

The fixes above work harder on hills we hold. On 044200Z the field's top cluster (king and #5: core 0.870 at FP ~1,232)
is an area we never reached (our best core at that FP 0.8659); we sat with #10-#14 (core ~0.862, FP ~1,140). Inject
(18,861 single random configs, best 81.40) cannot find narrow hills: a config shows its hill's potential only once
climbed, and the opening's ranking does not predict summits (the winning region started third).

| Part | |
|---|---|
| Screening | after the opening, `tr.sh_regions` (16) regions open at the opening's results far apart in calls, going well below the top 4 |
| Fair share | while screening, `tr.sh_floor` (80%) of the picks go to the regions in equal turns (Thompson sampling would favour the current leader) |
| Halving | every `tr.sh_step_min` (4, scaled) the better half by best + `tr.sh_gain_weight` (1) x gain since opened is kept (Jamieson & Talwalkar 2016; Hyperband): 16 -> 8 -> 4 by about minute 12-13 |
| Then | normal v12 with the survivors (`tr.regions` 4); restarts, lines, reopen and the early endgame as above |
| Off | `tr.sh_regions` <= `tr.regions` |

## Screening for the king's hill: indel slot and far starts (2026-10-11)

Every setting the platform allows has now been covered without reaching the king's hill: our hill's neighbourhood,
the boundary-pinned settings, the three jump families (test/jump_family/README.md) and the fixed settings
(`sample_ploidy` 2 and `dont_use_soft_clipped_bases` False are the platform's only allowed values, for everyone).
So the hill is inside our space, and the screen has to start on it and keep it. Two things worked against that:
every screening start came from the top of the opening's ranking (our usual hill and its neighbours), and the
halving ranked on score alone, while the king's edge is the indel part of core at equal FP (18/18 rounds) and a
young hill of that kind can trail on score.

| Part | |
|---|---|
| Far starts | `tr.sh_far_starts` (4) of the `tr.sh_regions` start at the best evaluated configs at least `tr.sh_far_setting` (0.2) in settings from every one of the opening's top `tr.sh_home` (20), and `tr.new_far` from the other starts |
| Indel part | per result from its class counts: the four indel terms of core with upstream's class rules (core = SNP part + indel part) |
| Indel slot | every halving keeps the `tr.sh_indel_keep` (1) regions with the best indel part among their configs with region FP at most the best config's x (1 + `tr.indel_fp_tol` 0.1) |
| Far slot | every halving keeps the `tr.sh_far_keep` (1) best far starts (best + gain) |
| Then | the reserved survivors are ordinary regions: the floor share keeps them funded, Thompson sampling decides the rest |
| Off | `tr.sh_indel_keep` 0, `tr.sh_far_starts` 0 |

## Screen ranking fixed (fleet 044200Z, 2026-10-11)

With `tr.sh_gain_weight` 1 the far starts (opened at 74.7-75.5) gained +5.9 to +7.0 in 8 minutes and took 3 of the 4
places; the best normal start (82.03 after 4 minutes, +0.08) was cut in the first halving. A low start gains points
almost automatically in its first passes, so gain since opened says nothing about the hill's top (the multi-start
climbs from far starts rose fast and then topped out at 82.4-82.6). Default now 0: rank by best; the indel slot and
the far slot still keep one indel-strong region and one far start each halving. The run was stopped at minute 12.

## Second fleet test, 044200Z (2026-10-11, 838 workers, all changes above + crashed-batch fix + legacy kings in the portfolio)

Final TRUE **82.661, rank 24 of 72** (king 82.777, gap 0.116). First v12 run: 82.715 (rank 14, 1,134 workers);
v10 with fixes 82.655; the search developer's v9.18 82.638.
- Screen (best-only ranking): kept r0 82.39, r5 82.28, r2 (indel slot) 82.17, r15 (far slot, opened at 74.4) 81.50.
- Loop: r0 reached 82.650 at minute 33 and held it; r5 82.61, r15 82.48, r2 82.35 by minute 47 (first run: live
  regions 82.46 at the same point, its 82.653 config sat in no region).
- Endgame on all four centres: +0.002 / +0.001 / +0.024 / 0; pick 82.6608 at minute 69.5. The first run's endgame
  jump (+0.179 from a second-tier region) did not recur.
- 0 configs lost to crashed batches (first run: ~3,400).
Reading: every search today on this round (40-start climbs, two legacy climbs, two v12 rounds) tops out at 82.6-82.7;
the 0.054 difference between the two v12 runs is inside the legacy project's measured run-to-run noise (0.08).
