# How the search was developed

The search proposes GATK HaplotypeCaller configs. Main submits. The evaluator runs them and scores them. This file is the history of that proposer: what each version was, what a round measured, and what the next version kept. The short score table is in [README.md](README.md). Code before v4 is in git, last together at `155b972`.

One comparison round anchors v1 through v4: **2026-10-02 09:29**, field best **89.912** of 66 submissions, 72-minute fleet runs unless noted. TRUE is the same config scored later on the released truth. One run's noise is about **0.04–0.08 TRUE**. A gap smaller than that, on one round, is not a result. A second round, **2026-10-06 01:05**, is where v4.1 was measured and v4.2 was designed. Its field best is 83.876 of 67. The two rounds are not on one scale.

The space never changed shape. About 25 published GATK parameters, three pinned (`emit_ref_confidence=NONE` because hap.py rejects GVCF output; the two adaptive-pruning floats because HaplotypeCaller never reads them), two constants, padding floored at 5. What remains is about **20 axes**. Eight are cheap inside one GATK run (genotyping ~1% extra, PairHMM ~12%). Twelve are expensive (assembly, activity, or their own run, ~55–100%). A family that differs only in a later stage shares one traversal. A result the search waits for must stay a small batch: packing thirty assembly variants onto one worker saves work and delivers the last result after ~16 run-times.

## What had to be true of every version

- Propose only. `client.pick_key` is the one way the search names a submission; main sends it.
- Rank by the validator's score. `common/scoring.py` wraps upstream. The search never recomputes the formula. From the forest onward, candidates are ranked on predicted **core** and **false positives per target** separately, then combined with the real v2 weights, because a calling threshold moves those two parts in opposite directions.
- Choose configs by score. Batching is a cost, not a reason to drop an axis.
- A failed config is an outcome (`ok=False`), never an exception.
- Judge a change on several rounds. The fake landscape and the simulator are for the controller. Replay and the fleet are the ground truth.

## 2026-09-30 — the part, then RF-UCB

`f829f47` made `search/` its own package. `EvalClient` moved out of main. A strategy is an async function of that client. The first strategy was **RF-UCB**: a random forest on the unit-cube encoding of a config, upper confidence bound, cheap genotyping variants proposed with each pick so they share a GATK run.

The same day, `ef99360` replaced the hand-written forest with scikit-learn `ExtraTrees` and pinned `emit_ref_confidence`. On the fake landscape, 200 evaluations, 6 seeds: mean best **93.2** against random **91.7**. The first real replay (round 20260928T092900Z, 20 minutes, 16 cores) evaluated **296** configs, TRUE **87.22**, rank **44/60**. Before the throughput fixes the same setup had managed 79 configs and rank 49.

`199bf59` renamed the client verbs from `submit*` to `propose*`, so "submit" only means the platform.

What this era fixed, and what it did not: the search could fill workers and attach cheap riders. It had no notion of a hill. A forest ranked absolute score per cost, which prefers a cheap genotyping tweak over an expensive axis that actually moves the calls.

## 2026-10-01 — fleets, and the pending cap

Search itself barely changed. The fleet did, and the search had to survive it. `0c866e0`: RF-UCB's pending cap had been a fixed 256, which left a 199-worker fleet at about 50 busy cores. The cap became `max(256, (family + 1) × workers)`. A 20-minute fleet replay on 20261001T044100Z (GIAB HG007, chr18, a 32-core node and a 180-core node) scored TRUE **78.47**, rank **52/68**, field best 79.77, local 78.65, pick gap 0, Spearman 1.0. The search was tracking its own answer set and still short of the field.

## 2026-10-02 to 10-04 — pattern search, then basin racing

Two strategies were tried beside the forest, then dropped as the default.

**Pattern search** (`pattern.py`, measured 2026-10-02 on k-NN landscapes built from truth-scored legacy rounds). Utility became predicted *improvement over the incumbent* per cost. Absolute score per cost had been letting any cheap variant win, so expensive axes almost never moved. Candidates included range edges (optima often sit at dangling-branch 20 or a small padding), pair moves, and a reserved share of full runs the forest could not refuse. Two parents stayed alive: highest core, and lowest false-positive rate inside a band of that core. On those stand-in landscapes, 14 workers, 45 minutes, 3 seeds: mean regret to the round's best **0.105**, against RF-UCB+TPE **0.17** with the same warm start and **0.56** without one. Most of the gain was the warm start.

**Basin racing** (`96920c6`, `basin.py`). Local score was used to decide which structural neighbourhood was worth a look, and was *not* used to pick a winner inside it. A submission changed only when a different structural signature led by **0.15**. That band came from round 20260928T092900Z: local 87.4 and 85.1 were different official bands (~88.16 vs ~86.04); inside 87.38–87.41 the official scores did not separate. Genotyping was not part of the signature. One family was one traversal; the other members only re-genotyped. Longer intensification inside one neighbourhood had not moved the official score.

Both taught the same limit. The local score is a biased stand-in for TRUE (on the later fleet rounds the bias is several tenths, sometimes a full point). Climbing it to a resolution finer than that bias does not change the rank. The 0.15 band was the first time the search refused to treat a small local gain as a new submission.

## 2026-10-04 — v1 and v2, one engine, asynchronous basins

The default became one streaming loop with no generation barrier. A round allows about 15–20 dependent steps (a GATK run is minutes under load) and as much width as the fleet has, so every step is as wide as the free workers, and an improvement is used when it arrives. The code of this engine is the `engine.py` / `basins.py` pair kept at `1602b81` behind `MINOS_ENGINE=v2`. v1 itself was overwritten the same week; the version table is the record of what it did.

**v1** (2026-10-04). Basin climbing, **24 riders per head**, first-accept. On 20261002T092900Z: **245k** configs, local 88.701, TRUE **89.645**, rank **34**. The rider count is the lesson. Genotyping riders cost ~1% and inflated the config count without adding information. 245k evaluations did not beat a later run with 103k.

**v2** (2026-10-04, still the basin engine). Best-of-sweep acceptance, valley tests, scouts adopted into basins. A sweep's centre moved to the **best** result once 80% of the sweep was back, for a gain of at least **0.005** (taken as the answer set's ranking resolution). A gain of **0.05** moved the centre at once. Other winners of the sweep were re-applied on the new centre (the legacy pipeline's accN: prefixes by gain, and pairs of the top four). A full-range probe of one axis ran at birth and when a basin stalled. Scouts were crossover and mutation of the basins' bests, ranked by a two-forest model (core and FP rate). A scout that was its own region could birth a basin. Opening was 12 minutes. Distinct regions had to be at least **4 live axes** more than 10% of range apart (`MIN_FAR`), not merely far in RMS. On the same round: **103k** configs, local 88.728, TRUE **89.736**, rank **24**. Best rank of the early basin era. Best-of-sweep beat first-accept, and stale centres stopped winning.

`498a490` (2026-10-05, landed on this engine before v3.2) added the pieces every later version still uses:

- `bam_caps`: `min_base_quality_score` above the BAM's highest base quality calls nothing. On one fleet run that was 409 zero scores.
- `AxisBook`: per-axis step, a cliff when a move crashes, and a dead axis parked after 8 identical scores, retested after 10 minutes.
- `far_axes`: the hill distance. RMS had been blurring "four knobs different" into one number.
- A designed opening wave over the dominant axes.

## 2026-10-05 — v3, v3.1, v3.2: isolate the hills

v2's scouts and basins could still wander into each other. v3 separated them.

**v3.** Isolated domains, a ParamILS-style climb, a funnel of 6 → 3 → 1 domains, shakes streamed onto spare workers. **32k** configs, local 88.661, TRUE **89.612**, rank **34**. Worse than v2. Random shakes on spare workers won nothing after the first 20 minutes.

**v3.1** (`1602b81`). The climb became iterated local search with a real mode switch. Coarse: one axis at ±1×, ±0.5×, ±2× its step, plus pairs of the strongest axes. A batch gaining under 0.05 switched the domain to **fine**: steps halved, flat axes left out, triples join the pairs, each exhausted sweep halves the scale again. Two flat batches froze the domain and its workers went elsewhere. Domains merged when they were one hill (by distance, or by a midpoint run) and the survivor inherited the loser's transplants, tested set, and finer steps. Nothing crossed between live domains. A structural move carried one genotyping sibling and one PairHMM sibling. **59k** configs, local 88.778, TRUE **89.714**, rank **25**. Fine steps paid late in the round. v2 stayed selectable with `MINOS_ENGINE=v2`.

**v3.2** (`155b972`, `engine3.py` + `domain.py`). The global search stopped ending at the opening.

- Wave 1: seeds, a portfolio, and a **Latin hypercube** over every axis (each axis stratified, pairs jointly covered), 8 genotyping riders.
- Wave 2: half the points where the model was confident, half where it was uncertain.
- A scout share for the whole round, ranked by predicted score plus uncertainty, forbidden inside a standing domain's hill.
- A freed slot (freeze, merge, or the funnel) was refilled from scout finds. Domains rolled.
- Division axes: the three live axes with the widest spread among the good configs. A config belonged to the nearest centre on those axes.

**116k** configs, local 88.771, TRUE **89.770**, rank **23**. Best TRUE and best rank on this round. The scout's finds all merged back into one hill. Read as a fact about the round, not only about the code: **one plateau, entered from many sides.**

## 2026-10-05 — v4: one controller, recombination, small batches

`0b5169e` deleted the v1–v3.2 strategies. One engine remained. The reasons were structural, from watching v3.2 and from a batching bug found while building v4.

Waited-for moves had been packed into big batches. A big batch runs on one worker, variants one after another. Thirty assembly variants saved ~45% of the work and delivered the last result after ~16 run-times, so the climb sat idle behind its own packing. v4 split `EvalJob.batch` into **small** and **big**. A result the search waits for (a climbing move, a merge child, a generation child, a midpoint) is small: one config per worker, at most 1.3 solo-runs of work, never riding in someone else's batch, sent to the fastest free node. Wave 0, the scout, and riders are big: at most 2.3 solo-runs, packed with anything that shares the traversal. Both limits are settings.

The round, as v4 first ran it:

| When | What |
|---|---|
| 0–15 min | Wave 0 (defaults, official seed, Latin hypercube, 4 threshold riders) while the answer set is built. Then up to 3 niched genetic generations. A generation launches once 60% of the previous is back, or when workers sit idle. It stops early when a generation gains under 0.05. |
| ~15 min | The hill registry founds one domain per hill, up to 6. Two configs are different hills at 4 or more live settings more than 10% of range apart. |
| 15–62 min | ParamILS in every domain. Winning settings are offered across domains every 5 minutes. The scout keeps 15% of the workers. |
| any time | Two domains on one hill merge. A merge evaluates recombination children: every combination of the settings where the centres differ, when 8 or fewer differ, otherwise 256 random combinations, plus mutated children. The best of parents and children is the new centre and inherits both tested sets. |
| 34–40 min | Scheduled merges, closest pairs first, down to 3 domains. |
| 50–62 min | The remaining domains merge into one. |
| 62–72 min | A final on that one point, then the **plateau-centre pick**: among configs within 0.05 of the best, the one whose neighbours score highest. Main submits that config, not the raw best. |

The controller is one pass over facts after every batch of results, and at least every 2 seconds. Domains do not decide their own lifecycle. They report a centre, a pool, and a progress rate (own gain per minute over the last 10 minutes). Workers are shared by that rate, with a floor so a new domain is not starved to zero.

Fleet result, same comparison round, hand-tuned during the run (`opening.max_min`, `climb.decisive`): **47k** configs, local 88.699, TRUE **89.711**, rank **25**. Inside the noise of v3.2's 89.770. Not a win. The machinery (registry, recombination, the pick, small batches) is what later versions edit.

The same day, three changes that are not a new climb but changed what the climb sees:

- `d82ab7f`. Every tunable number moved into `common/tuning.py`. `tuning.json` on main overrides it and is re-read live. A live setting applies at its next use. A round setting (`domains.max`, `riders.wave0`, the submission intervals) waits for the next round. Tests never read the file. The run manifest records the settings at the start and every change, so a hand-tuned run can be told apart afterwards. The dashboard draws `client.notes["engine"]` every 2 seconds.
- `3cc3808`, from the first v4 fleet run. A 98 GB node had been given 73 workers, ran out of memory, and looked 5× slow. Workers now keep a RAM reserve. Node speed is learned from small groups only, because big-batch cost-table error had made healthy nodes look 2–3× slow. A node slower than 3× the fleet median gets big batches only. The engine waits at most 60 seconds for nodes to join.
- `79fe84b`. Wave 0 takes up to 40% of its points from a **portfolio**: `search/portfolio/legacy.json`, 749 configs from 79 legacy rounds (2026-09-06 to 09-18), the round's read-length profile first (A = 250 bp, B = everything else), spaced at least 0.05 apart in mean normalised distance. Reviewed runs append their own near-bests to `$DATA_DIR/portfolio/ours.jsonl`, which is read first. Opening hard limit 15 minutes. The 3-to-1 merge starts at minute 50. Old scores are not reused. A new round has a new donor and region; only the configs are seeds.

## 2026-10-05 — v4.1, written from the first v4 fleet run

`e2e0df4`. The first v4 fleet run had reached perturbation never, one domain from minute 41, and 86% of the workers on a scout whose finds could no longer found a domain. v4.1 changed the climb and the find rules. The generation opening and the evidence merges stayed.

Climb:

- A fresh domain **scans every setting over its whole range** first (the legacy pipeline's injection: on its rounds, +0.2 and 7× more top-30 configs per evaluation than the opening wave). Default 9 points per numeric axis.
- A coarse sweep that comes back without a gain is followed by a **wide** all-pairs round of the strongest axes, capped (`climb.wide_max`).
- Perturbation does not stop. Shake 2, then 3, then 5, then 8 settings, one local round around the best shaken point, keep a result that is better **or equal** (ParamILS). An exact tie is not a move during a sweep: identical scores mean the setting changed nothing here.
- Parallel ParamILS **chains** on the workers the main line leaves idle, up to `climb.max_chains` (80 at the time). A chain result better than the domain by the acceptance resolution moves the centre.

Finds:

- Off-hill configs at least as good as the weakest domain, best first.
- A valley test (one run halfway to the nearest centre) may found a domain when a slot is free, or when the stage has fewer domains than its minimum (4 before minute 34, 2 before minute 50).
- With every slot full, a find 0.05 better may replace the weakest domain.
- Anything else close to the nearest domain is merged into it, all round long.
- A find is tested once. A point on a merging domain's hill is handed to the domain born from the merge.

The hill rule, once the ranker is fitted, counts only settings whose importance is at least `domains.importance_share` of the mean. v4.1 ran that share at **0.5**. The scout drops from 15% to 5% when scheduled merges begin. In the final, that share scores neighbours of the pick candidates. The schedule scales with the search's real length, node-wait included.

## 2026-10-06 01:05 — the v4.1 fleet run

Round **20261006T010500Z**. 719 workers. Schedule scale 0.99 (a 71-minute search). 287 portfolio seeds, profile B, plus a hypercube: wave 0 was 719 points. This is the run v4.2 was written from. Log: `test/runs/replay_20261006T010500Z_0452.log`.

What the controller did:

- Generation 1 launched at t+6.4 min, 654 points, while wave 0 was 579/719 back.
- Generation 2 launched at t+13.1 min, 666 points. Founding at t+14.8 min had **2 of those 666 back**. The third wave was paid for and then abandoned.
- Six domains were founded, local 82.85 to 83.16.
- Within two minutes, four of them were paused inside merges of 179 children. A merging domain clears its pool and proposes nothing. 179 is `max(children_min, 0.25 × workers)` after dedup.
- One full-range scan was logged in the whole round (domain 0, 161 configs, at t+26). No wide round, no fine sweep, no perturbation round.
- 33 domains were born. 27 merges ran. A domain was often merged again in the same controller pass that had just born it. Several of those merges resolved at 0 or 4 children back: the 179 jobs had been launched and the parent was already gone.
- At t+43 a merge fired with **"0 settings apart"**. Importance share 0.5 had dropped every axis the ranker called unimportant, so two centres became one hill.
- At t+52 the status line showed **540 scouts** in flight. The 5% quota applies only while a domain has a move to send. When every domain is merging, the fill loop's fallback gives every idle worker to the scout.

Where the local score actually moved:

| When | Best local | Source |
|---|---|---|
| t+15 min | 83.17 | Opening. |
| t+31 min | 83.73 | One axis: `min_base_quality_score` on domain 12, about **+0.41**. |
| t+40 min | 83.90 | A recombination child. The parent was paused for the merge. |
| t+43 min | 83.95 | A chain, once (`max_reads_per_alignment_start` with `min_dangling_branch_length`). |
| t+47 min | 84.10 | Another recombination child. |
| t+53 min | 84.13 | The scheduled 3-to-1 merge. One domain left. |
| t+66 min | 84.14 | A scout point during the final. |

Ten centre moves were single-axis. Three were chains. The rest of the +0.97 from the opening were merge children. The final's plateau pick selected the best config itself: neighbour mean 84.136 against best 84.138.

Review: submitted TRUE **83.554**, rank **37/67**, field best **83.876**, gap **0.323**. Pick gap on TRUE 0.001. Spearman 1.0 on the ten configs that were re-scored, because they are the same plateau. Local sat about **0.58 above TRUE**. Another ten minutes inside this plateau cannot close 0.32 on the field. Either the opening never scored the field's basin, or the answer set prefers this one. Those are different bugs. The measurement that separates them is one score of the field-best config against this round's answer set. It has not been done.

Read against the project's own noise rule: the +0.41 from `min_base_quality_score` and the +0.37 from recombination children are real on this log. The chain's +0.05 is at the noise floor. The pick did nothing.

## 2026-10-06 — v4.2, the current defaults

`04f0211`. Written from that run and from an A/B of the switches on the stand-in landscape. Every change is a setting. The v4.1 behaviour is one switch away. Tests: 194 passed, including the v4.1 path with the switches off. **No fleet round has run v4.2 yet.**

| Setting | v4.2 default | What it answers in the 01:05 run | What it risks |
|---|---|---|---|
| `opening.steady` | 1 | Generation 2 had 2/666 results at founding. Children now go to idle workers as they free up, `opening.generations × workers` in all. | The early stop watches the **global** best. A flat leader ends breeding while other islands have barely started. |
| `opening.islands` | 1 | Soft niching had been collapsing onto one hill. Children never mate across islands; a child on another island's best is redrawn; founding takes one domain per island. | Untested on a fleet. `candidates` keeps only the top 600 scores once 1,200 exist, and breeding reads that dict, so a weaker island can disappear while the leader is flat. The 01:05 log crossed that cap. |
| `domains.found_min_far` | 6 | Founded centres were close enough to be merged immediately. | |
| `domains.importance_share` | 0 | 0.5 produced "0 settings apart" and fused the round into one domain. The hill rule counts every live setting again. `tuning.json` on this machine sets the same value. | Duplicate climbs stay apart until the clock, including axes that do not change the calls. |
| `merges.on_evidence` | 0 | Evidence merges paused the domain that had just gained +0.41. Domains merge only on the schedule. An exhausted domain retires. | |
| `merges.retire` | 1 | Scheduled merges no longer run children. The stronger domain stays; the weaker's differing settings become transplant moves (each alone, plus pairs of the first four). | The +0.37 after 83.73 **came from recombination children**. Transplants cannot express a combination of five settings. This is the switch the 01:05 log does not support. |
| `merges.find_children` | 60 | A scout find is still recombined with the nearest centre, and that domain **keeps climbing**. A better child moves it. | Those children go out as big batches at priority 6, ahead of the scan. One find can still preempt the climb. There is no "one find-merge at a time" guard. The 01:05 run re-merged domains in the same tick they were born. |
| `climb.reapply_full` | 1 | A sweep's other winners are re-applied alone, then as a running sum, then as pairs of the top four. This is v3.2's accN, restored. | |
| `climb.max_chains` | 300 | Chains are meant to fill the workers the main line leaves idle. The main line is refilled first and goes out before chain moves. | The 01:05 run's only chain gain was +0.05, then the domain was merged away. v3 had already measured that random shakes on spare workers win nothing after 20 minutes. Chains accept at 0.005 and a chain recentre clears the main-line pool, including a scan. |
| `climb.eager_wide_domains` | 2 | With two domains or fewer, the first coarse sweep includes the wide all-pairs round, so one domain can fill the fleet. | On a fresh domain the first sweep is the scan, and the code then forces the wide round because `wide_done` is clear, before the coarse ±1/±0.5/±2 singles. The docstring and the code disagree about when the wide round runs. |

Before the opening, v4.2 reads the BAM and drops values of three settings that GATK cannot tell apart (`search/bamgroups.py`, the same rules as `evaluator/equivalence.py`). `min_base_quality_score` and `base_quality_score_threshold` are proposed only at the group representatives. `max_reads_per_alignment_start` stops at the most reads sharing one start; below the published minimum the axis is one value. The PairHMM threshold's list is the finest one, at the lowest legal mapping quality, and each config is then snapped to the group its own mapping quality defines. A seed that carries some other integer is snapped onto the representative, so it is the same GATK run the executor would have deduplicated.

The reduction is the product of those three integer counts, and it is not one figure for every round. Counted at mapping quality 59, where a best config on 01:05 sits:

| Round | `min_base_quality_score` | `base_quality_score_threshold` | `max_reads_per_alignment_start` | Product |
|---|---|---|---|---|
| 20261006T010500Z | 41 → 6 | 45 → 5 | 276 → 50 | 1/339 |
| 20261002T115300Z | 41 → 24 | 45 → 20 | 276 → 9 | 1/118 |
| 20261002T092900Z | 41 → 33 | 45 → 37 | 276 → 1 | 1/417 |

What the search actually proposes is a bit larger than that row, because the PairHMM list has to stay valid at every mapping quality. On 01:05 that list is 26 values, the base-quality cap removes the empty group above 30, and the product is 1/78. One five-hundredth was the rough estimate. None of these three rounds reaches it.

One cap collapses after a run, and only upward (`search/reductions.py`). The callset carries the haplotype cap the assembly actually needed, and only when no k-best search hit its cap. The floor is read from the callset of the same fidelity as the score. A later proposal above that floor is rewritten to a cap that already ran at or above the floor. Every cap at or below the floor stays, and a later report never lowers the floor. If a later report raises the floor, a stored cap that now sits below it is replaced by the cap from that report, so a higher cap is never rewritten onto a smaller assembly. `max_alternate_alleles` stays every integer from 1 to 20: the VCF lists alleles present in the genotype, so a count taken from it can miss an allele a tighter cap would drop. The calling threshold stays the full range 30–100. On the measured curves the score rises and then falls, so two equal ends can have the maximum between them.

Across the 177 stored rounds (23 Sep 2026–6 Oct 2026) the BAM removes at most one of the 20 free dimensions: `max_reads_per_alignment_start`, and only when the busiest start has fewer than 25 reads. 106 rounds go from 20 to 19. 71 stay at 20. `min_mapping_quality_score` stays 0–60 on every round; every mapping quality is present. The three-axis integer grid (the product above, at the finest PairHMM list) has median 1/430. The smallest factor in the 177 is about 1/14, the largest about 1/3900.

The seven published floats, from the GATK 4.5 source:

| Float | In the search | Resolution |
|---|---|---|
| `adaptive_pruning_initial_error_rate`, `pruning_lod_threshold` | Fixed at the GATK defaults. Adaptive pruning is off, so neither is read. | One value each. |
| `contamination_fraction_to_filter` | Still the full range 0–0.05. | Exact. GATK removes `(int)(evidenceCount * fraction)` reads. On 01:05 the region depth tops out at 224, so `[0, 0.05]` is 762 runs if every evidence count up to that depth can occur. Not installed: the evidence count is known only after assembly. |
| `standard_min_confidence_threshold_for_calling` | Full range 30–100. | The sorted phred confidences of the sites that assembly produced. On the current-best config the best score is constant from 77.988 to 78.480. |
| `active_probability_threshold` | Full range 0.001–0.05. | The Gaussian-smoothed activity at each base. Those values move with the heterozygosity prior. |
| `heterozygosity`, `indel_heterozygosity` | Full ranges. | `Math.log10` of the value is the Hardy-Weinberg prior. The callset changes when a posterior crosses the calling threshold, a genotype flips, or a base crosses the activity threshold. |

The sampler already writes floats to six decimal places. That is the sampler's rounding. GATK parses the double as written.

The README's minute-by-minute table still describes v4.1 (three generations, evidence merges, recombination on the schedule). The defaults above are what a round runs today.

## What is in the tree now

| Module | Since | Role |
|---|---|---|
| `client.py`, `base.py` | 09-30 | `EvalClient`, `Outcome`, the random baseline. |
| `space.py` | 09-30, caps 10-05 | Unit-cube encoding, pinned settings, floors, BAM caps, cost tiers. |
| `axes.py` | 10-05 | Steps, cliffs, parked axes, `far_axes`, the Latin hypercube. |
| `families.py` | 10-05 | Genotyping riders on global points. |
| `ranker.py` | 10-05 | Two ExtraTrees, core and FP rate. Orders scout candidates. Never accepts a move. |
| `portfolio.py` | 10-05 | Opening seeds. Legacy 749, then our reviewed runs. |
| `hills.py` | v4 | The one distance rule. Ownership, founding, the scout's exclusion zone, and merges all ask it. |
| `genetic.py` | v4, islands in v4.2 | Niched generations, island generations, merge children. |
| `climber.py` | v4, scan and chains in v4.1 | One domain's ParamILS. Reports facts. Does not found, merge, or retire itself. |
| `engine.py` | v4 | Opening, controller, scout, allocation, the pick. |
| `bamgroups.py` | v4.2 | The three BAM axes, reduced to the values the round can tell apart. |
| `reductions.py` | v4.2 | The haplotype cap, collapsed upward after a finished run. |

`rf_ucb.py`, `tpe.py`, `pattern.py`, `basin.py`, `basins.py`, `domain.py`, and `engine3.py` are gone from the tree. They are in git through `155b972`.

## What each measured step is worth keeping

Kept, because a round showed it:

- Best-of-sweep, not first-accept (v1 → v2).
- A gain under 0.005 is not a move. A gain of 0.05 may move the centre immediately (v2, still the defaults).
- Hill distance in axes, not RMS (v2).
- Fine steps late, and inheriting the tested set and the finer steps at a merge (v3.1).
- A Latin hypercube that actually covers every axis, and a scout that is forbidden inside a standing hill (v3.2).
- Small batches for anything the search is waiting on (v4).
- The plateau-centre pick, as a guard against answer-set blur. On 20261006T010500Z it selected the best config; the plateau was already flat.
- The full-range scan at birth (v4.1). On 01:05 it almost never ran, and the one-axis move that did run was the largest single gain of the round.
- Steady breeding into idle workers (v4.2), because a generation launched at t+13 was thrown away at t+15.
- Not pausing a domain for an evidence merge (v4.2). That pause is why the scan did not run.

Not kept, or not supported by the measurement:

- 24 riders per head (v1). Counts without information.
- Streamed random shakes as the way to use spare workers (v3). Rank went backwards.
- Treating v4's 89.711 against v3.2's 89.770 as a reason to change the climb. The gap is inside 0.08.
- Importance share 0.5 (v4.1). It erased the hill rule on 01:05.
- Replacing recombination with transplants (v4.2's `merges.retire`) **as if the 01:05 run had asked for it**. The run's second large gain was the children.

## Open, as of v4.2

1. v4.2 has not been fleet-tested. The next replay of 20261006T010500Z is the comparison: TRUE 83.55 is the number to beat, and 83.88 is the field.
2. A new domain still does not have a protected share until its first scan is back. Gain-weighted allocation and a priority-6 find-merge both cut in front of that scan. The 01:05 run finished one scan.
3. The scout fallback is uncapped when no domain has a move ready. That is the 540-scout line.
4. A chain can recentre at 0.005 and abort a scan. The decisive bar is 0.05. The scan is the higher-value work.
5. Island breeding and the top-600 candidate cap disagree. The cap was written to bound memory. It also deletes the other islands.
6. Score the field-best config of 20261006T010500Z on that round's answer set before spending another version on exploiting the local plateau.
7. `contamination_fraction_to_filter` has an exact finite resolution (762 runs on 01:05 if every evidence count up to depth 224 can occur). It is still the continuous range. The evidence count is known only after assembly, so a grid taken from raw depth is a safe superset and is not in the search.
8. The other four live floats stay continuous. Their plateaus are the data's own probabilities, and a preset step can fall between two of them.
