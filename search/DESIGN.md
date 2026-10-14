# Search engine v9: design (2026-10-08)

The search engine is now ours to build (the search developer asked for it on 2026-10-08). v9 is a new engine written
from scratch and lives in `search/`. The previous engine (v4.22) is kept in `search/v422/` as a frozen backup only:
it is not selectable, and its settings are gone from `common/tuning.py` (git commit `ad32f75` has it runnable).

Status of each part:
- **AGREED** with the user (2026-10-08).
- **PROPOSED**: my recommendation, waiting for a yes.
- **OPEN**: still to discuss.

---

## Why a new engine (the evidence it is built on)

- **The opening sets the hill; the rest only fine-tunes.** In 16 recent fleet replays, the best basin after the
  opening was already 0.20–1.21 TRUE below the field's #1. The remaining ~55 minutes added 0.0–0.6. The winning basin
  ranked 1st by opening score in only 3 of 16 runs.
- **Selection is not the problem.** The pick gap is 0–0.07, and local ranks configs like TRUE (Spearman 0.82–0.998
  on reviewed configs; 0.999 over 56k configs on R2).
- **The winners are on another hill, not a better point on ours.** On R2 (`20261006T010500Z`) the winners have
  core 0.8745 at FP 1,165. None of our 46,571 TRUE-scored configs reaches both. The legacy project concluded the same in
  September (`/minos/legacy/107_history.md`: "a different basin, not a better point on ours").
- **What paid in the current engine:**
  - recombination children between hills mid-climb: +0.55 (v4.1);
  - full-range scans: the largest contributor in recent runs;
  - single-setting steps: +0.17.

  What didn't: ParamILS chains (never beyond +0.05), the wide all-pairs round (0 moves in 4,143 configs), the scout
  (switched off).
- **Resolution matters.** Winning moves are 1–4 units. The engine's 12%-of-range steps skip them. Full integer
  resolution found +0.111 on R1 in 4 moves (`test/pairsweep`).
- **Interactions exist but are small-group.** R1 was additive near the top. On R2, a read filter
  (`min_mapping_quality_score`) broke an additive gain. Legacy found a real pair (calling threshold × mismapping rate,
  +0.036 together while each lost alone). Moves of 10+ settings rarely pay (user's experience; legacy:
  model/GA proposals moving 12–17 settings could not make the 1-setting endgame moves).
- **The local N is wrong.** The answer set fixes N at 548; true N was 270–671 in recent rounds (local − TRUE
  −1.2 … +5.1). A wrong N sets the wrong exchange rate between core and false positives. N is now the setting
  `answer_set.n`.

Sources:
- `test/fingerprint/README.md`, `test/pairsweep/README.md`, `test/additive/README.md` (our experiments);
- `/minos/legacy` (read only);
- `/minos/minos-miner/test/` and `search/DEVELOPMENT.md` (the search developer's workspace, read only).

---

## Step A: candidates

**A1. BAM analysis. AGREED.**
- Use only **exact** reductions from the BAM (`search/bamgroups.py`): values proven to give the same GATK run.
- Never cap a range from belief. The legacy project capped `min_base_quality_score` at 30 from one 100 bp round and
  then never found the mbq-36 family on 148 bp rounds.
- No crash model. About 10% of opening configs crash, and the rate drops sharply once the search works near good
  configs. Failed configs are dropped and never retried (GATK is deterministic).

**A2. Portfolio. AGREED.**
- Starting configs come from the portfolio: legacy 749 plus ours, the round's read-length profile first.
  Measured value: a running start, +0.25 at the start, about +0.16 by the end. Configs transplanted from other rounds
  are not answers: they score 1.2–2.5 below finals.
- **Dedupe against the whole portfolio,** not per round. Legacy had 22–37% near-duplicates.
- **A replay never seeds configs from its own round.** v4.2's "83.83 at minute 5" on R2 was most likely that
  leakage.
- **Three groups:**
  - **Both:** the Pareto front on (core, FP), i.e. configs no other config beats on both.
  - **Core:** highest core within 0.01 of the best core.
  - **FP: fewest false positives** among configs whose core is within 0.01 of the best (AGREED 2026-10-08).

  Note: the current code's "fp" group keeps the **most** false positives (`engine.remember("fp", fp)` sorted
  highest-first; portfolio `FP_BAND` 0.25 below the maximum).

  Why fewest:
  1. The winners sit at the low-FP end of the frontier (R2: same core, about 50 fewer FPs).
  2. Choosing parents by core and FP counts does not depend on N, which is off by up to ±50%.
  3. Score-only selection collapses the population (legacy: 90% of members at `min_dangling_branch_length`
     19–20).

  Limit: crossing a core parent with an FP parent often lands in between ("a core-rich block always brings its FP",
  legacy). The groups mainly keep both ends of the frontier alive. Unmeasured, so it is a setting.

**A3. GA, one generation (2026-10-08).** Wave 0, then generation 1 breeds **within** each group.
- Generations 2 and 3 are not run. On `20261008T152900Z` generation 1 raised the best by 0.155 and generation 2
  by nothing; generation 3 added 0.037, inside one-pair noise. On `20261008T175300Z` generation 1 added 0.184 and
  generation 2 had not moved the best when the clock cut it.
- **Diverse parents.** Legacy crossover collapsed from 792 new children to 123 when later generations bred only
  from the top.
- Mutation: **30% of children get 2 settings redrawn over the full range** (AGREED 2026-10-08).
- Step C still breeds the cross (core × fp) after the climb.

**A3b. Generation 1 is waited out (2026-10-08).** Wave 0 is bounded by `candidates.max_min` (20 minutes, scaled):
on a heavy region it took about 10 minutes to come back. Generation 1 then gets `candidates.generation_min`
(12 minutes, scaled) of its own and is waited until 80% of it is back. The wave-0 clock does not abandon it.
`climb.end_min` is 55, so Step C has about 15 minutes before the pick. Ten minutes is one fleet wave with nothing
left to climb.

**A3c. Parents are the top 20 (v9.3).** Generation 1 builds the both / core / fp groups from the 20 highest scores. A config 0.10 or more below the best does not get a seat, so the core extreme that sat at 83.93 on `20261009T044200Z` is not one of the three climbers. Fewer than three seats run when fewer than three configs are that close and far apart.

**A3d. Sixty parents, kept far apart (v9.7).** The parent list holds 60 scores. A config is added only when it is at least `candidates.min_distance` from every parent already on the list (callset distance, or the mean setting difference when a sketch is missing). A higher score on a hill already in the list is skipped. The both / core / fp groups are built inside those 60, and Step B breeds from them. Twenty was one hill: on `20261009T093000Z` and `20261009T104200Z` the children came back too close to the climbers.

**A4. Pick candidates far apart. AGREED.**
- **Count:** 1 candidate per 500 workers for Step B.
- **Distance:** callset distance (AGREED 2026-10-08), i.e. how many variants one config calls that the other doesn't, weighted
  by the score's class weights. Identical callsets are distance 0 whatever their settings.
  - Every settings-based distance failed (unused settings and cliffs distort it). Recent runs founded 12–24 centres
    "≥ 6 settings apart", all within 0.5 local on one plateau.
  - Two plateau configs disagreed on about 200 called variants.
  - Cost: each node attaches a ~1 KB MinHash sketch of the callset (per class) to every result. This needs an
    evaluator change.

---

## Step B: climb each candidate to its peak. AGREED

The most important step. Find each candidate's local peak as fast and correctly as possible. Local only: Step C
handles far points.

**B1. The grid.**
- Settings move by **δ = their resolution**, independent of range:
  - integers: 1;
  - `min_base_quality_score` / `base_quality_score_threshold`: the next BAM group;
  - enums and bools: the neighbouring value;
  - floats: the agreed grids (`test/pairsweep/sweep.py: values`).
- Each gridded setting takes **−δ, 0, +δ**. The grid is **every combination** over a group of 5 settings: 3⁵ = 243
  points, fewer when the BAM collapses a setting (R1: `max_reads_per_alignment_start` has one value, so 3⁴ = 81).
  A true global grid over 20 settings (3²⁰ ≈ 3.5 billion) cannot run, and moves of 10+ settings rarely pay or are
  sums of small-group moves.
- **Expensive groups** (settings that cannot share a GATK run, so each point is a full run):
  - odd iterations: the 5 read-pass settings (`min_base_quality_score`, `min_mapping_quality_score`,
    `assembly_region_padding`, `heterozygosity`, `max_reads_per_alignment_start`);
  - even iterations: 5 of the 7 assembly/activity settings (`min_pruning`, `min_dangling_branch_length`,
    `recover_all_dangling_branches`, `max_num_haplotypes_in_population`, `active_probability_threshold`,
    `min_assembly_region_size`, `max_assembly_region_size`), the 5 rotating.
  - Why the read-pass group: these are explored least today (every value is a full run), and 3 of the 4 settings
    that separated the legacy winner from our configs were in it (mapping quality, max reads, mbq).
- **v9.2 riders.** Every grid point carries the 4 genotyping settings at ±δ inside its own GATK run (~1% each,
  so they fit and the run stays one run). PairHMM (±δ), every pair of the cheap settings, and the rest of each
  cheap ladder stay on the **centre only**, and are sent only after every main point of the grid is already
  running, packed into as few runs as one big batch holds. On `20261009T033000Z` PairHMM on every point missed
  the 1.3-run cap and each leftover started the heavy work again: one read grid took 12.6 minutes and three
  centres were run 1,197, 835 and 456 times. The grid itself is about one GATK run of wall time (~3 minutes
  on this fleet) plus the slow tail.
- **Cost per iteration:** one run per main point (about 160 on this BAM, 243 when every axis can step both ways),
  plus the centre's packed cheap batches on workers those points leave free.

**B1b. Lines (added 2026-10-08 after the first fleet test).** Each grid wave also sends every grid setting alone over
the rest of its ladder (`search/lines.py`): the gaps step growth jumped over and rungs beyond the nearest tried value,
as many as the climber's share of the fleet leaves room for; every rung of the cheap settings rides on the centre's
run. Why: on `20261008T175300Z` (true answer set, rank 40 of 72) `min_mapping_quality_score` zig-zagged. At climber 1's
peak (38) both neighbours lost (37 −0.063, 39 −0.036) while 42–48 won (+0.06 to +0.14); at climber 0's peak (46),
47 lost, 48 won (+0.028), 49 lost. A ±1 grid takes any zig-zag for a peak; only Step C found these, 20–30 minutes
later. A line winner moves the centre like any point but starts no step growth.

**B2. The move.**
- The **best grid point** (a move of up to 5 settings together, with its riders) becomes the new centre.
- **No combining** of winners from different points or groups: they were never evaluated together, and combining
  messes up the steps (R2: two additive winners plus mapping quality gave +0.004 instead of +0.042).
- **No separate x2/x3 wave** (corrected 2026-10-08): x2/x3 is the step growth of B3, applied in the next grid.

**B3. Step growth (across iterations).** v9.1, 2026-10-08.
- A setting whose move just won uses a larger step in the same direction in the next grid. The opposite direction
  stays at δ.
- A long ladder (30 rungs or more: mapping quality, reads, assembly sizes, haplotypes, padding, and the other wide
  axes) grows **1δ, 2δ, 3δ, 5δ**, and stays at 5δ while it keeps winning. A short ladder still grows 1δ, 2δ, 3δ, 4δ…
- When a grown step stops winning, that setting returns to δ.
- Example: A+1 wins the grid over (A, B, C); the next grid is around (A+1, B, C) with A at −1, 0, +2 and B, C at
  −1, 0, +1.

**B4. Randomness. AGREED: 20 random jumps.** When the last read-pass iteration and the last assembly iteration have
both gained less than 0.005 (a stall), 20 random local jumps from the centre: 2–4 settings, each 2–5 rungs (a bool
flips). A jump gaining at least 0.005 becomes the centre and the grid resumes; otherwise the centre is the peak.

**B5. Stop.** When an iteration's best gain is **< 0.005 local** (the resolution at which our answer set still
orders configs) for both grid kinds in a row and the jumps find nothing. A smaller gain still moves the centre.
v9.1: that peak does not start Step C. The seat climbs the next Step A result (the same distance rule as the
pick) until `climb.end_min`. Step C starts only after every seat has left Step B.

**B7. Waiting (built 2026-10-08, `climb.wait_riders`, default 0).** A wave decides once its grid points are back;
riders still running count at the climber's next decision. Waiting for every rider made each grid wave 5.8 minutes
instead of 3.0 in the simulator (500 workers): the riders that do not fit a small batch go to big batches of up to
2.3 runs. Same climb, 9 iterations instead of 6 in Step B.

**B6. Time.** One iteration = one grid wave, about one run time plus stragglers.
A candidate still climbing at `climb.end_min` (minute 55, scaled) stops. Its seat does not start another hill.
Step C starts then.

**B8. Free workers run the GA (2026-10-08).** A value is one GATK run, one worker. Whenever a worker is free and a
climbing candidate still has an unsent grid or line value, that value is sent. The line is the gaps and one step
past the nearest tried rung; it does not thicken until the fleet is full. A hill that finishes hands its workers
to the next Step A result. Every worker still free runs the GA, through the whole climb and through the final
fine-tune. During Step B that is generation 1's within-group breed. During the fine-tune it is the core × fp
cross. The GA is the lowest priority: main points first, then the centre's cheap changes, then the GA.

**B9. The GA's best goes to the lowest climber (v9.3).** When a GA child is the best score, and no climber is already on that config, the climber with the lowest score drops its unfinished grid and climbs the GA config. Runs already out still score. That climber's score becomes the GA score, so the next GA best is handed to whichever climber is then lowest. On the finished 03:30 round the GA led the grids by 0.03 and once jumped 0.07, then the fine-tune passed it; the grids never climbed the GA's config.

**B10. A flat seat gives its workers away (v9.3).** A gain clears the flat mark on the settings in that grid only. The last two grids under `climb.accept_min`, on a seat 0.10 or more behind the best seat, stop that climber. The next Step A candidate takes the seat. On `20261009T044200Z` the third seat opened 0.52 below the best, one later gain of 0.025 forgot its flat read grids, and it held a third of the workers until minute 54. No backup candidate climbed.

**B11. A climber's config is not climbed twice (v9.4).** A config that matches a climber's centre, differs from it in at most two settings, or has the same callset (callset distance under `candidates.min_distance`) is not sent and is not handed to another climber. The climber's own grid still runs. On `20261009T055400Z` the three centres differed only in max alternate alleles and haplotypes, scored the same callset, and each kept a full GATK run.

**B16. Four settings apart (v9.8).** The setting-count rule is four. A config that differs from a seated centre in fewer than four searched settings is the same hill: it is not sent, not handed to another climber, and not used to fill an empty seat. The climber's own grid still runs. On `20261009T130600Z` climber 1 and climber 2 differed in three settings and called the same variants.

**B17. Step B mutation rises from 0.3 to 0.7 (v9.8).** Step A and Step C still mutate `candidates.mutate_frac` of their children (0.3). During Step B that share starts at 0.3 and rises in a straight line to 0.7 at `climb.end_min`. A mutated child still moves `candidates.mutate_axes` settings by one rung.

**B18. Lowest best, and five climbers (v9.9).** A GA child that is the best score is climbed by the climber whose own runs have the lowest best. A run that is already back counts even while that climber's centre is still the older config. The handoff waits until the GA score beats every climber's best. `candidates.workers_per_candidate` is 194, so a 972-worker fleet seats five climbers at once. On the chr22 climb the three climbers together held about 160 cores and the GA held about 540.

**B21. Four climbers (v9.11).** `candidates.workers_per_candidate` is 243, so a 972-worker fleet seats four climbers. On chr15 the fifth seat started 0.20 below the Step A best and its first grid took about 7 minutes.

**B22. A higher GA score moves the lowest climber immediately (v9.12).** Step B checks on every pass. When the highest GA score is above the lowest climber's own best, that climber drops its unfinished grid and climbs the GA config. It does not wait for a read grid or an assembly grid to end. The same-hill rule still applies, and a climber already on that config stays there.

**B23. The next free GA still moves the lowest climber (v9.13).** The handoff walks GA scores from the top. A config already a climber's centre is skipped. The lowest climber takes the first remaining score that is higher than its own best, unless that config is the same hill as another climber. On chr13 v9.12 the overall GA best was already seated, so a free GA 0.12 above the lowest climber was left unclimbed.

**B24. Step B opens with GA only (v9.14).** The first `climb.ga_open_min` minutes (8, scaled) send the whole fleet to the GA. No opening grid. The window ends sooner when a free GA score beats the lowest climber's best. Then each seat starts its read grid from that centre. On the finished rounds a G score beat Step A in 1–3 minutes (6–11 on a long region) while the first grid was still occupying the fleet.

**B25. Two GA children cannot share a hill (v9.15).** A GA child is not sent when it differs in fewer than four searched settings from another GA child already sent, or when it is the same callset. The dim rule already applied to climber centres; it did not compare GA children to each other, so one hill filled the dashboard's top three.

**B26. Dim-aware elitist leftover GA. Built, v9.17.** Step B leftover (including the GA-open window) still breeds on every free worker. Half of those children (`candidates.elite_frac` = 0.5) are **elite-anchored**; the other half stay the current within-group pairs.

A new best is already a parent on the next leftover batch (the 60 far parents rebuild when the best score changes). It is only one of 60, so it is a parent about 1/60 of the time. Elite-anchored makes it one parent of every other child, which is keep-best / elitist recombination on an expensive evaluator: the elite is never re-queued (it is already scored); it only donates genes.

How an elite-anchored child is built:

1. Parent A is the current best scored config. Parent B is drawn from the far parent groups, not A, the same both / core / fp plan as today.
2. Uniform crossover as today.
3. The child must differ in **at least four searched settings** from A. If crossover left it on A's hill, copy further settings from B until it does, then one-rung mutate one extra movable axis. That is the dim rule written into the operator: a clone of the elite is refused by v9.15, so the operator must leave the hill or the fleet fills with near-duplicates.
4. `accept` is unchanged: same hill as a seated climber, closer than `candidates.min_distance`, a tied callset, or the same hill as another GA child already sent, is not sent.

Not used: a second Step A generation; copying elites into the queue; replace-worst of a small population; NSGA-II. Step A generation 1 and Step C's two generations stay as they are. A plateau still sends one-rung probes first; leftover workers after that use this mix.

The risk is the same as replace-worst takeover in the literature (about 5–6 population-sizes of evaluations): too many children of one elite. The dim rule and `elite_frac` 0.5 are the cap. If elite-anchored children are almost all refused, the far stream still fills the fleet. Measure as always: TRUE, rank, gap; noise on one pair is 0.04–0.08. Oracle tests run v9.16 first (B27 only), then v9.17 (B28 compass-then-cube, B27 far seat, B26 elite leftover).

**B28. Compass then a cube on the live axes. Built, v9.17, ships in v9.18 with B29.** Each Step B iteration used to send the 5-setting cube (up to 243 full GATK runs) and wait for every head (`climb.quorum` 1.0). That was 5–11 minutes per iteration on the finished climbs. v9.17 sends the **compass** first: each of the five grid settings at −δ and +δ, the others held (~10 runs). An axis is **live** when a face scores strictly above the centre. Zero or one live axis ends the iteration: the best point of the compass, the lines, and the riders becomes the centre, as before. Two or more live axes send the **3^m cube on those axes only**, skipping the centre and the faces already scored (2 axes → 8 new points). The move is the best of compass + cube + lines. Winners from different faces are not added without that joint eval (B2 / R2). A smaller gain than `climb.accept_min` still moves the centre (B5); `accept_min` still marks a stall. Lines and the centre's cheap packs stay as they are. The full 243-point cube is not sent.

**B29. Leftover cores keep running GA. Built, v9.18.** Compass leaves most of the fleet free. Those cores are supposed to run leftover GA. Three things left them idle on v9.17:

1. B25 dimmed a new GA child against **every GA child already sent this round**. After ~10k children almost every new kid was “the same hill as another GA child.” Exact keys stay in `fresh`; tied values stay in the tie table. v9.18 dims only against **climber centres and GA children still in flight**. A finished GATK run is not a seat.
2. The 60 far parents rebuilt whenever 50 more configs had scored, even when the best had not moved. Each rebuild was about a minute of callset distance with the loop stuck, so cores drained. v9.18 rebuilds that list **only when the best score changes**.
3. Leftover GA ran only on a tick where no climber closed a wave. v9.18 sends leftover GA on the **same pass** as the compass, after every climber head and centre pack is queued, for every remaining `idle_slots`.

Heads still go first. Then enough GA children to take every free core.

**B53. The compass stays on the highest score. Built, v9.30.** Every earlier version paid for one of these and then gave it up. The resolution grid (B1) and the compass-then-cube (B28) find a one-setting move and then check that two winning faces work together. The joint cube has to be scored; announcing it and closing the wave drops the interaction (B2). The centre stays through that open grid and moves to the best point the grid scored (B51). The next compass starts on the highest score (B37, B45), and a new fingerprint gets the one-neighbour valley first (B52). One climber (the fifth seat and the four-seat valley both spent the round on centres the fine-tune then ignored). Free workers breed (B8, B29, B49). The line does not sit in the gate (B50 and the five-rung widen in B41 filled the fleet and the compass never closed). A quiet window does not silence the compass (B24, B46). Step B repeats that compass until the fine-tune, instead of one view and a crash. The fine-tune is the same compass on the best for the last eight minutes. Quadratic dummies, callset piles, exact BAM groups, and the portfolio stay.

**B52. Separate mountains. Built, v9.29.** Step A fits a quadratic in the two strongest directions of the nearest samples around each apparent peak. The crest of a downward bowl, snapped to legal rungs, is that peak's dummy; an upward bowl sends nothing. One climber sits on the highest fingerprint. The other piles wait. A fingerprint the GA finds above that climber is checked only against that seat. A scored config lies on that callset path when both legs are shorter than the direct distance and the legs add up to at most that distance plus 0.02. A score on the path at or below the lower centre minus 0.10 is the valley, and the climber moves onto the higher fingerprint. The other pile waits. When the middle of the path (0.25 to 0.75) has no scored config, one midpoint is sent. A midpoint inside the band, below the lower centre and above lower minus 0.10, gets one quarter from the higher centre, and that is the last run. A quarter still inside the band is no valley. No valley, and the callsets are farther than 0.10, sends the lower fingerprint away before the ascent; its parents are not absorbed. Callsets within 0.10 meet: the higher pile's radius grows to cover the lower parents, the parent cap rises from 32 to 64, 96, and 128, and the lower climber leaves. A midpoint whose callset is not between the two ends leaves both in place. A fingerprint the GA finds above the lowest climber uses this same check. A dip moves that lowest climber onto the new fingerprint. A pair already decided is not reopened. Step B breeds children that copy the centre and change one to three settings taken from the parents. The grid is a compass on those settings: one rung either side, every other setting held. A face that does not beat the centre is not walked further. Two or more live faces send their joint configs. The centre stays until that grid is finished, then moves to the best point the grid scored. The fine-tune's gate is the compass alone. A child that only opened the view does not become the centre. A finished grid with nothing higher stops that climber. The fine-tune does the same: a higher score waits until the open grid finishes, then the next grid starts from the round's best. The midpoint, the view, and the fine-tune are small batches. While they are out, every free worker breeds inside the piles, and those children stay in flight when the batch returns.

**B51. The grid writes the setting. Built, v9.28.** A new high inside a pile was replacing that pile's fingerprint and the climber's centre, and the grid restarted from that config. The settings the grid had been walking were dropped. While a grid is in flight, a higher score joins the pile and the fingerprint stays. The centre moves when the wave finishes, to the setting the grid found. While the GA still owns the fleet, the highest score within 0.10 still becomes the centre, so the grid starts from it.

**B50. A free worker walks the line. Built, v9.27.** On `20261010T235400Z` the fleet sat near 400 of 1134 workers. The planned grid is one rung, and the workers it does not fill were left idle. Those workers walk further up the same line. The wave still closes on the planned heads.

**B49. Elite probes are not the same hill. Built, v9.27.** A probe copies the best config and changes one to three settings. Two such probes differ from each other in fewer than four settings, and a child within four settings of one already in flight was refused. On `20261010T235400Z` each pass sent a handful and left out about 20,000, so the fleet fell to a few dozen busy workers and the same pass ran again. During the elite window a probe of the best is sent. An exact duplicate is still refused.

**B48. Elite probes, then the grid. Built, v9.26.** Elite reproduction is there to find which settings of the best config to change, and then stop. Two things kept it going until the fine-tune. A child within four settings of a climber was refused, so the children that actually probe the best config never left. And every candidate that took a climber restarted the four quiet minutes, so the fleet never stayed on the grid. While the GA owns the fleet, every child copies the best config and takes one to three settings from a mate. Those children are sent. The quiet window starts once, when the opening hands the fleet over, and a later candidate does not restart it. After it, Step B does not breed.

**B47. A peaked centre is not free. Built, v9.26.** The occupied-centre guard only counted a climber who was still climbing. On `20261010T224200Z` climber 3 peaked on `5b94dd3c`. That config was still the highest member of a candidate, so every new climber was moved onto it, joined climber 3, and stopped. The seat was filled and the swap ran again. A centre held by a peaked climber is now skipped the same way as one a climber is still climbing.

**B46. Early Step B gives the climbers' workers to the GA. Built, v9.25.** While a GA child or a candidate is still moving a climber, a grid is a blind search: the next update throws it away. For four scaled minutes after the last such move the climbers send nothing, compass included, and the fleet breeds. When that goes quiet, the grid is one rung past the compass and does not widen. The fine-tune is never held back.

**B45. The fine-tune restarts from a higher score. Built, v9.24.** Step C reads the best config once and climbs it. On `20261010T213000Z` the fine-tune opened on 85.6164. Climber 204's last grid was still in flight, and 85.6251 came back three minutes later. The fine-tune stayed on 85.6164 and peaked. The pick submitted the higher config with about 40 seconds left, and that config was never climbed. The fine-tune now leaves its grid and starts again from the round's best whenever a higher score arrives.

**B44. A far config takes the lowest climber without waiting for a peak. Built, v9.24.** A config at least 0.10 from every fingerprint, above the lowest climbing pile, and above the lowest candidate, replaces that candidate. The swap then moves the climber. The pile it leaves waits as a candidate. On the finished pile rounds the climbers were still marked stopped when the clock ended, so a new fingerprint that waited for a peak never got a seat.

**B43. A candidate does not take an occupied centre. Built, v9.24.** The swap moves the lowest climber onto a candidate's highest config that no climber already holds. On `20261010T213000Z` pile 6's highest config was climber 0's centre. Each swap landed on that centre, the new climber joined and stopped, and the empty seat was filled and swapped again. Climbers 1 through 145 were that loop.

**B42. Step C does not breed. Built, v9.23.** The two generations before the fine-tune bred only the first pile, and on the finished pile rounds they did not move the best. That time stays with the climbers. A peaked seat still goes to the best candidate. The last eight scaled minutes are one climb of the best config. Idle cores in that climb do not breed, so a GA child cannot pull the fine-tune off the winner. Leftover GA during Step B is unchanged.

**B41. The second half of Step B widens the grid. Built, v9.23.** The first half keeps the short line: one rung past the compass, so most cores stay on the GA. From the midpoint of Step B until the fine-tune, each side of a grid setting takes five rungs. The climber does more of the ladder, and the GA takes the cores that leaves.

**B40. A higher candidate takes the lowest climber. Built, v9.23.** When a candidate's highest score is above the lowest climbing pile, that climber moves to the candidate and develops from its highest config. The pile it left waits as a candidate. A candidate that does not beat the lowest pile waits until a climber reaches the top.

**B39. Six piles, two of them candidates. Built, v9.22.** The search keeps six fingerprints. Four have a climber. The other two wait with no climber. When a climber reaches the top, the higher candidate takes that seat. The finished pile stays in the history, and its slot takes the next fingerprint that is at least 0.10 from every pile.

**B38. A dead read grid is not repeated. Built, v9.22.** The read grid is the same five settings every odd turn. Once each of them has failed since its own last gain, the next grid is the next assembly group. The checklist still finishes only when every read setting and every assembly setting is dead, and the jump round is unchanged.

**B37. The climber's centre is the highest score inside the pile. Built, v9.22.** A config within 0.10 of a pile becomes that climber's centre when its score is higher, even when the config also sits within 0.10 of another pile. The fingerprint moves only when the new centre would stay at least 0.10 from every other pile. A config farther than 0.10 from the fingerprint does not become the centre.

**B36. The highest score within 0.10 is the reference. Built, v9.21.** A config within 0.10 of a pile replaces that pile's reference when its score is higher, and that pile's climber develops from it. A config farther than 0.10 can open a pile. A move that would sit inside another pile does not replace the reference.

**B35. Mating pool of 32. Built, v9.20.** A pile breeds from up to 32 parents that differ in four or more settings. On the measured round the flat pile had 33 eligible parents inside the 0.05 ball, and the 32nd was 0.049 behind the centre. Eight kept only the first 0.015 of that plateau.

**B34. Fingerprint history. Built, v9.20.** Every fingerprint a pile stands on is logged, and it stays logged after the pile moves on. A peaked fingerprint is marked closed. A config within 0.10 of any remembered fingerprint, including a closed one, does not open a pile and is not seated. The pile's own live trail is the climb still in progress. The same pile, and a similar one, is not explored again.

**B33. A peaked fingerprint is not climbed again. Built, v9.20.** Once a climber reaches the peak of a pile's fingerprint, that fingerprint is not seated again. A different config inside the same ball is seated only when it beats that peak by at least 0.005. The workers that pile would have used stay on the piles that are still climbing, and on the search for a fingerprint at least 0.10 from every pile.

**B32. A higher returned run takes the centre. Built, v9.20.** A config this climber sent that comes back above the centre becomes the centre, including one that lands after a GA swap or after the grid has already closed. The pile's fingerprint moves with it when the 0.10 gap to every other pile still holds. A return that would sit inside another pile does not move the centre.

**B31. Four callset piles. Built, v9.19.** After wave 0, scored configs form up to four piles. A pile's fingerprint is its best config's callset sketch. Distance 0.05 or less joins the nearest pile, and a higher score replaces that fingerprint when the new centre stays at least 0.10 from every other pile. Distance from 0.05 to 0.10 opens nothing. Distance 0.10 or more, and within 1.0 of the best, opens a pile. One climber per pile. Inside a pile the mating pool is at most eight configs that differ in four or more settings; fewer than two parents mutates the pile's best. Half of Step B's leftover GA breeds inside the piles. The other half mutates six settings, looking for a fingerprint that is at least 0.10 from every pile. A climber does not take a step that leaves its ball or closes the 0.10 gap. The pick is the highest local score.

**B30. A far seat stays a far hill. Built, v9.18.** On v9.17 the far seat opened 0.20 behind, then `adopt_ga_best` gave it every free GA because it was the lowest climber. Those children were the best hill. Both seats ended on the same score; the far basin was gone. v9.18 marks that seat `origin=far` from the pick (and from a later far replacement). A GA child on a near climber's hill (four settings, or callset / setting distance) is taken only by a **near** seat. A far seat takes a GA child only when that child is far from every near centre. Near leftover GA still elite-anchors (B26). The far seat keeps climbing its own basin.

**B27. One seat is the best far hill. Built, v9.16, ships in v9.17 with B26.** After GA-open, Step B still seats `workers // candidates.workers_per_candidate` climbers. All but one are the current pick: best scores, pairwise far, within `SEAT_GAP` 0.10 of the best. The remaining seat is the best scored config that is a different hill from every one of those (four searched settings and callset / setting distance), even when it is more than 0.10 behind. It still has to be within 1.0 of the best, the same cap as an empty-seat fill (B15). If every opening climber is still the same hill after GA-open, the lowest takes that far config.

If nothing scored qualifies, scouts for that seat start at the far ends of the 60 parents (the best core parent and the fewest-FP parent), not at the global best centre. Each scout still moves six settings by about a quarter of each ladder. B15's empty-seat scouts from the best centre stay as they are for a later empty seat.

Leftover GA pairing stays as in v9.15 (no elite-anchored children). Success is a climber whose callset is far from the Step A best and a TRUE gap that moves by more than 0.1 across three rounds. If that seat stays ~0.3 behind and never catches up, the far basin was junk and the seat rule drops.

**B20. A grid end hands the lowest climber to the GA (v9.10).** When a read grid or an assembly grid ends in Step B, the climber whose own runs have the lowest best takes the highest GA score if that score is higher. The GA child does not have to be the best score in the round, and it does not have to beat the other climbers. A jump wave does not hand off. One climber already on that config stays there.

**B19. Step C breeds twice (v9.9).** After the hills finish, Step C runs two GA generations, then the fine-tune. Both generations cross the core parents with the false-positive parents. Generation 1 mutates 0.3 of its children. Generation 2 mutates 0.6. The first generation is given half of the time before the fine-tune's 8 minutes, and the second breeds from the scores that came back, using whatever time is left. The fine-tune is still the climb of the single best config.

**B12. A plateau is probed one rung at a time (v9.4).** Once two configs share the best score, Step B's GA stops breeding. A setting that already varies inside that score is the same callset and is not sent. Each setting that still has one value is probed by the climber that holds the score: one rung, every other setting held. The same score makes that rung inert. A loss of 0.10 closes that direction. A smaller loss closes that rung and the next probe steps one further. A mutation, before a plateau and inside the fine-tune cross, moves one rung rather than a random rung of the whole ladder. On `20261009T055400Z` the idle GA scored about 14,000 runs after minute 20 and the best stayed 87.2721; a random calling threshold and a random padding came back 0.011 and 0.078 short.

**B15. An empty seat takes the best far config (v9.7).** When a climber leaves and the Step A backups are all close to a climber that is still seated, the seat goes to the best scored config that is far from every centre and within 1.0 of the best. Callset distance decides when both configs have a sketch. A config within two settings of a centre is left out. When nothing scored qualifies, up to two waves of scouts are sent. A scout starts at the best centre and moves six settings by about a quarter of each ladder, so it clears the seat distance, and the rest of that config stays. On `20261009T093000Z` the second peak found no far Step A backup, and the GA children being bred were the neighbourhood of the same hill.

**B14. A plateau does not idle the fleet (v9.6).** The one-rung probes of a shared best score still go out once, on the lowest climber that holds it. The workers those probes do not use keep breeding, including after the probes are finished. A child is sent only when it is at least `candidates.min_distance` from every seated climber (callset distance, or the mean setting difference when the child has no sketch yet). A child that repeats a tied value, or that sits within two settings of a seated climber, is still not sent. On `20261009T093000Z` the shared score 86.6074 stopped breeding at minute 35 and about 900 workers sat idle until the climbers on that score peaked.

**B13. A tied value is not sent again (v9.5).** A value becomes tied when two scored configs share an exact score and that value is the one setting still unexplained between them. A later config is not sent when every value it changes has already tied against that scored config. An untried value of the same setting is still sent. On `20261009T070600Z` half the scores repeated a callset. Refusing every other value of a tied setting also refused climbs of 0.05 to 0.11, so only the value that tied is refused.

---

## Step C: GA, then one fine-tune. v9.1, 2026-10-08

Starts only after Step B has left every seat. It does not climb a new hill while Step B is still running.

**C0.** One GA generation, bred the same way as a cross (core × fp, and the Pareto group × either), over every
config Step A and Step B scored. The fine-tune keeps 8 minutes; the generation uses what is left before that.
The best scored config is then fine-tuned with Step B's climb until the pick. A second generation does not fit
in the 15 minutes.

The single-setting pass below is not scheduled. Step B's lines use the same rung plan. A grid only tries a setting
at −δ and +δ (or the grown step) from the centre, so along each single setting it leaves two kinds of untried values:

- **Case 1, gaps:** values step growth jumped over. The climb went 3 → 4 → 6 → 9, 9 lost, so 6 is the peak; 7 and 8
  were never tried.
- **Case 2, beyond:** values past the nearest tried value. The climb stopped at 3 because 4 was worse, but 8 may be
  better than 3.

**C1. What has been tried** comes from the results, not from bookkeeping: for each setting, the configs that differ
from the peak in that setting only.

**C2. One wave per peak:** every untried rung between the peak and the nearest tried value on each side (case 1), and
`far.coarse` = 25 rungs spread over the rest of each side (case 2), made denser while the wave still fits the fleet's
share for that peak. The cheap settings take every rung (they share the peak's GATK run).

**C3. Then:**
- a point that beats the peak by 0.005 becomes the centre, Step B climbs from there, and Step C runs again;
- otherwise, a second high point on a setting's profile (higher than its neighbours, a lower point between it and the
  peak, within `far.second_margin` = 0.1 of the peak) gets every rung around it tried; the best 3 per peak
  (`far.landings`) get a short climb of 2 grids (`far.landing_iterations`), and one that ends above the peak becomes
  the centre, Step B climbs on, Step C runs again;
- otherwise the peak is final.

Time left after that goes to climbing from the best config off every hill climbed so far, until the pick.

Not built: the same over pairs of settings (only if single settings turn out to find things). Dropped: recombining
peaks with each other (legacy closed it: 0 of 114 children beat their parents).

Simulator, 500 workers: the first pass at minute 50 sent 663 single-setting points and found
`max_assembly_region_size` → 266 (+0.052) that Step B's grid had not reached; Step B then climbed on to 95.47.

---

## Structure (built)

- `search/engine.py`: A -> B -> C, the schedule (written for 72 minutes, scaled to the round), the pick, and
  `client.notes["engine"]` (version, phase, basins: one row per climber with `first_key`/`best_key` for the review).
- `search/candidates.py` (Step A), `search/climb.py` (Step B), `search/far.py` (Step C).
- `search/steps.py`: each setting's ladder at its resolution (δ = one rung).
- `search/distance.py`: callset distance from the sketches the evaluator attaches to every CallSet
  (`evaluator/sketch.py`, `CallSet.sketch`); mean normalised setting distance when a sketch is missing.
- Shared, kept: `client.py`, `space.py`, `bamgroups.py`, `portfolio.py`, `axes.py` (`latin_hypercube`), `base.py`.
- `search/v422/`: the previous engine, frozen backup, not selectable.
- Interfaces kept: main calls `search.searcher(ranges, workers)`; the engine drives `EvalClient`.
- Settings: the `candidates`, `climb` and `far` groups in `common/tuning.py`.
- The dashboard draws this schedule: wave 0, generation 1, Step B until minute 55, Step C's generation and fine-tune, and the idle-core GA. v9.18 is v9.16 + v9.17 plus leftover fill (B29) and a far seat that stays far (B30). Its GA child keeps each verified group with one parent: base-quality floor with the gap penalty, the activity threshold, the mapping-quality floor, and the PCR model; the calling threshold with the mismapping rate and the base-quality threshold; the maximum region size with the padding and the minimum region size. Every other setting is chosen on its own, and the one-rung mutation moves only those.

## The full system (built, v9.30)

Built in v9.30. One climber finishes a compass on the highest score. The other cores look for a fingerprint at least 0.10 away, then probe settings the compass has not tried. A core with neither of those waits. The finished rounds say the winner is usually another callset, and the fine-tune adds less than the gap to that winner.

**Two lanes, one pick.** The climb lane is one climber. The hill lane looks only for a fingerprint at least 0.10 away that might beat that climber. The pick is still the highest local score. There is no second selector.

**Work has an order, and a core with nothing in that order waits.**

1. The open compass, then its cube once two faces beat the centre. This batch is small and it finishes.
2. One valley probe, when a fingerprint above the climber has no scored point on the path.
3. Hill-lane children. Each one is a cross of two piles. Fourteen knob pairs in the 748 saved starting configs have Pearson |r| at least 0.20, and those pairs resolve into three groups that a child takes whole from one parent: base-quality floor with the gap penalty, the activity threshold, the mapping-quality floor, and the PCR model; the calling threshold with the mismapping rate and the base-quality threshold; the maximum region size with the padding and the minimum region size. Every other knob was set on its own and comes from either parent. There is no extra mutation. A child must land at least 0.10 in callset distance from the climber and from every hill-lane child still in flight. A child that fails that test is not sent. With only one pile, four to six of the unlinked knobs move by about a quarter of the ladder, and a linked knob stays as the parent set it.
4. Climb-lane probes. Each one copies the centre and changes one to three settings the compass has not tried. They rank the next axes. They do not move the centre, and they do not hold the wave.
5. Stop. Filling the rest with near-duplicates is how the compass stopped finishing. An idle core is the correct state when no child can change the pick by more than the noise.

**Opening.** The portfolio range goes out first and the climber sits as soon as that wave is back. Step A presents the scored configs within 0.50 of the best score. That range stays open through Step B. While the compass pinpoints, one small batch at a time keeps sampling: first a short dummy pass, then the rest of the hypercube. The range is logged again whenever its peak or its size changes. The hill lane and the valley stay. The population mates inside the range.

**The climb, from the end of the opening until eight scaled minutes before the pick.** Up to four axes. An axis already tried near this centre is ordered by the largest one-setting score change seen there. An axis never tried uses the read group, then the assembly group, and a dead axis is skipped. The compass is one rung either side. Sixteen cores stay idle so those heads start together, and the wave closes when the heads return. For the last 12 minutes before the pick, when the climb is the work that can still move the score, each side also tries the next rung. That extra rung does not enter the joint cube. A face inside 0.05 of the centre is not live. Two or more live faces send only their joint configs, and that cube is part of the gate: the wave stays open until those runs return. A crashed run has no score, so that face is not live, and the crash does not end Step B. The centre stays for the whole batch, then moves to the best point that beat it. That score becomes the pile's fingerprint. A pile within 0.10 of the new fingerprint is merged into it, members and all. The next compass starts on the highest score within 0.10. The cores the compass is not holding run the pile GA: half a cross inside the piles, which may land on the climber, and half six unlinked settings aimed at least 0.10 from every pile. A linked group stays with one parent. The compass is not held while that GA runs. A higher score on a new fingerprint takes the valley first. A flat group marks those axes dead until the centre moves, and that move reopens them. A neighbour that tied the centre does not take a slot from a setting with no neighbour yet. Step B stays seated until the fine-tune clock.

**The last eight minutes.** The same compass, on the highest score, still reaching two rungs out. The hill lane stops. A new fingerprint cannot take the seat with too little time left to climb it. A higher score on the same mountain still becomes the next centre once the open grid is back.

**What stays out.** Four climbers. Generation 2. A quiet window that sends no compass. A line, or five rungs, inside the gate. Moving the centre while a grid is in flight. A second copy of the compass that reads config keys as if they were configs.

**How to tell it worked.** Three rounds, TRUE, rank, and the gap to first. On each round the log shows a finished compass and, when two faces won, a finished cube. The hill lane reports its best far fingerprint even when the valley drops it. A climb that gains under 0.10 while the hill lane never produced a fingerprint 0.10 away is a failed opening, and that is the number that matters.

## Evaluation (PROPOSED)

- Replays on several released rounds, v9 vs v4.22. The same config count matters more than the same wall time
  (legacy).
- Run-to-run noise is 0.04–0.08 TRUE, so one pair of runs on one round decides nothing: at least 3 rounds.
- Report: final TRUE, rank, gap to #1, and per step (A: best candidate TRUE; B: gain per candidate; C: gain).
