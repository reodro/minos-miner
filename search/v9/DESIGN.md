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

**A3. GA, 3 generations. AGREED.**
- Generations 1–2 breed **within** each group. Generation 3 crosses Core × FP and Both × the others.
- **Diverse parents in every generation.** Legacy crossover collapsed from 792 new children to 123 when later
  generations bred only from the top.
- Mutation: **30% of children get 2 settings redrawn over the full range** (AGREED 2026-10-08;
  the current `genetic.MUTATE_FRAC` / `MUTATE_AXES`, also legacy pipeline v2's rule).
- Measured: v4's three generations gained +0.22 / +0.10 / +0.20 local and beat extra random sampling by about 0.75.

**A3b. Every generation runs (2026-10-08).** A generation launches once 80% of the previous one is back, and the
candidates are picked only once 80% of the last one is back. `candidates.max_min` (30 minutes, scaled) is a safety
limit, not the schedule. Why: on `20261008T175300Z` real GATK waves took 4.5–10 minutes (a 14.1 Mb region); the old
15-minute limit picked the candidates 15 seconds after generation 2 launched and skipped generation 3.

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
- **Cheap riders on every grid point,** batched into the same GATK run: the 4 genotyping settings (~1% of a run
  each) and the 4 PairHMM settings (~12% each) at ±δ. About +1 run per grid point. Gives every expensive × cheap pair.
- **Cheap pairs on the centre's own run:** all pairs of the 8 cheap settings at ±δ (28 pairs × 8 ≈ 224 riders,
  about 15 runs). Covers cheap × cheap synergies (legacy: calling threshold × mismapping rate).
- **Cost per iteration:** about 500 run-equivalents on R2 (243 × ~2), about 160 on R1.

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

**B3. Step growth (across iterations).**
- A setting whose move just won uses **2δ** in the same direction in the next grid, then **3δ, 4δ, … for as long
  as it keeps winning**. The opposite direction stays at δ.
- When a grown step stops winning, that setting returns to δ (fine resolution near the peak is kept).
- Example: A+1 wins the grid over (A, B, C); the next grid is around (A+1, B, C) with A at −1, 0, +2 and B, C at
  −1, 0, +1.

**B4. Randomness. AGREED: 20 random jumps.** When the last read-pass iteration and the last assembly iteration have
both gained less than 0.005 (a stall), 20 random local jumps from the centre: 2–4 settings, each 2–5 rungs (a bool
flips). A jump gaining at least 0.005 becomes the centre and the grid resumes; otherwise the centre is the peak.

**B5. Stop.** When an iteration's best gain is **< 0.005 local** (the resolution at which our answer set still
orders configs) for both grid kinds in a row and the jumps find nothing. A smaller gain still moves the centre. The
peak then goes to Step C.

**B7. Waiting (built 2026-10-08, `climb.wait_riders`, default 0).** A wave decides once its grid points are back;
riders still running count at the climber's next decision. Waiting for every rider made each grid wave 5.8 minutes
instead of 3.0 in the simulator (500 workers): the riders that do not fit a small batch go to big batches of up to
2.3 runs. Same climb, 9 iterations instead of 6 in Step B.

**B6. Time.** One iteration = one grid wave, about one run time plus stragglers.
A candidate still gaining at about minute 60 hands its best to Step C.

---

## Step C: what Step B's grid did not look at. AGREED 2026-10-08, built

Runs on each climber as soon as Step B leaves it at its peak (or at minute 50, `climb.end_min`), while other
climbers may still be climbing (changed 2026-10-08: before, every peak waited for the last climber). Step B's grid only ever tries a setting at −δ and +δ (or the grown step) from the centre, so along each
single setting it leaves two kinds of untried values:

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
- The dashboard still shows v4 panels; it is updated later (not part of this step).

## Evaluation (PROPOSED)

- Replays on several released rounds, v9 vs v4.22. The same config count matters more than the same wall time
  (legacy).
- Run-to-run noise is 0.04–0.08 TRUE, so one pair of runs on one round decides nothing: at least 3 rounds.
- Report: final TRUE, rank, gap to #1, and per step (A: best candidate TRUE; B: gain per candidate; C: gain).
