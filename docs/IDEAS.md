# Ideas backlog

Not planned work - ideas to consider after the fixed plan (evaluator -> answer_set -> simple RF_UCB -> real miner).
Each entry: what, why it might help, what decides whether to build it. Legacy (friend's) ideas are unverified until
we measure them ourselves.

## evaluator

- **E4 threaded variant execution.** One batch does the shared prework (read pass + regions, ~60 s at load), then
  runs each variant's assembly -> PairHMM -> genotyping on a thread pool inside the same JVM. 1 base + 9 assembly
  variants: ~140 s wall instead of ~780 s, same ~780 core-s (small capped batches would cost ~960 core-s).
  Needs: per-thread RNG overrides, per-thread reference reader (and aligner/annotation engine if not thread-safe),
  one ordered task queue per variant with traversal draws synced before each region, more memory (regions held until
  all variants finish). Separate processes/machines are impractical (region state is too big to ship).
  Decide after search exists: worth it only if assembly/activity-heavy families are common.
- **Family cost cap + splitting.** Estimate batch time = solo + sum of variant tier costs; split families that exceed
  ~4-5 min; spread pieces over idle workers. Build together with family-aware EvalClient (step 3).
- **Early divergence abort.** Stop a variant the moment it diverges and queue its solo run immediately instead of at
  batch end. Decide from G2 (how often rounds diverge).
- **Smarter batch main.** Divergence depends on the batch's main config (config b: 0/12 diverged on 055300Z where
  stock had 8/13); the executor could learn groupings that stay exact on the current round.
- **S1 additive scorer.** Per-cluster score contributions so a callset's score is a sum over clusters.
- **E1 fleet executor + E3 node self-test.** One agent per sub-VPS (pull jobs in chunks, heartbeat, requeue on node
  loss), score on the node, nodes download round data directly (main only announces the round).

## answer_set

- **Read-diff against parent BAMs** (legacy, "normal" rounds): round BAM = downsampled parent BAM with edited reads,
  so injections are recoverable by diffing reads. Those rounds have not occurred since 09-12; build only if they
  return. Same trick might work for 1000G rounds if they are built from the public 30x CRAMs.
- **Round-kind detector.** Read length + read-name format tell in seconds which kind a round is.
- **Re-test "the answer set cannot be purified"** (legacy claim: injected and background variants are
  indistinguishable by read statistics). Measure ourselves before relying on it.
- **Upgrade the answer set during a live round** (parked by the user).

## search

- Algorithms reported for winners (only RF_UCB confirmed): QGA, Policy Ladder, QLock, restarts, InjectDrain / GA
  crossover. Legacy found coordinate sweeps ~7x more productive than waves, and that the king's basin can need
  >= 2 axes moved together.
- **Family-aware proposals:** explore genotyping/PairHMM params densely around good points (~1-12% of a run each),
  spend assembly/activity moves sparingly (~55-60%).
- **Portfolio seeds** (legacy): diverse, known-to-run configs from past rounds as starting points.

## main (real miner)

- **Fallback chain** (legacy): exact answer set -> partial -> static config unsearched; never score zero.
- **Early safe submission** of a known-good config, then keep improving; top miners submit ~t+71.6 of 72 min.
- **Round logs**: local score vs real score and rank per round, feeding backtests.
