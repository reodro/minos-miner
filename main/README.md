# main

The miner process: watch for a round, fetch it, build the answer set, run the search, keep the best config
submitted. Wires the other parts together through their public APIs.

Public API: `mine_round`, `run_forever`, `Parts`, `MinedRound`, `UpstreamPlatform`, `PlatformSubmitter`,
`DryRunSubmitter`, `run_round`/`run_round_async`, `LoopSettings`.

## One round (`miner.mine_round`)

1. **Answer set** (`provider.build`). If it fails, or no donor is identified -> submit the first valid **baseline**
   and stop (searching without a truth is noise). A round never ends without a valid submission.
2. **Search** (`loop.run_round_async`): the baselines (highest priority, so they run first) and all strategies start
   together, so every worker is busy from the first minute; baselines are tried in order until one is submitted (or a
   better config is), so a round never ends empty. Strategies run until `deadline - stop_before_deadline_s`.
3. **Submissions**: a better config is resubmitted at most every `resubmit_interval_s` (300 s live - on-chain
   commitments are rate-limited); the best config found is submitted once more just before the loop stops.
4. Scoring (one upstream hap.py run per config) runs in a background thread, so GATK workers stay fed while
   hap.py runs (with a
   fleet, the nodes score).
5. Round end: after the final submission the executor ends the round (`end_round`): its queued jobs are dropped and
   its running GATK containers killed, locally or on every fleet node, so the next round starts with every worker free.

## Status and post-round review (`main.status`, `main.review`)

- **Status** (`StatusBoard`): every few seconds `<data_dir>/status/live.json`; at round end one `round` line in
  `status/rounds.jsonl` with three sections (`answer_set`, `evaluator`, `search`), and `status/runs/<run_id>/`
  (`manifest.json`: final + baseline config, every config's params and local score; `answer_set.json.gz`).
  A run is one mining of one round: `<tag>-<start epoch>` (replays repeat tags). The dashboard reads all of this.
- **Review** (`Reviewer`): next to `live`, a background thread checks every 5 min:
  - truth released (one round after the deadline) -> `review` line: answer-set accuracy vs the released truth
    (germline exact/sites, synthetic, N), TRUE scores of the final, the baseline and the top 8 by local score,
    Spearman(local, true), pick gap, gain over baseline. VCFs on main are reused; others (fleet runs) are re-run
    here with GATK, solo and in parallel (one worker per config, no batching: a batch runs its variants one after
    another, ~2.5-3x a solo run) - that is what main's spare CPU is for.
  - round finalized (a few rounds later) -> `rank` line: rank among the real submissions, gap to the field's best.
  - `replay` reviews at once (its truth is on disk); `python -m main review <run_id|tag>` reviews one run by hand.

## Live platform (`miner.UpstreamPlatform`)

Follows the upstream miner (`minos_subnet/neurons/miner.py`):
- signed `POST /v2/round-status` and `/v2/submit-config` (tool config = `{"tool", "version", "gatk_options"}`);
- when network config says `config_commitment_enabled` (it does today): a salted commitment is recorded in the
  local ledger **before** publishing it on chain, then submitted with the config;
- **never pays a resubmission fee**: if fees are enabled and the free allowance is used, the submission is not sent.

**Untested**: needs `bittensor`, `httpx` and a registered wallet, none of which exist on the dev box yet.
Try `python -m main live --wallet W --hotkey H --demo` first.

## Commands

```bash
.venv/bin/python -m main replay <round_tag> --minutes 20 [--strategy rf_ucb|random]
.venv/bin/python -m main live --wallet W --hotkey H [--demo]
# either one with a fleet: main only proposes, sub-VPS nodes run GATK and score (see evaluator/README.md)
.venv/bin/python -m main replay <tag> --fleet            # sub-VPSs from fleet.json (repo folder), controlled over ssh
.venv/bin/python -m main review <run_id|tag>             # review one stored run now (answer-set accuracy, TRUE, rank)
```

With `--fleet`, `mine_round` publishes each round (truth files, BAM location) to the nodes once the answer set
exists; when main stops, every node agent stops. The review then re-runs the configs it scores on main (their
VCFs are on the nodes).

## Open work

- Test the live path (demo mode, then one real round with a known-good static config). The 5-minute review
  watcher only runs inside `live`, so it has only been exercised through `replay` and `review` so far.
- Submission timing: tune `resubmit_interval_s` / `stop_before_deadline_s` against the chain's commitment limit.

## Ideas (not planned)

- **Portfolio of known-good configs** as baselines: a third of stored configs crashed on unseen regions in the
  legacy project; the first submission must be one that runs.
- Top miners submit around t+71.6 of 72 minutes (legacy measurement): use the whole window, safely.
