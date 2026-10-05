# main

The miner process: watch for a round, fetch it, run the search while the answer set is built, keep the best config
submitted. Wires the other parts together through their public APIs.

Main builds nothing of another part: it asks `search.searcher(...)`, `evaluator.local_executor(...)` /
`evaluator.fleet_executor(...)` / `evaluator.scorer()` and `answer_set.builder()` for the pieces of a round (`make_parts`).

Public API: `mine_round`, `run_forever`, `Parts`, `MinedRound`, `UpstreamPlatform`, `PlatformSubmitter`,
`DryRunSubmitter`, `run_round`/`run_round_async`, `LoopSettings`.

## One round (`miner.mine_round`)

1. **Publish** (fleet only): the round (BAM location, region) goes to the nodes before anything else; they fetch the
   BAM and run the first wave while main works.
2. **Search** (`loop.run_round_async`) starts at once: the baselines (highest priority, so they run first) and the
   search's opening wave fill every worker from the first minute. Baselines are tried in order until one is
   submitted (or a better config is), so a round never ends empty. The search runs until
   `deadline - stop_before_deadline_s`.
3. **Answer set** (`provider.build`) is built in a thread meanwhile (minutes on a live round); results made before it
   exists wait and are scored the moment it does - on main and on the nodes (`publish_answer_set`, again on every
   newer version). If it fails, or no donor is identified -> the round's work is dropped and the first valid
   **baseline** is submitted, unsearched (searching without a truth is noise). A round never ends without a valid
   submission.
4. **Submissions**: a better config is resubmitted at most every `resubmit_interval_s` (300 s live - on-chain
   commitments are rate-limited); the best config found is submitted once more just before the loop stops.
5. Scoring (one upstream hap.py run per config) runs off the loop in `evaluator.ScoringPool`'s processes - one per
   GATK worker, as a fleet node scores its own results - so GATK workers stay fed while hap.py runs (with a fleet,
   the nodes score and main scores nothing).
6. Last minutes: at round start main tells the executor when the final submission happens
   (`set_finish_by(deadline - stop_before_deadline_s - score_margin_s)`), so it starts no GATK work that could not be
   scored by then.
7. Round end: after the final submission the executor ends the round (`end_round`): its queued jobs are dropped and
   its running GATK containers killed, locally or on every fleet node, so the next round starts with every worker free.

## Status and post-round review (`main.status`, `main.review`)

- **Status** (`StatusBoard`): every few seconds `<status_dir>/live.json`; at round end one `round` line in
  `<status_dir>/rounds.jsonl` with three sections (`answer_set`, `evaluator`, `search`), and `<status_dir>/runs/<run_id>/`
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

**Setup** (main only): `.venv/bin/pip install -r requirements-live.txt` (bittensor 10.3.1 + httpx, upstream's pins;
fleet nodes never need them). Wallets are machine-wide in `~/.bittensor/wallets/`, managed with `btcli` (installed
machine-wide with pipx, `bittensor-cli==9.21.1`). The miner needs the hotkey and the coldkey's public part
(`coldkeypub.txt`) only - the coldkey itself stays offline; rewards go to the coldkey that registered the hotkey.

**Untested on the platform**: the wallet and libraries are installed, but no live round has run yet.
Try `python -m main live --wallet W --hotkey H --demo` first.

## Commands

```bash
.venv/bin/python -m main replay <round_tag> --minutes 20                    # runs the search part's searcher, as live does
.venv/bin/python -m main live --wallet W --hotkey H [--demo]
# either one with a fleet: main only proposes, sub-VPS nodes run GATK and score (see evaluator/README.md)
.venv/bin/python -m main replay <tag> --fleet            # sub-VPSs from fleet.json (repo folder), controlled over ssh
.venv/bin/python -m main review <run_id|tag>             # review one stored run now (answer-set accuracy, TRUE, rank)
```

With `--fleet`, `mine_round` publishes each round (BAM location, region) to the nodes first and the loop publishes
the answer set when it exists; when main stops, every node agent stops. The review then re-runs the configs it
scores on main (their VCFs are on the nodes).

## Open work

- Test the live path (demo mode, then one real round with a known-good static config). The 5-minute review
  watcher only runs inside `live`, so it has only been exercised through `replay` and `review` so far.
- Submission timing: tune `resubmit_interval_s` / `stop_before_deadline_s` against the chain's commitment limit.

## Ideas (not planned)

- **Portfolio of known-good configs** as baselines: a third of stored configs crashed on unseen regions in the
  legacy project; the first submission must be one that runs.
- Top miners submit around t+71.6 of 72 minutes (legacy measurement): use the whole window, safely.
