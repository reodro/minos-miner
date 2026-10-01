# main

The miner process: watch for a round, fetch it, build the answer set, run the search, keep the best config
submitted. Wires the other parts together through their public APIs.

Public API: `mine_round`, `run_forever`, `Parts`, `MinedRound`, `UpstreamPlatform`, `PlatformSubmitter`,
`DryRunSubmitter`, `run_round`/`run_round_async`, `LoopSettings`.

## One round (`miner.mine_round`)

1. **Answer set** (`provider.build`). If it fails, or no donor is identified -> submit the first valid **baseline**
   and stop (searching without a truth is noise). A round never ends without a valid submission.
2. **Search** (`loop.run_round_async`): baselines are evaluated first and the first one that runs is submitted;
   then all strategies run concurrently until `deadline - stop_before_deadline_s`.
3. **Submissions**: a better config is resubmitted at most every `resubmit_interval_s` (300 s live - on-chain
   commitments are rate-limited); the best config found is submitted once more just before the loop stops.
4. Scoring runs in a background thread in packed batches, so GATK workers stay fed while hap.py runs.
5. Round end: after the final submission the executor ends the round (`end_round`): its queued jobs are dropped and
   its running GATK containers killed, locally or on every fleet node, so the next round starts with every worker free.

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
```

With `--fleet`, `mine_round` publishes each round (truth files, BAM location) to the nodes once the answer set
exists; when main stops, every node agent stops. Replay cannot print the TRUE score then (the final VCF is on a
node).

## Open work

- Test the live path (demo mode, then one real round with a known-good static config).
- Round logs: local vs true score and rank per round, feeding backtests.
- Submission timing: tune `resubmit_interval_s` / `stop_before_deadline_s` against the chain's commitment limit.

## Ideas (not planned)

- **Portfolio of known-good configs** as baselines: a third of stored configs crashed on unseen regions in the
  legacy project; the first submission must be one that runs.
- Top miners submit around t+71.6 of 72 minutes (legacy measurement): use the whole window, safely.
