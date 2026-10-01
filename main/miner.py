"""The miner process: watch for a round, fetch it, build the answer set, search, submit. Owner: main.

    run_forever(source, make_parts, reference_for)   live: poll the platform, mine every new round (needs a wallet)
    mine_round(round, parts, ...)     one round with the fallback chain, used by live and replay

Fallback chain (a round must never end without a valid submission):
  1. a known-good baseline is submitted as soon as it has run (main.loop does this first)
  2. answer set failed, or no donor identified -> the baseline alone, unsearched (search without a truth is noise)
  3. otherwise search until the deadline, resubmitting at most every `resubmit_interval_s` (chain commitments are
     rate-limited), and the best config found is submitted once more just before the loop stops.

Live platform access (UpstreamPlatform) follows the upstream miner (minos_subnet/neurons/miner.py): signed
round-status and submit-config calls, and the on-chain config commitment when the network config enables it.
It NEVER pays a resubmission fee: if a submission would need payment it is not sent. It needs bittensor + httpx
and a wallet, none of which exist on the dev box yet - so it is untested; test it in demo mode first.
"""
from __future__ import annotations

import asyncio
import json
import logging
import os
import subprocess
import threading
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Callable, List, Optional, Protocol, Tuple

from answer_set.base import AnswerSetProvider
from common.params import ParamRanges
from common.settings import SETTINGS
from common.store import ResultStore
from common.types import Config, Regime, Round
from evaluator.base import Executor, Scorer
from main.loop import LoopSettings, RoundOutcome, run_round_async
from main.submitter import Submitter
from search.base import Strategy

log = logging.getLogger(__name__)


# ---- the pieces one round needs ------------------------------------------------------------------------------
@dataclass
class Parts:
    provider: AnswerSetProvider
    executor: Executor
    scorer: Scorer
    strategies: List[Strategy]
    submitter: Submitter
    ranges: ParamRanges
    baselines: List[Config]                       # known-good configs, best first; [0] is the static fallback
    settings: LoopSettings = field(default_factory=LoopSettings)
    store: Optional[ResultStore] = None
    status: Optional[object] = None               # main.status.StatusBoard: what the dashboard shows


@dataclass
class MinedRound:
    round_id: str
    mode: str                                     # "searched" | "static" (fallback) | "failed"
    outcome: Optional[RoundOutcome] = None
    note: str = ""


def mine_round(round: Round, parts: Parts, deadline_ts: Optional[float], bam_url: Optional[str] = None,
               bai_url: Optional[str] = None) -> MinedRound:
    if parts.status:
        parts.status.round_started(round, deadline_ts)
    mined = _mine_round(round, parts, deadline_ts, bam_url, bai_url)
    if parts.status:
        parts.status.round_ended(mined)
    return mined


def _mine_round(round: Round, parts: Parts, deadline_ts: Optional[float], bam_url: Optional[str] = None,
                bai_url: Optional[str] = None) -> MinedRound:
    """One round with the fallback chain. Never raises. With a fleet executor, the round (truth files, BAM location)
    is published to the nodes once the answer set exists."""
    try:
        answer_set = parts.provider.build(round)
    except Exception as e:  # noqa: BLE001 - the round still gets the static baseline
        log.exception("answer set failed for %s", round.round_id)
        return _static(round, parts, f"answer set failed: {e}")
    if parts.status:
        parts.status.answer_set(answer_set)
    if answer_set.regime == Regime.UNKNOWN or not answer_set.germline:
        return _static(round, parts, "no donor identified: searching without a truth would be noise")
    publish = getattr(parts.executor, "publish_round", None)
    if publish:
        try:
            publish(round, answer_set, bam_url, bai_url)
        except Exception as e:  # noqa: BLE001 - no fleet for this round: still submit something
            log.exception("publishing %s to the fleet failed", round.round_id)
            return _static(round, parts, f"fleet publish failed: {e}")

    class Built:                                  # the loop builds once; hand it the answer set we already have
        def build(self, _round):
            return answer_set

        def update(self, a, callsets):
            return parts.provider.update(a, callsets)

    try:
        out = asyncio.run(run_round_async(round, Built(), parts.executor, parts.scorer, parts.strategies,
                                          parts.submitter, parts.ranges, parts.settings, parts.baselines,
                                          store=parts.store, deadline_ts=deadline_ts,
                                          on_tick=parts.status.tick if parts.status else None))
        return MinedRound(round.round_id, "searched", out)
    except Exception as e:  # noqa: BLE001
        log.exception("round loop failed for %s", round.round_id)
        if _submitted(parts, round):              # a better config already went out: the baseline must not replace it
            return MinedRound(round.round_id, "searched", note=f"round loop failed after submitting: {e}")
        return _static(round, parts, f"round loop failed: {e}")


def _submitted(parts: Parts, round: Round) -> bool:
    return any(rid == round.round_id for rid, _ in getattr(parts.submitter, "sent", []))


def _static(round: Round, parts: Parts, why: str) -> MinedRound:
    log.warning("%s: static fallback (%s)", round.round_id, why)
    for c in parts.baselines:
        if not parts.ranges.errors(c) and parts.submitter.submit(round, c):
            return MinedRound(round.round_id, "static", note=why)
    return MinedRound(round.round_id, "failed", note=why + "; no baseline could be submitted")


# ---- live: the platform ------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class LiveRound:
    round_id: str
    region: str
    deadline_ts: Optional[float]                  # submission_end_time as epoch seconds
    bam_url: str
    bai_url: str
    num_mutations: Optional[int] = None           # N, if the platform ever reports it live (Optional upstream)


class RoundSource(Protocol):
    def current(self) -> Optional[LiveRound]:
        """The round open for submissions right now, or None."""


def fetch_round(live: LiveRound, reference_for: Callable[[str], str], rounds_dir: Optional[str] = None) -> Round:
    """Download the round's BAM + index (presigned URLs) and describe it as a Round."""
    chrom, start, end = Round.parse_region(live.region)
    tag = Round.tag_of(live.round_id)
    d = Path(rounds_dir or SETTINGS.rounds_dir) / tag
    d.mkdir(parents=True, exist_ok=True)
    for url, name in ((live.bam_url, "input.bam"), (live.bai_url, "input.bam.bai")):
        dest = d / name
        if not dest.exists():
            subprocess.run(["curl", "-sSfL", "--retry", "3", "-o", str(dest) + ".part", url], check=True, timeout=900)
            os.replace(str(dest) + ".part", dest)
    (d / "round.json").write_text(json.dumps({"round_id": live.round_id, "region": live.region}))
    deadline = datetime.utcfromtimestamp(live.deadline_ts).isoformat() + "Z" if live.deadline_ts else None
    return Round(live.round_id, chrom, start, end, bam=str(d / "input.bam"), reference=reference_for(chrom),
                 deadline=deadline, n_targets=live.num_mutations)


def run_forever(source: RoundSource, make_parts: Callable[[Round], Parts], reference_for: Callable[[str], str],
                poll_s: float = 20.0) -> None:
    """Live loop: every new open round is fetched and mined once. If the BAM cannot be fetched, the round still gets
    the static baseline (submitting needs no BAM)."""
    done = set()
    while True:
        try:
            live = source.current()
        except Exception as e:  # noqa: BLE001 - the platform being down is not a reason to stop
            log.warning("round status failed: %s", e)
            live = None
        if live and live.round_id not in done:
            done.add(live.round_id)
            try:
                rnd = fetch_round(live, reference_for)
            except Exception as e:  # noqa: BLE001
                log.exception("fetching round %s failed", live.round_id)
                chrom, start, end = Round.parse_region(live.region)
                rnd = Round(live.round_id, chrom, start, end, bam="", reference=reference_for(chrom))
                mined = _static(rnd, make_parts(rnd), f"round download failed: {e}")
            else:
                mined = mine_round(rnd, make_parts(rnd), live.deadline_ts, live.bam_url, live.bai_url)
            log.info("round %s: %s %s", mined.round_id, mined.mode, mined.note)
        time.sleep(poll_s)


class UpstreamPlatform:
    """Round status + submission through the upstream client (untested: needs bittensor, httpx and a wallet).

    Mirrors neurons/miner.py: signed POST /v2/round-status and /v2/submit-config; when network config says
    config_commitment_enabled, a salted commitment is recorded locally and published on chain first (rate-limited
    by the chain - a skipped commitment still submits, like upstream). Never pays: if resubmission fees are enabled
    and the free allowance is used, the submission is not sent.
    """

    def __init__(self, wallet_name: str, hotkey: str, netuid: int = 107, network: str = "finney", demo: bool = False):
        import sys
        sys.path.insert(0, SETTINGS.minos_subnet_dir)
        import bittensor as bt  # noqa: F401 - imported here so the rest of the project runs without it
        from utils.platform_client import MinerPlatformClient, PlatformConfig
        self.bt = bt
        self.wallet = bt.wallet(name=wallet_name, hotkey=hotkey)
        self.subtensor = None if demo else bt.subtensor(network=network)
        self.netuid, self.demo = netuid, demo
        self.client = MinerPlatformClient(self.wallet.hotkey, PlatformConfig(base_url=SETTINGS.platform_url), demo=demo)
        self.submits: dict = {}
        # the upstream client is async and submit() is called from inside main's running event loop, so every
        # platform call runs on this object's own loop in a helper thread
        self._loop = asyncio.new_event_loop()
        threading.Thread(target=self._loop.run_forever, daemon=True, name="platform-loop").start()

    def _call(self, coro, timeout: float = 120.0):
        return asyncio.run_coroutine_threadsafe(coro, self._loop).result(timeout)

    def current(self) -> Optional[LiveRound]:
        s = self._call(self.client.get_round_status())
        if not s.get("has_active_round") or s.get("status") != "open":
            return None
        end = s.get("submission_end_time")
        deadline = datetime.fromisoformat(str(end).replace("Z", "+00:00")).timestamp() if end else None
        return LiveRound(s["round_id"], s["region"], deadline, s["bam_presigned_url"], s["bam_index_presigned_url"],
                         s.get("num_mutations"))

    def submit(self, round: Round, config: Config) -> bool:
        from utils import config_commit
        net = self._call(self.client.get_network_config()) or {}
        used = self.submits.get(round.round_id, 0)
        free = net.get("free_submissions_per_round", 1)                # upstream submission_payment's key
        if net.get("resubmission_fee_enabled") is True and used >= int(1 if free is None else free):
            log.warning("%s: resubmission would need a fee; not submitting (fees are never paid)", round.round_id)
            return False
        tool_config = dict(config.to_submission())        # {"tool", "<tool>_options"} like upstream, plus:
        try:
            from utils.config_loader import get_tool_version   # upstream helper, as its miner uses
            tool_config["version"] = get_tool_version(config.tool)
        except Exception as e:  # noqa: BLE001
            log.warning("tool version unavailable (%s); submitting without it", e)
        commitment = block = nonce = None
        if net.get("config_commitment_enabled") is True and not self.demo:
            commitment, block, nonce = self._commit(round.round_id, config.tool, tool_config, config_commit)
        res = self._call(self.client.submit_config(round_id=round.round_id, tool_name=config.tool,
                                                   tool_config=tool_config, config_commitment=commitment,
                                                   commitment_block=block, config_nonce=nonce))
        ok = bool(res.get("success"))
        if ok:
            self.submits[round.round_id] = used + 1
        return ok

    def _commit(self, round_id, tool, tool_config, config_commit):
        from utils import bt_compat
        hotkey = self.wallet.hotkey.ss58_address
        nonce = config_commit.new_nonce()
        commitment = config_commit.compute_commitment(netuid=self.netuid, round_id=round_id, hotkey=hotkey,
                                                      tool_name=tool, tool_config=tool_config, nonce=nonce)
        # the nonce is persisted BEFORE publishing: a commitment that cannot be opened is worse than none
        config_commit.CommitmentLedger().record({"round_id": round_id, "hotkey": hotkey, "netuid": self.netuid,
                                                 "tool_name": tool, "nonce": nonce, "commitment": commitment,
                                                 "block": None, "tool_config": config_commit.submission_config(tool_config)})
        block = None
        try:
            ok, reason, block = bt_compat.commit(self.subtensor, wallet=self.wallet, netuid=self.netuid,
                                                 data=config_commit.chain_payload(round_id, commitment))
            if not ok:
                log.warning("%s: chain commitment skipped (%s); submitting anyway", round_id, reason)
        except Exception as e:  # noqa: BLE001 - like upstream, a failed commitment does not block the submission
            log.warning("%s: chain commitment failed (%s)", round_id, e)
        return commitment, block, nonce


class PlatformSubmitter:
    """main.submitter.Submitter over a live platform (UpstreamPlatform or anything with .submit)."""

    def __init__(self, platform):
        self.platform = platform
        self.sent: List[Tuple[str, str]] = []             # (round_id, config_key) accepted, like DryRunSubmitter

    def submit(self, round: Round, config: Config) -> bool:
        try:
            ok = self.platform.submit(round, config)
            if ok:
                self.sent.append((round.round_id, config.key))
            return ok
        except Exception as e:  # noqa: BLE001 - a failed submission is reported, not raised into the loop
            log.warning("submission failed for %s: %s", round.round_id, e)
            return False
