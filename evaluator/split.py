"""Split runner: one full-region GATK run as k padded region pieces at once, merged into one callset. Owner: evaluator.

Why: a full-region HaplotypeCaller run is single-threaded where it spends its time (activity detection, assembly), so
extra cores make a step wider, not shorter. Run as k pieces on k workers at once, a config comes back several times
sooner: the search gets more sequential steps per round, paid for in cores (evaluator.split_k; a node then offers
workers // k slots).

Exactness - the random-stream chain. GATK draws from one global random stream through the whole region (QD jitter and,
with contamination_fraction_to_filter > 0, the contamination read removal in genotyping). A piece that simply starts
fresh draws different numbers from its second piece on: 2026-10-21 (test/shard_check, k = 8) 3 of 14 configs with
contamination came out up to 0.036 off, every differing call in pieces 2-8, none in piece 1. Measured since
(test/rng_seek): a piece whose generators are moved to the whole run's state at the start of its part is exact, and
each piece's own draw counts give that state. So every piece runs with the patch's split chain (gatk_patch
HaplotypeCallerEngine, -Dminos.split.*): it walks its region as usual, keeps each own-part region's likelihoods for the
main config and every genotyping / PairHMM batch variant, then waits for the generator states the piece before it hands
on, genotypes its own part again with them (genotyping is ~1% of a run), and hands its states on. Each config's
<output>.split.vcf.gz is exact by construction; the patch checks it (each region's read-pass draws replayed on the
config's own stream, and the region at the cut seen alike by both pieces) and reports "exact" or "diverged". A config
not reported exact by every piece runs once more as one whole run, so a split callset is always the validator's.

Pieces: the region [start, end] is cut into k equal parts; piece i runs on its part plus `pad` bp either side (clipped
to the region, as the validator's GATK sees only the region's reads), and the merge keeps each piece's records whose
POS lies in its own part. Batches: only genotyping and PairHMM variants share a split batch (the chain re-genotypes
them from stored likelihoods); assembly and activity variants run on their own.
"""
from __future__ import annotations

import logging
import time
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from common.types import FULL, CallSet, Config, Round

log = logging.getLogger(__name__)
DEFAULT_PAD = 2000
CHAINABLE = ("genotyping", "pairhmm")      # variant tiers the patch's chain re-genotypes


def pieces(round: Round, k: int, pad: int) -> List[Tuple[str, int, int]]:
    """(region string, own part start, own part end exclusive) of each of k pieces, 1-based like the round."""
    length = round.end - round.start + 1
    out = []
    for i in range(k):
        lo = round.start + i * length // k
        hi = round.start + (i + 1) * length // k
        out.append((f"{round.chrom}:{max(round.start, lo - pad)}-{min(round.end, hi - 1 + pad)}", lo, hi))
    return out


def merge(parts: List[str], cores: List[Tuple[int, int]], out: Path) -> int:
    """Write each part's records whose POS lies in its own core to `out` (bgzipped, indexed); returns the count."""
    import pysam
    out.unlink(missing_ok=True)
    Path(str(out) + ".tbi").unlink(missing_ok=True)
    head = pysam.VariantFile(parts[0])
    n = 0
    with pysam.VariantFile(str(out), "wz", header=head.header) as w:
        for path, (lo, hi) in zip(parts, cores):
            with pysam.VariantFile(path) as f:
                for rec in f:
                    if lo <= rec.pos < hi:
                        w.write(rec)
                        n += 1
    head.close()
    pysam.tabix_index(str(out), preset="vcf", force=True)
    return n


def _unlink(path) -> None:
    if path:
        Path(path).unlink(missing_ok=True)
        Path(str(path) + ".tbi").unlink(missing_ok=True)


def read_status(path: Path) -> Dict[str, Tuple[str, str]]:
    """A piece's chain report: {"main" | "v<n>" (the n-th variant): (its output's file name or "-", "exact" |
    "diverged")}; {} when the chain did not finish (unsupported, cancelled, error, or no file)."""
    try:
        text = path.read_text()
    except OSError:
        return {}
    out = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) == 3:
            out[parts[0]] = (parts[1], parts[2])
    return out


class SplitRunner:
    """Wraps a GatkRunner: FULL runs and batches are split into k chained pieces; every other fidelity (a sub-region,
    WHOLE) and every other method go to the inner runner unchanged."""

    def __init__(self, inner, k: int, pad: int = DEFAULT_PAD):
        self.inner, self.k, self.pad = inner, max(1, int(k)), max(0, int(pad))
        self.stats = {"configs": 0, "exact": 0, "whole_reruns": 0}

    def __getattr__(self, name):
        return getattr(self.inner, name)

    def _split(self, fidelity: str) -> bool:
        return self.k > 1 and fidelity == FULL

    def _tier(self, main: Config, variant: Config, fidelity: str) -> Optional[str]:
        from evaluator.gatk import VARIANT_COST
        cost = self.inner.variant_cost(main, variant, fidelity)
        if cost is None:
            return None
        return next((t for t, c in VARIANT_COST.items() if abs(c - cost) < 1e-9), None) if cost > 0 else "genotyping"

    def variant_cost(self, main: Config, variant: Config, fidelity: str = FULL) -> Optional[float]:
        """Split batches hold only variants the chain re-genotypes; assembly and activity variants are not batchable."""
        cost = self.inner.variant_cost(main, variant, fidelity)
        if cost is None or not self._split(fidelity):
            return cost
        return cost if self._tier(main, variant, fidelity) in CHAINABLE else None

    def _cap(self, round: Round, config: Config, regions: List[str]) -> None:
        """The merged run's haplotype cap: the largest of its pieces' (each reproduces its piece), unknown if any is."""
        caps = [self.inner.kbest_needed(round.round_id, config.key, r) for r in regions]
        with self.inner._warm_lock:
            key = (round.round_id, config.key, FULL)
            if caps and all(c is not None for c in caps):
                self.inner._kbest[key] = max(caps)
            else:
                self.inner._kbest.pop(key, None)

    def _chain(self, round: Round, configs: List[Config]) -> List[CallSet]:
        """configs[0] and its (chainable) variants as k chained pieces; one exact FULL callset per config."""
        from evaluator.gatk import DIVERGED
        ps = pieces(round, self.k, self.pad)
        chain = uuid.uuid4().hex[:12]
        folder = self.inner._out(round, configs[0], ps[0][0]).parent      # the round's callsets: /data/output in GATK
        cancel = folder / f"split_{chain}.cancel"
        status = [folder / f"split_{chain}.p{i}.status" for i in range(len(ps))]
        state = [folder / f"split_{chain}.p{i}.state" for i in range(len(ps))]

        def piece(i: int) -> List[CallSet]:
            region, lo, hi = ps[i]
            props = {"minos.split.part": f"{lo}:{hi}",
                     "minos.split.main.out": f"/data/output/{self.inner._out(round, configs[0], region).name}.split.vcf.gz",
                     "minos.split.status": f"/data/output/{status[i].name}",
                     "minos.split.state.out": f"/data/output/{state[i].name}",
                     "minos.split.cancel": f"/data/output/{cancel.name}"}
            if i > 0:
                props["minos.split.state.in"] = f"/data/output/{state[i - 1].name}"
            self.inner.extra.props = props
            try:
                if len(configs) == 1:
                    return [self.inner.run(round, configs[0], region)]
                return self.inner.run_batch(round, configs, region)
            finally:
                self.inner.extra.props = None

        t0 = time.time()
        results: List[Optional[List[CallSet]]] = [None] * len(ps)
        with ThreadPoolExecutor(len(ps), thread_name_prefix="split") as pool:
            futures = {pool.submit(piece, i): i for i in range(len(ps))}
            for f in as_completed(futures):
                i = futures[f]
                try:
                    results[i] = f.result()
                except Exception as e:  # noqa: BLE001 - a piece that raised is a failed piece
                    results[i] = [CallSet(round.round_id, c.key, "", "", 0.0, ok=False, error=str(e)[-300:],
                                          fidelity=ps[i][0]) for c in configs]
                if any(not cs.ok and cs.error != DIVERGED for cs in results[i]):
                    cancel.write_text("cancelled\n")          # the pieces after it stop waiting for its states
        wall = time.time() - t0
        reports = [read_status(p) for p in status]
        out = []
        for j, config in enumerate(configs):
            self.stats["configs"] += 1
            parts = [results[i][j] for i in range(len(ps))]
            split_vcfs = [str(self.inner._out(round, config, region)) + ".split.vcf.gz" for region, _, _ in ps]
            key = "main" if j == 0 else f"v{j}"
            want = ["-" if j == 0 else self.inner._out(round, config, region).name for region, _, _ in ps]
            crashed = next((i for i, p in enumerate(parts) if not p.ok and p.error != DIVERGED), None)
            exact = crashed is None and all(reports[i].get(key) == (want[i], "exact") for i in range(len(ps)))
            if crashed is not None:
                p = parts[crashed]
                out.append(CallSet(round.round_id, config.key, "", "", wall, ok=False, fidelity=FULL,
                                   error=f"{p.error} (split piece {crashed + 1}/{len(ps)} {ps[crashed][0]})"))
            elif exact:
                target = self.inner._out(round, config, FULL)
                merge(split_vcfs, [(lo, hi) for _, lo, hi in ps], target)
                self._cap(round, config, [r for r, _, _ in ps])
                out.append(self.inner._success(round.round_id, config, target, wall, FULL))
                self.stats["exact"] += 1
            else:
                log.info("split chain of %s not exact (%s): one whole run", config.key,
                         ", ".join(f"p{i + 1} {reports[i].get(key, ('', 'no report'))[1]}" for i in range(len(ps))))
                self.stats["whole_reruns"] += 1
                out.append(self.inner.run(round, config, FULL))
            for p, s in zip(parts, split_vcfs):
                _unlink(p.vcf)
                _unlink(s)
        for f in status + state + [cancel]:
            f.unlink(missing_ok=True)
            Path(str(f) + ".part").unlink(missing_ok=True)
        return out

    def run(self, round: Round, config: Config, fidelity: str = FULL) -> CallSet:
        if not self._split(fidelity):
            return self.inner.run(round, config, fidelity)
        return self._chain(round, [config])[0]

    def run_batch(self, round: Round, configs: List[Config], fidelity: str = FULL, warm: bool = False) -> List[CallSet]:
        """A split batch of the head and its chainable variants; any other variant runs on its own, split."""
        if not self._split(fidelity):
            return self.inner.run_batch(round, configs, fidelity, warm)
        head = configs[0]
        keep = [head] + [c for c in configs[1:] if self._tier(head, c, fidelity) in CHAINABLE]
        by_key = {c.key: cs for c, cs in zip(keep, self._chain(round, keep))}
        return [by_key[c.key] if c.key in by_key else self._chain(round, [c])[0] for c in configs]
