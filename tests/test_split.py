"""evaluator/split.py: a full-region run as k padded region pieces, merged; and main's whole-run check before it
submits a config with contamination (main/loop.py)."""
import asyncio
import os
import random
import tempfile
import threading
import unittest
import unittest.mock
from pathlib import Path

import pysam

from common.tuning import TUNING
from common.types import FULL, WHOLE, AnswerSet, CallSet, Config, Regime, Round, Variant
from evaluator.split import SplitRunner, merge, pieces

ROUND = Round("2026-01-01T00:00:00+00:00", "chr1", 1001, 2000, "", "")
CONTIG = "chr1"


def write_vcf(path: Path, positions) -> str:
    header = pysam.VariantHeader()
    header.contigs.add(CONTIG, length=10_000)
    header.add_sample("S")
    header.formats.add("GT", 1, "String", "Genotype")
    with pysam.VariantFile(str(path), "wz", header=header) as w:
        for pos in sorted(positions):
            rec = w.new_record(contig=CONTIG, start=pos - 1, stop=pos, alleles=("A", "C"))
            rec.samples["S"]["GT"] = (0, 1)
            w.write(rec)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    return str(path)


def positions_of(vcf: str):
    with pysam.VariantFile(vcf) as f:
        return [r.pos for r in f]


class FakeInner:
    """What SplitRunner needs from a GatkRunner, with the patch's split chain: a run writes the 'true' calls inside its
    region (plus one padding artefact at each inner cut, which the merge must drop); with chain properties it also
    writes <output>.split.vcf.gz and its status line ("exact" unless the config is in `diverge`)."""

    def __init__(self, root: Path, calls, caps=None, fail_region=None, diverge=None, no_report=None):
        self.root, self.calls, self.caps = root, sorted(calls), caps or {}
        self.fail_region, self.diverge, self.no_report = fail_region, diverge or set(), no_report
        self._warm_lock, self._kbest, self.runs = threading.Lock(), {}, []
        self.extra = threading.local()

    def _region(self, round, fidelity):
        if fidelity in (FULL, WHOLE):
            return round.start, round.end, False
        a, b = fidelity.split(":")[1].split("-")
        return int(a), int(b), True

    def _out(self, round, config, fidelity):
        out = self.root / round.tag / f"{config.key}_{abs(hash(fidelity)) % 10**8}.vcf.gz"
        out.parent.mkdir(parents=True, exist_ok=True)
        return out

    def _host(self, container_path):
        return self.root / ROUND.tag / Path(container_path).name

    def _walk(self, round, config, fidelity):
        lo, hi, piece = self._region(round, fidelity)
        pos = [p for p in self.calls if lo <= p <= hi]
        if piece:
            pos += [p for p, cut in ((lo + 1, lo > round.start), (hi - 1, hi < round.end)) if cut and p not in pos]
        out = self._out(round, config, fidelity)
        vcf = write_vcf(out, pos)
        props = getattr(self.extra, "props", None)
        if props:
            write_vcf(Path(str(out) + ".split.vcf.gz"), pos)
        cap = self.caps.get(fidelity)
        if cap is not None:
            with self._warm_lock:
                self._kbest[(round.round_id, config.key, fidelity)] = cap
        return vcf

    def _report(self, configs, fidelity, names):
        props = getattr(self.extra, "props", None)
        if not props or fidelity == self.no_report:
            return
        lines = [f"{'main' if i == 0 else f'v{i}'}\t{n}\t{'diverged' if c.key in self.diverge else 'exact'}"
                 for i, (c, n) in enumerate(zip(configs, names))]
        self._host(props["minos.split.status"]).write_text("\n".join(lines) + "\n")

    def run(self, round, config, fidelity=FULL):
        self.runs.append((config.key, fidelity))
        if fidelity == self.fail_region:
            return CallSet(round.round_id, config.key, "", "", 1.0, ok=False, error="rc=3 Graph must have ref source",
                           fidelity=fidelity)
        vcf = self._walk(round, config, fidelity)
        self._report([config], fidelity, ["-"])
        return CallSet(round.round_id, config.key, vcf, "x", 1.0, fidelity=fidelity)

    def run_batch(self, round, configs, fidelity=FULL, warm=False):
        out = []
        for c in configs:
            self.runs.append((c.key, fidelity))
            out.append(CallSet(round.round_id, c.key, self._walk(round, c, fidelity), "x", 1.0, fidelity=fidelity))
        self._report(configs, fidelity, ["-"] + [self._out(round, c, fidelity).name for c in configs[1:]])
        return out

    def variant_cost(self, main, variant, fidelity=FULL):
        if "assembly" in variant.params:
            return 0.55
        return 0.011

    def kbest_needed(self, round_id, key, fidelity=FULL):
        with self._warm_lock:
            return self._kbest.get((round_id, key, fidelity))

    def _success(self, round_id, config, out, dt, fidelity):
        return CallSet(round_id, config.key, str(out), "x", dt, fidelity=fidelity,
                       haplotype_cap=self.kbest_needed(round_id, config.key, fidelity))

    def branch_key(self, config, fidelity=FULL):
        return "skeleton:" + fidelity


def cfg(i):
    return Config("gatk", {"heterozygosity": 0.001 * (i + 1)})


class PiecesTest(unittest.TestCase):
    def test_cores_cover_the_region_once_and_pieces_are_padded_and_clipped(self):
        ps = pieces(ROUND, 4, 50)
        self.assertEqual([(lo, hi) for _, lo, hi in ps], [(1001, 1251), (1251, 1501), (1501, 1751), (1751, 2001)])
        self.assertEqual([r for r, _, _ in ps],
                         ["chr1:1001-1300", "chr1:1201-1550", "chr1:1451-1800", "chr1:1701-2000"])


class MergeTest(unittest.TestCase):
    def test_each_piece_keeps_only_its_own_core(self):
        d = Path(tempfile.mkdtemp())
        a = write_vcf(d / "a.vcf.gz", [1005, 1240, 1260])         # 1260 lies in piece b's core
        b = write_vcf(d / "b.vcf.gz", [1245, 1260, 1499, 1600])   # 1245 is a's; 1600 lies outside both cores
        n = merge([a, b], [(1001, 1251), (1251, 1501)], d / "m.vcf.gz")
        self.assertEqual(n, 4)
        self.assertEqual(positions_of(str(d / "m.vcf.gz")), [1005, 1240, 1260, 1499])
        self.assertTrue((d / "m.vcf.gz.tbi").exists())


class SplitRunnerTest(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.calls = [1010, 1249, 1250, 1251, 1252, 1600, 1999]

    def left(self, keep=()):
        return sorted(p.name for p in (self.root / ROUND.tag).iterdir() if str(p) not in keep and not p.name.endswith(".tbi"))

    def test_the_chained_pieces_merge_to_the_whole_run_and_leave_nothing_behind(self):
        inner = FakeInner(self.root, self.calls)
        runner = SplitRunner(inner, 4, 50)
        cs = runner.run(ROUND, cfg(0))
        self.assertTrue(cs.ok)
        self.assertEqual(cs.fidelity, FULL)
        self.assertEqual(positions_of(cs.vcf), self.calls)
        self.assertEqual(sorted(f for _, f in inner.runs), sorted(r for r, _, _ in pieces(ROUND, 4, 50)))
        self.assertEqual(self.left(keep=(cs.vcf,)), [])
        self.assertEqual(runner.stats, {"configs": 1, "exact": 1, "whole_reruns": 0})

    def test_each_piece_gets_its_chain_and_only_later_pieces_wait_for_a_state(self):
        inner = FakeInner(self.root, self.calls)
        seen = []
        run = inner.run
        inner.run = lambda r, c, f=FULL: (seen.append((f, dict(inner.extra.props or {}))), run(r, c, f))[1]
        SplitRunner(inner, 3, 50).run(ROUND, cfg(0))
        props = {f: p for f, p in seen}
        ps = pieces(ROUND, 3, 50)
        self.assertNotIn("minos.split.state.in", props[ps[0][0]])
        self.assertEqual(props[ps[1][0]]["minos.split.state.in"], props[ps[0][0]]["minos.split.state.out"])
        self.assertEqual(props[ps[2][0]]["minos.split.part"], f"{ps[2][1]}:{ps[2][2]}")
        self.assertIsNone(getattr(inner.extra, "props", None))

    def test_a_config_not_reported_exact_runs_once_whole(self):
        for kw in ({"diverge": {cfg(0).key}}, {"no_report": pieces(ROUND, 4, 50)[2][0]}):
            inner = FakeInner(self.root, self.calls, **kw)
            runner = SplitRunner(inner, 4, 50)
            cs = runner.run(ROUND, cfg(0))
            self.assertTrue(cs.ok)
            self.assertEqual(inner.runs[-1], (cfg(0).key, FULL))
            self.assertEqual(runner.stats["whole_reruns"], 1)

    def test_whole_and_sub_region_runs_are_not_split(self):
        inner = FakeInner(self.root, self.calls)
        runner = SplitRunner(inner, 4, 50)
        whole = runner.run(ROUND, cfg(0), WHOLE)
        self.assertEqual((whole.fidelity, positions_of(whole.vcf)), (WHOLE, self.calls))
        runner.run(ROUND, cfg(1), "chr1:1001-1100")
        self.assertEqual(inner.runs, [(cfg(0).key, WHOLE), (cfg(1).key, "chr1:1001-1100")])

    def test_off_with_one_piece(self):
        inner = FakeInner(self.root, self.calls)
        SplitRunner(inner, 1).run(ROUND, cfg(0))
        self.assertEqual(inner.runs, [(cfg(0).key, FULL)])

    def test_a_crashed_piece_fails_the_config_says_which_and_cancels_the_chain(self):
        ps = pieces(ROUND, 4, 50)
        inner = FakeInner(self.root, self.calls, fail_region=ps[2][0])
        cs = SplitRunner(inner, 4, 50).run(ROUND, cfg(0))
        self.assertFalse(cs.ok)
        self.assertIn("Graph must have ref source", cs.error)
        self.assertIn("split piece 3/4", cs.error)
        self.assertEqual(self.left(), [])

    def test_the_haplotype_cap_is_the_largest_pieces_or_unknown(self):
        ps = pieces(ROUND, 2, 50)
        inner = FakeInner(self.root, self.calls, caps={ps[0][0]: 6, ps[1][0]: 9})
        self.assertEqual(SplitRunner(inner, 2, 50).run(ROUND, cfg(0)).haplotype_cap, 9)
        inner = FakeInner(self.root, self.calls, caps={ps[0][0]: 6})
        self.assertIsNone(SplitRunner(inner, 2, 50).run(ROUND, cfg(1)).haplotype_cap)

    def test_a_batch_chains_genotyping_variants_and_runs_assembly_variants_alone(self):
        asm = Config("gatk", {"heterozygosity": 0.002, "assembly": 1})
        inner = FakeInner(self.root, self.calls, diverge={cfg(2).key})
        runner = SplitRunner(inner, 3, 50)
        out = runner.run_batch(ROUND, [cfg(0), cfg(1), asm, cfg(2)])
        self.assertEqual([c.config_key for c in out], [cfg(0).key, cfg(1).key, asm.key, cfg(2).key])
        for c in out:
            self.assertTrue(c.ok)
            self.assertEqual(positions_of(c.vcf), self.calls)
        self.assertEqual(sum(1 for k, f in inner.runs if k == asm.key), 3)        # its own 3 pieces
        self.assertIn((cfg(2).key, FULL), inner.runs)                             # diverged: one whole run
        self.assertIsNone(runner.variant_cost(cfg(0), asm))
        self.assertEqual(runner.variant_cost(cfg(0), cfg(1)), 0.011)

    def test_other_methods_reach_the_inner_runner(self):
        self.assertEqual(SplitRunner(FakeInner(self.root, []), 4).branch_key(cfg(0)), "skeleton:full")


class WholeCheckBeforeSubmitTest(unittest.TestCase):
    """With split runs on, a config with contamination is submitted only on its WHOLE run's score."""

    def test_a_contaminated_config_waits_for_its_whole_run(self):
        from evaluator import FakeRunner, FakeScorer, SimExecutor
        from main import DryRunSubmitter, LoopSettings, run_round
        from search.axes import latin_hypercube
        from search.space import Space
        from tests.helpers import ranges
        r = ranges()
        runner = FakeRunner(r, seed=0)
        base = r.defaults("gatk")
        pool = [c for c in latin_hypercube(Space(r, "gatk"), random.Random(1), 400)
                if float(c.params.get("contamination_fraction_to_filter") or 0) > 0]
        good = max(pool, key=runner.quality)
        self.assertGreater(runner.quality(good), runner.quality(base))
        seen = []

        async def searcher(client):
            seen.append(client)
            await client.propose(good)
            for _ in range(200):
                if (good.key, WHOLE) in client.latest:
                    break
                await asyncio.sleep(0)
            await asyncio.sleep(0)

        class Stub:
            def build(self, rnd):
                return AnswerSet(rnd.round_id, 1, Regime.GIAB, "HG002", [Variant("chr20", 5, "A", "G", "0/1")], [],
                                 600, {})

            def update(self, a, c):
                return a

        sub = DryRunSubmitter()
        rnd = Round("2026-01-01T00:00:00+00:00", "chr20", 1, 1000, "", "")
        with unittest.mock.patch.dict(TUNING.values, {"evaluator.split_k": 8}):
            out = run_round(rnd, Stub(), SimExecutor(runner, 8, duration=lambda job: 30.0), FakeScorer(runner),
                            searcher, sub, r, LoopSettings(stop_before_deadline_s=60, resubmit_interval_s=0),
                            baselines=[base], deadline_ts=30 * 60)
        client = seen[0]
        self.assertIn((good.key, WHOLE), client.latest)
        self.assertIn(good.key, [k for _, k in sub.sent])
        when = {k: (t, s) for t, k, s in out.submissions}
        self.assertEqual(when[good.key][1], client.latest[(good.key, WHOLE)].score)


if __name__ == "__main__":
    unittest.main()
