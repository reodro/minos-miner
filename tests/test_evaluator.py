import os
import unittest
from pathlib import Path
from unittest import mock

from common.settings import SETTINGS
from common.types import FULL, Config, Round
from evaluator import SAFE_ANNOTATION_EXCLUSIONS, GatkRunner
from tests.helpers import ranges

ROUND_DIR = os.path.join(SETTINGS.rounds_dir, "20260928T032900Z")


class GatkCommandTest(unittest.TestCase):
    """The runner's command must be the validator template's command plus only the verified -AX exclusions."""

    def test_command_is_validator_command_plus_safe_exclusions(self):
        import sys
        sys.path.insert(0, SETTINGS.minos_subnet_dir)
        from templates import gatk as template
        r = Round("r", "chr20", 1, 1000, bam="/x/input.bam", reference="/y/chr20.fa")
        cfg = ranges().defaults("gatk")
        ours = GatkRunner(threads=1, memory_gb=4).command(r, cfg, FULL, Path("/z/out.vcf.gz"))

        captured = {}
        def fake_run(cmd, **kw):
            captured["cmd"] = cmd
            raise FileNotFoundError      # stop right after the command is built
        with mock.patch.object(template.subprocess, "run", side_effect=fake_run), \
                mock.patch.object(template.Path, "exists", return_value=True), \
                mock.patch.object(template.os, "cpu_count", return_value=16):
            template.variant_call(Path("/x/input.bam"), Path("/y/chr20.fa"), Path("/z/out.vcf.gz"), r.region,
                                  {"gatk_options": cfg.params, "threads": 1, "memory_gb": 4})
        theirs = captured["cmd"]
        strip = lambda c: [x for i, x in enumerate(c) if not (i > 0 and c[i - 1] == "--name")]  # noqa: E731
        extra = []
        for a in SAFE_ANNOTATION_EXCLUSIONS:
            extra += ["-AX", a]
        self.assertEqual(strip(ours), strip(theirs) + extra)

    def test_safe_set_keeps_quality_affecting_annotations(self):
        for kept in ("QualByDepth", "DepthPerAlleleBySample", "RMSMappingQuality"):
            self.assertNotIn(kept, SAFE_ANNOTATION_EXCLUSIONS)

    def test_rejects_other_tools(self):
        with self.assertRaises(ValueError):
            GatkRunner().run(Round("r", "chr20", 1, 2, "", ""), Config("bcftools", {}), FULL)


@unittest.skipUnless(os.environ.get("RUN_SLOW") and os.path.exists(os.path.join(ROUND_DIR, "input.bam")),
                     "set RUN_SLOW=1 (needs round 20260928T032900Z, reference, docker; ~2 min)")
class GatkEndToEndTest(unittest.TestCase):
    def test_stock_config_reproduces_published_validator_score(self):
        from answer_set import OracleProvider
        from backtest.data import load
        from evaluator import HappyScorer
        r = load("20260928T032900Z", reference=os.path.join(SETTINGS.data_dir, "reference/chr22/chr22.fa"))
        cs = GatkRunner(threads=4).run(r, ranges().defaults("gatk"))
        res = HappyScorer().score(cs, OracleProvider().build(r))
        self.assertAlmostEqual(res.score / 100, 0.7665917881, places=10)   # a published validator score


if __name__ == "__main__":
    unittest.main()
