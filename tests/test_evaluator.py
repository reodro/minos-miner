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
        # container name and network are docker settings, not GATK arguments
        strip = lambda c: [x for i, x in enumerate(c) if x != "--network"  # noqa: E731
                           and not (i > 0 and c[i - 1] in ("--name", "--network"))]
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
        r = load("20260928T032900Z", reference=SETTINGS.reference_fasta("chr22"))
        cs = GatkRunner(threads=4).run(r, ranges().defaults("gatk"))
        res = HappyScorer().score(cs, OracleProvider().build(r))
        self.assertAlmostEqual(res.score / 100, 0.7665917881, places=10)   # a published validator score


class VariantFileTest(unittest.TestCase):
    """The patch reads variant columns by position: 1-4 genotyping, 5-8 PairHMM, 9-14 assembly, 15-17 activity."""

    def test_columns_and_defaults(self):
        from evaluator.gatk import variant_line
        f = variant_line(Config("gatk", {"min_pruning": 3, "recover_all_dangling_branches": True,
                                         "max_assembly_region_size": 450}), "/data/output/v.vcf.gz").rstrip("\n").split("\t")
        self.assertEqual(len(f), 18)
        self.assertEqual(f[0], "/data/output/v.vcf.gz")
        self.assertEqual(f[1:5], ["30.0", "0.000125", "6", "0.0"])
        self.assertEqual(f[5:9], ["18", "CONSERVATIVE", "10", "45"])
        self.assertEqual(f[9:15], ["3", "4", "true", "128", "0.001", "2.302585"])
        self.assertEqual(f[15:18], ["0.002", "50", "450"])

    def test_activity_params_batch_together_but_padding_does_not(self):
        r = GatkRunner(patch_jar="/nonexistent.jar")
        base = Config("gatk", {"min_pruning": 2})
        self.assertEqual(r.branch_key(base), r.branch_key(Config("gatk", {"min_pruning": 2, "active_probability_threshold": 0.005,
                                                                         "min_assembly_region_size": 30})))
        self.assertNotEqual(r.branch_key(base), r.branch_key(Config("gatk", {"min_pruning": 2, "assembly_region_padding": 50})))


class WorkersRuleTest(unittest.TestCase):
    def test_gatk_workers_floor_95_percent(self):
        from common.settings import default_workers
        with mock.patch.dict(os.environ, {"MINOS_WORKERS": ""}):
            self.assertEqual([default_workers(c) for c in (2, 8, 16, 32, 36, 180)], [1, 7, 15, 30, 34, 171])


class VariantCostTest(unittest.TestCase):
    def test_most_expensive_changed_stage_sets_the_cost(self):
        from evaluator.gatk import VARIANT_COST
        r = GatkRunner(patch_jar="/nonexistent.jar")
        base = Config("gatk", {"min_pruning": 2})
        self.assertEqual(r.variant_cost(base, base), 0.0)
        gt = Config("gatk", {"min_pruning": 2, "standard_min_confidence_threshold_for_calling": 40.0})
        self.assertEqual(r.variant_cost(base, gt), VARIANT_COST["genotyping"])
        mixed = Config("gatk", {"min_pruning": 3, "pair_hmm_gap_continuation_penalty": 12})
        self.assertEqual(r.variant_cost(base, mixed), max(VARIANT_COST["assembly"], VARIANT_COST["pairhmm"]))
        self.assertEqual(r.variant_cost(base, Config("gatk", {"min_pruning": 2, "active_probability_threshold": 0.004})),
                         VARIANT_COST["activity"])
        self.assertIsNone(r.variant_cost(base, Config("gatk", {"min_pruning": 2, "assembly_region_padding": 50})))
        self.assertIsNone(GatkRunner().variant_cost(base, gt))          # no patch jar: nothing batches


class SpliceRepairTest(unittest.TestCase):
    """Repaired callset = repair run's records (computed regions) + batch records of the regions it skipped."""

    def test_keeps_batch_records_only_in_skipped_regions(self):
        import gzip, tempfile
        from evaluator.gatk import splice_repair
        hdr = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        rec = lambda pos, tag: f"chr1\t{pos}\t.\tA\tG\t50\t.\t{tag}\n"
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            with gzip.open(d / "batch.vcf.gz", "wt") as f:
                f.write(hdr + rec(15, "batch") + rec(150, "batch-wrong") + rec(420, "batch"))
            with gzip.open(d / "repair.vcf.gz", "wt") as f:
                f.write(hdr + rec(160, "repair"))
            (d / "regions").write_text("chr1:1-100|true|3|a\t0\t0\t0\t0\nchr1:101-300|true|9|b\t0\t1\t0\t0\n"
                                       "chr1:301-500|true|4|c\t1\t0\t0\t0\n")
            (d / "repair.log").write_text("S\tchr1:1-100|true|3|a\nC\tchr1:101-300|true|9|z\nS\tchr1:301-500|true|4|c\n")
            splice_repair(d / "batch.vcf.gz", d / "regions", d / "repair.vcf.gz", d / "repair.log", d / "out.vcf.gz")
            with gzip.open(d / "out.vcf.gz", "rt") as f:
                body = [l.split("\t")[1] + ":" + l.split("\t")[7].strip() for l in f if not l.startswith("#")]
            self.assertEqual(body, ["15:batch", "160:repair", "420:batch"])
            self.assertTrue((d / "out.vcf.gz.tbi").exists())


class GvcfNotBatchedTest(unittest.TestCase):
    def test_gvcf_configs_run_alone(self):
        r = GatkRunner(patch_jar="/nonexistent.jar")
        self.assertIsNone(r.branch_key(Config("gatk", {"emit_ref_confidence": "GVCF"})))
        self.assertIsNotNone(r.branch_key(Config("gatk", {"emit_ref_confidence": "NONE"})))


if __name__ == "__main__":
    unittest.main()
