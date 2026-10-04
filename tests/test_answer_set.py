import gzip
import os
import unittest

from common.settings import SETTINGS
from common.types import Regime, Round, Variant


class K2OutputTest(unittest.TestCase):
    """K2Provider reads the builder's files as written; scoring is pointed at those exact files."""

    def _write(self, d, name, lines):
        with gzip.open(os.path.join(d, name), "wt") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")
            f.writelines(line + "\n" for line in lines)

    def test_read_output(self):
        import json
        import tempfile
        from pathlib import Path
        from answer_set.k2.provider import read_output
        with tempfile.TemporaryDirectory() as d:
            self._write(d, "truth.vcf.gz", ["chr21\t100\t.\tA\tG\t60\tPASS\t.\tGT\t0/1",
                                            "chr21\t200\t.\tCA\tC\t60\tPASS\tSYNTHETIC\tGT\t1/1",
                                            "chr21\t300\t.\tTAA\tT,TA\t60\tPASS\tSYNTHETIC\tGT\t1/2"])
            self._write(d, "mutations.vcf.gz", [f"chr21\t{p}\t.\tA\tG\t60\tPASS\tSYNTHETIC\tGT\t0/1" for p in (1, 2, 3)])
            with gzip.open(os.path.join(d, "scores.tsv.gz"), "wt") as f:
                f.write("pos\tref\talt\tp\n200\tCA\tC\t0.9100\n300\tTAA\tT\t0.8000\n300\tTA\tT\t0.7000\n")
            Path(d, "report.json").write_text(json.dumps({"mode": "GIAB:HG004", "model": "giab track g006", "threshold": 0.5}))
            r = Round("2026-09-30T03:29:00+00:00", "chr21", 1, 1000, bam="/x.bam", reference="/ref.fa")
            a = read_output(r, Path(d))
            self.assertEqual((a.regime, a.donor, a.n_targets), (Regime.GIAB, "HG004", 3))
            self.assertEqual(a.germline, [Variant("chr21", 100, "A", "G", "0/1")])
            self.assertEqual([(c.variant.pos, c.p) for c in a.synthetic], [(200, 0.91), (300, 0.7)])   # 1/2: weaker allele
            self.assertEqual(a.notes["truth_vcf"], os.path.join(d, "truth.vcf.gz"))
            self.assertEqual((a.notes["region"], a.notes["reference"]), ("chr21:1-1000", "/ref.fa"))
            r2 = Round(r.round_id, "chr21", 1, 1000, bam="/x.bam", reference="/ref.fa", n_targets=611)
            self.assertEqual(read_output(r2, Path(d)).n_targets, 611)                                # platform N wins

    def test_regime_of(self):
        from answer_set.k2.provider import regime_of
        self.assertEqual(regime_of("1000G:HG03559"), (Regime.KG, "HG03559"))
        self.assertEqual(regime_of("1000G:NA18878 WITHOUT DRAGEN (download failed, reads only)"), (Regime.KG, "NA18878"))
        self.assertEqual(regime_of("unrecognised (reads only)"), (Regime.UNKNOWN, None))


class BuilderFactoryTest(unittest.TestCase):
    def test_names(self):
        from answer_set import BUILDER_NAMES, K2Provider, builder
        self.assertIsInstance(builder(), type(builder(BUILDER_NAMES[0])))           # None -> the default
        self.assertIsInstance(builder("k2"), K2Provider)
        with self.assertRaises(ValueError):
            builder("nope")

    def test_env_override(self):
        from unittest import mock
        from answer_set import K2Provider, builder
        from answer_set.builders import SearchableOracle
        with mock.patch.dict(os.environ, {"ANSWER_SET_BUILDER": "oracle"}):
            self.assertIsInstance(builder(), SearchableOracle)
            self.assertIsInstance(builder("k2"), K2Provider)                          # an explicit name still wins
        with mock.patch.dict(os.environ, {"ANSWER_SET_BUILDER": ""}):
            self.assertIsInstance(builder(), K2Provider)                              # empty: the default

    @unittest.skipUnless(os.path.exists(os.path.join(SETTINGS.rounds_dir, "20261001T164100Z", "truth.vcf.gz")),
                         "needs released round 2026-10-01 16:41")
    def test_oracle_is_searchable(self):
        from answer_set.builders import SearchableOracle
        from backtest.data import load
        a = SearchableOracle().build(load("20261001T164100Z"))
        self.assertEqual(a.regime, Regime.GIAB)                                       # main searches only if not UNKNOWN
        self.assertTrue(a.germline and a.synthetic)


class K2RobustnessTest(unittest.TestCase):
    """A weaker build must never be silent: retries, input manifest, warnings."""

    def test_retry(self):
        from answer_set.k2.util import retry
        calls, waits = [], []

        def flaky():
            calls.append(1)
            if len(calls) < 3:
                raise OSError("connection reset")
            return "ok"
        self.assertEqual(retry(flaky, 3, 10, sleep=waits.append), ("ok", ["OSError: connection reset"] * 2))
        self.assertEqual(waits, [10, 20])
        with self.assertRaises(RuntimeError) as e:
            retry(lambda: (_ for _ in ()).throw(OSError("down")), 2, 0, sleep=lambda s: None)
        self.assertIn("failed after 2 attempts: OSError: down", str(e.exception))

    def test_retry_bug_is_not_retried(self):
        from answer_set.k2.util import retry
        calls = []

        def bug():
            calls.append(1)
            return None < 1                                   # the TypeError readdiff hit on an unmapped read
        with self.assertRaises(RuntimeError) as e:
            retry(bug, 3, 10, sleep=lambda s: self.fail("waited before a retry of a bug"))
        self.assertEqual(len(calls), 1)
        self.assertIn("a bug, not retried", str(e.exception))

    def test_retry_budget(self):
        from answer_set.k2.util import retry
        now, waits = [0.0], []

        def slow():
            now[0] += 100                                     # every try takes 100 s and fails
            raise OSError("timeout")
        with self.assertRaises(RuntimeError) as e:
            retry(slow, 5, 10, sleep=lambda s: (waits.append(s), now.__setitem__(0, now[0] + s)), budget_s=250,
                  clock=lambda: now[0])
        self.assertEqual(waits, [10, 20])                     # tries end at 100, 210, 330; a 30 s wait: past 250
        self.assertIn("out of time after 3 attempts", str(e.exception))

    def test_bcftools_timeout_removes_container(self):
        import subprocess
        from unittest import mock
        from answer_set.k2 import util
        calls = []

        def run(args, **kw):
            calls.append(args)
            if args[:2] == ["docker", "run"]:
                raise subprocess.TimeoutExpired(args, kw["timeout"])
            return subprocess.CompletedProcess(args, 0)
        with mock.patch.object(util.subprocess, "run", run):
            with self.assertRaises(RuntimeError) as e:
                util.bcftools("bcftools view x", [], timeout_s=5)
        self.assertIn("timed out after 5 s", str(e.exception))
        name = calls[0][calls[0].index("--name") + 1]
        self.assertEqual(calls[1], ["docker", "rm", "-f", name])

    def test_original_reads_fetch_is_bounded(self):
        """A stuck download fails the fetch within its budget instead of blocking the round."""
        import subprocess
        import tempfile
        from pathlib import Path
        from unittest import mock
        from answer_set.k2 import giab_origin
        with tempfile.TemporaryDirectory() as d:
            Path(d, "HG001.novaseq.pcr-free.50x.dedup.grch38.bam.bai").write_bytes(b"x")
            now = [0.0]

            def run(args, **kw):
                now[0] += kw["timeout"]                       # every piece hangs until killed
                raise subprocess.TimeoutExpired(args, kw["timeout"])
            with mock.patch.object(giab_origin, "CACHE", Path(d)), \
                    mock.patch.object(giab_origin.subprocess, "run", run), \
                    mock.patch.object(giab_origin.time, "monotonic", lambda: now[0]), \
                    mock.patch.object(giab_origin.time, "sleep", lambda s: now.__setitem__(0, now[0] + s)), \
                    mock.patch.object(giab_origin.cf, "ThreadPoolExecutor",
                                      lambda n: __import__("concurrent.futures").futures.ThreadPoolExecutor(1)):
                with self.assertRaises(RuntimeError):
                    giab_origin.fetch("HG001", "chr1", 1000, 2000, chunks=2, budget_s=400)
            self.assertLessEqual(now[0], 400 + giab_origin.PIECE_TIMEOUT_S)
            self.assertEqual([x.name for x in Path(d).iterdir()], ["HG001.novaseq.pcr-free.50x.dedup.grch38.bam.bai"])

    def test_readdiff_skips_unmapped_reads(self):
        """An unmapped round read (placed at its mate, no reference_end) that is an exact piece of an original read
        broke the whole read comparison (TypeError sorting read ends) in 1 of 28 rounds."""
        import random
        import tempfile
        import pysam
        from answer_set.k2 import readdiff
        rnd = random.Random(7)
        seq = "".join(rnd.choice("ACGT") for _ in range(151))
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100000}]}

        def bam(path, reads):
            with pysam.AlignmentFile(path, "wb", header=header) as o:
                for name, s, start, unmapped in reads:
                    r = pysam.AlignedSegment()
                    r.query_name, r.query_sequence, r.reference_id, r.reference_start = name, s, 0, start
                    r.query_qualities = pysam.qualitystring_to_array("I" * len(s))
                    if unmapped:
                        r.flag = 4
                    else:
                        r.cigarstring, r.mapping_quality = f"{len(s)}M", 60
                    o.write(r)
            pysam.index(path)
        with tempfile.TemporaryDirectory() as d:
            bam(f"{d}/orig.bam", [("o1", seq, 5000, False)])
            bam(f"{d}/round.bam", [("r1", seq[10:138], 5010, False), ("r2", seq[20:148], 5020, True)])
            found, stats = readdiff.recover(f"{d}/round.bam", f"{d}/orig.bam", "chr1", 4000, 6000, workers=1)
        self.assertEqual(found, {})
        self.assertEqual(stats, {"unedited": 1})

    def test_manifest(self):
        import tempfile
        from pathlib import Path
        from answer_set.k2.inputs import check, record
        with tempfile.TemporaryDirectory() as d:
            root, man = Path(d), Path(d, "inputs.json")
            a, b = Path(d, "a.vcf.gz"), Path(d, "b.bed")
            a.write_text("AAAA"); b.write_text("BB")
            self.assertEqual(check([a], man, root), [f"no input manifest {man} (python -m answer_set.k2.inputs record)"])
            self.assertEqual(record([a, b], man, root), 2)
            self.assertEqual(check([a, b], man, root, full=True), [])
            a.write_text("AAAAA")                                                     # other size
            os.utime(b, ns=(1, 1))                                                    # touched
            self.assertEqual(check([a, b, Path(d, "c")], man, root),
                             ["changed: a.vcf.gz", "modified since recorded: b.bed", "missing: c"])
            b.write_text("XX"); os.utime(b, ns=(1, 1))                                 # same size, new content
            self.assertEqual(check([b], man, root, full=True), ["changed: b.bed"])

    def test_warnings(self):
        from answer_set.k2.provider import warnings_of
        self.assertEqual(warnings_of({"inputs": "ok", "readdiff": {"edited": 5}}), [])
        w = warnings_of({"inputs": ["changed: giab/HG002.vcf.gz"], "readdiff": {"error": "failed after 3 attempts"}})
        self.assertEqual(len(w), 2)
        self.assertIn("HG002", w[0]); self.assertIn("read comparison failed", w[1])
        self.assertIn("needed 2 tries", warnings_of({"readdiff": {"failed_tries": ["OSError: x"]}})[0])


class K2HapgenoTest(unittest.TestCase):
    """The local haplotype genotyper on a synthetic 1 kb reference with a 12-A homopolymer at 501-512."""

    def setUp(self):
        import random
        import tempfile
        import pysam
        self.tmp = tempfile.TemporaryDirectory()
        d = self.tmp.name
        rnd = random.Random(3)
        body = "".join(rnd.choice("CGT") for _ in range(1000))   # no A: the homopolymer is the only A run
        self.ref = body[:500] + "A" * 12 + body[512:]
        with open(f"{d}/chrT.fa", "w") as f:
            f.write(">chrT\n" + self.ref + "\n")
        pysam.faidx(f"{d}/chrT.fa")
        self.fasta = f"{d}/chrT.fa"
        self.header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrT", "LN": len(self.ref)}]}

    def tearDown(self):
        self.tmp.cleanup()

    def _bam(self, path, reads):
        import pysam
        with pysam.AlignmentFile(path, "wb", header=self.header) as o:
            for name, seq, start in sorted(reads, key=lambda r: r[2]):
                r = pysam.AlignedSegment()
                r.query_name, r.query_sequence, r.reference_id, r.reference_start = name, seq, 0, start
                r.cigarstring, r.mapping_quality = f"{len(seq)}M", 60
                r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
                o.write(r)
        pysam.index(path)

    def _round(self, originals_carry_alt):
        """8 round reads with an extra A in the homopolymer (each cut from its own original) + 8 reference reads."""
        alt = self.ref[:506] + "A" + self.ref[506:]
        rounds, origs = [], []
        for i in range(16):
            o0 = 400 + 3 * i                                   # original: 151 bp; round read: its bases 10..137
            carry = i % 2 == 0
            src = alt if carry and originals_carry_alt else self.ref
            mark = lambda x: x[:20] + ("A" if x[20] != "A" else "C") + x[21:]   # the original's own sequencing error
            o_seq = mark(src[o0:o0 + 151])                      # (outside the window), so each read has one source
            r_src = mark((alt if carry else self.ref)[o0:o0 + 151])
            origs.append((f"o{i}", o_seq, o0))
            rounds.append((f"r{i}", r_src[10:138], o0 + 10))
        d = self.tmp.name
        self._bam(f"{d}/round.bam", rounds)
        self._bam(f"{d}/orig.bam", origs)
        return f"{d}/round.bam", f"{d}/orig.bam"

    def test_minos_style_spelling(self):
        import pysam
        from answer_set.k2 import hapgeno
        fa = pysam.FastaFile(self.fasta)
        s, e = 470, 560
        ref = self.ref[s - 1:e]
        hap = ref[:20] + ref[21:40] + "T" + ref[40:]              # a deletion and an insertion 20 bp apart
        ev = hapgeno._events_of(fa, "chrT", s, ref, hap)
        self.assertEqual(len(ev), 2)
        self.assertTrue(all(len(k[1]) != len(k[2]) for k in ev))  # two indels, not a run of mismatches

    def test_edited_insertion_is_genotyped(self):
        from answer_set.k2 import hapgeno
        bam, orig = self._round(originals_carry_alt=False)
        res = hapgeno.genotype(bam, self.fasta, "chrT", [(501, 502)], threads=1, orig_bam=orig)
        ins = {k: v for k, v in res.items() if k[0] != "window"}
        self.assertEqual(list(ins), [(500, self.ref[499], self.ref[499] + "A")])   # left-aligned 1-bp insertion
        v = ins[(500, self.ref[499], self.ref[499] + "A")]
        self.assertGreater(v["hg_post"], 0.99)
        self.assertLess(v["hg_hom"], 0.01)                        # het: half the reads
        self.assertEqual((v["hg_edit"], v["hg_edit_ori"], v["hg_unedit"]), (8.0, 8.0, 0.0))

    def test_donor_variant_is_not_an_edit(self):
        from answer_set.k2 import hapgeno
        bam, orig = self._round(originals_carry_alt=True)        # the originals already carry the extra A
        res = hapgeno.genotype(bam, self.fasta, "chrT", [(501, 502)], threads=1, orig_bam=orig)
        v = res[(500, self.ref[499], self.ref[499] + "A")]
        self.assertEqual((v["hg_edit"], v["hg_unedit"]), (0.0, 8.0))


class K2HapwinTest(unittest.TestCase):
    """Window decisions: the genotyper's whole-window answer replaces the per-variant one where trusted."""

    def setUp(self):
        import tempfile
        from answer_set.k2.util import Bed
        self.tmp = tempfile.TemporaryDirectory()
        with open(f"{self.tmp.name}/c.bed", "w") as f:
            f.write("chrT\t0\t1000\n")
        self.bed = Bed(f"{self.tmp.name}/c.bed", "chrT")
        self.s, self.ref = 101, "ACGTTGCA" + "A" * 10 + "CGTGCATG"         # window 101-126, 10 A at 109-118
        self.ins = (108, "A", "AA")                                     # one more A (left-aligned)
        alt = self.ref[:8] + "A" + self.ref[8:]
        self.w = {"ref": self.ref, "pair": (self.ref, alt), "post": 0.99, "post2": 0.01, "lr2": 9.0, "hg_anchored": 14.0,
                  "hg_n_haps": 2.0, "hg_win_edit_other": 0.0, "ref_reads": 7,
                  "called": [{"events": [self.ins], "hom": False, "reads": 7, "edit": 7, "edit_ori": 7, "unedit": 0,
                              "edit_q": 36.0}]}

    def tearDown(self):
        self.tmp.cleanup()

    def test_added_where_trusted(self):
        import numpy as np
        from answer_set.k2 import hapwin
        donor_rec = (300, "C", "T")                                     # elsewhere: untouched
        recs = [(300, "C", "T", "hom", False, 2.0)]
        answer = {(r[0], r[1], r[2]): (r[3], r[4]) for r in recs}
        c = hapwin.candidates([(self.s, self.s + len(self.ref) - 1, self.w)], answer, {}, {donor_rec: "hom"}, self.bed)
        self.assertEqual(len(c), 1)
        self.assertEqual(c[0][4], {self.ins: "het"})
        self.assertEqual(c[0][5].shape, (len(hapwin.FEATS),))
        y, ok = hapwin.label(c, {self.ins: ("het", True), donor_rec: ("hom", False)})
        self.assertEqual((list(y), ok), ([1], [False]))                 # called pair right, answer wrong: a fix
        out, n = hapwin.apply(recs, c, np.array([0.9]), 0.4)
        self.assertEqual(n, 1)
        self.assertIn((108, "A", "AA", "het", True, 0.9), out)
        self.assertIn(recs[0], out)
        self.assertEqual(hapwin.apply(recs, c, np.array([0.2]), 0.4), (recs, 0))   # below the cut-off: unchanged

    def test_donor_records_never_removed(self):
        from answer_set.k2 import hapwin
        d = (112, "A", "G")                                             # the person's own record inside the window
        recs = [(112, "A", "G", "het", False, 2.0)]
        answer = {(r[0], r[1], r[2]): (r[3], r[4]) for r in recs}
        c = hapwin.candidates([(self.s, self.s + len(self.ref) - 1, self.w)], answer, {}, {d: "het"}, self.bed)
        self.assertEqual(c, [])                                         # the call lacks it: the window is not eligible

    def test_donor_genotype_kept(self):
        import numpy as np
        from answer_set.k2 import hapwin
        recs = [(108, "A", "AA", "hom", False, 2.0)]                    # published hom; the genotyper says het + more
        w = dict(self.w, called=[{"events": [self.ins, (120, "G", "T")], "hom": False, "reads": 7, "edit": 7,
                                  "edit_ori": 7, "unedit": 0, "edit_q": 36.0}],
                 pair=(self.ref, self.ref[:8] + "A" + self.ref[8:19] + "T" + self.ref[20:]))
        answer = {(r[0], r[1], r[2]): (r[3], r[4]) for r in recs}
        c = hapwin.candidates([(self.s, self.s + len(self.ref) - 1, w)], answer, {}, {self.ins: "hom"}, self.bed)
        out, _ = hapwin.apply(recs, c, np.array([0.9]), 0.4)
        self.assertIn(recs[0], out)                                     # still hom, still the person's record
        self.assertEqual(sum(1 for r in out if r[:3] == self.ins), 1)


class K2SpellingTest(unittest.TestCase):
    """Spelling of replaced windows: options include Minos's insertion + deletion pairs; the chooser's pick is written."""

    def setUp(self):
        import random
        import tempfile
        import pysam
        self.tmp = tempfile.TemporaryDirectory()
        rnd = random.Random(5)
        body = "".join(rnd.choice("CG") for _ in range(400))
        self.ref_full = body[:200] + "ATATATATATATATATAT" + body[218:]     # an AT repeat at 201-218
        with open(f"{self.tmp.name}/chrS.fa", "w") as f:
            f.write(">chrS\n" + self.ref_full + "\n")
        pysam.faidx(f"{self.tmp.name}/chrS.fa")
        self.fa = pysam.FastaFile(f"{self.tmp.name}/chrS.fa")

    def tearDown(self):
        self.tmp.cleanup()

    def test_pair_spellings_find_minos_pair(self):
        from answer_set.k2 import hapwin
        s = 190
        ref = self.ref_full[s - 1:s + 39]
        pair = [(195, ref[5], ref[5] + "T"), (210, "TA", "T")]        # an insertion and a deletion 15 bp apart
        hap = hapwin._apply(ref, s, pair)
        from answer_set.k2 import hapcand
        got = [sorted(x) for x in hapwin.pair_spellings(self.fa, "chrS", s, ref, hap)]
        self.assertTrue(got)
        self.assertTrue(all(hapwin._apply(ref, s, g) == hap for g in got))          # every option: the same DNA
        self.assertIn(sorted(hapcand.left_align(self.fa, "chrS", *k) for k in pair), got)   # Minos's own pair

    def test_overlapping_spellings_are_dropped(self):
        from answer_set.k2 import hapwin
        s, ref = 190, self.ref_full[189:229]
        hap = hapwin._apply(ref, s, [(195, ref[5], ref[5] + "T"), (210, "TA", "T")])
        for ev, _ in hapwin.spelling_options(self.fa, "chrS", s, ref, hap):
            ev = sorted(ev)
            self.assertTrue(all(y[0] > x[0] + len(x[1]) - 1 for x, y in zip(ev, ev[1:])))   # no overlapping records
            self.assertEqual(hapwin._apply(ref, s, ev), hap)                                # the called DNA
        self.assertFalse(hapwin._clean([(200, "TAT", "T"), (201, "A", "AG")], s, ref, hap))  # insertion inside a deletion

    def test_respell_writes_the_chosen_spelling(self):
        import numpy as np
        from answer_set.k2 import hapwin

        class PreferPairs:                                            # a chooser that always picks a pair spelling
            def predict_proba(self, X):
                q = (X[:, hapwin.SPELL_FEATS.index("scheme")] == 99).astype(float)
                return np.c_[1 - q, q]
        s = 190
        ref = self.ref_full[s - 1:s + 39]
        hap = hapwin._apply(ref, s, [(195, ref[5], ref[5] + "T"), (210, "TA", "T")])
        ev = hapwin.respell(self.fa, "chrS", s, {"ref": ref, "pair": (hap, hap)}, {}, PreferPairs())
        self.assertEqual(len(ev), 2)
        self.assertEqual(set(ev.values()), {"hom"})                   # both copies carry it
        self.assertEqual(hapwin._apply(ref, s, list(ev)), hap)        # same DNA as called

    def test_own_records_keep_their_spelling(self):
        import numpy as np
        from answer_set.k2 import hapwin
        own = (195, "C", "CT")
        recs = [(195, "C", "CT", "het", False, 2.0)]
        cand = (190, 230, {}, [(own, "het")], {own: "het"}, np.zeros(len(hapwin.FEATS)))
        out, n = hapwin.apply(recs, [cand], np.array([0.9]), 0.4, respell_fn=lambda c: {(201, "A", "G"): "het"})
        self.assertEqual((out, n), (recs, 1))                         # not respelled, the person's record kept


class K2WriteTest(unittest.TestCase):
    """The person's own records are written exactly as Minos's truth writes them: NIST's REF, ALT order and genotype
    order (2|1 -> 2/1). The validator classes a two-allele record by the allele its genotype names first."""

    def test_donor_genotype_order_kept(self):
        import gzip
        import tempfile
        from unittest import mock
        import pysam
        from answer_set.k2 import build
        with tempfile.TemporaryDirectory() as d:
            ref = "A" * 100 + "CGTGTGTGTGTGTGT" + "A" * 100
            with open(f"{d}/chrW.fa", "w") as f:
                f.write(">chrW\n" + ref + "\n")
            pysam.faidx(f"{d}/chrW.fa")
            raw = {101: ("CGTGTGTGTGTGTGT", ("CGT", "C"), "2/1")}            # NIST 2|1: a 12-bp and a 14-bp deletion
            from answer_set.k2.util import trim
            recs = [trim(101, "CGTGTGTGTGTGTGT", a_) + ("het", False) for a_ in ("CGT", "C")]   # as read_alleles gives them
            recs.append((50, "A", "G", "het", True))                             # + an inserted SNP
            with mock.patch.object(build, "ref_path", lambda c: __import__("pathlib").Path(f"{d}/chrW.fa")):
                build.write_vcf(f"{d}/t.vcf.gz", "chrW", recs, raw)
            lines = [l.rstrip("\n").split("\t") for l in gzip.open(f"{d}/t.vcf.gz", "rt") if not l.startswith("#")]
        self.assertEqual([l[1] for l in lines], ["50", "101"])                # sorted
        self.assertEqual(lines[1][3:5] + [lines[1][9]], ["CGTGTGTGTGTGTGT", "CGT,C", "2/1"])
        self.assertEqual(lines[0][7], "SYNTHETIC")


class K2AssembleTest(unittest.TestCase):
    """Model scores -> truth records: donor background copied, threshold, <= 2 alleles per position, BED rule."""

    def test_giab_round(self):
        from answer_set.k2.assemble import assemble_records
        donor = {(100, "A", "G"): "het", (150, "C", "CT"): "hom"}
        keys = [(100, "A", "G"), (200, "T", "C"), (300, "G", "A"), (300, "G", "GA"), (300, "G", "T"),
                (400, "A", "C"), (500, "C", "G")]
        p_any = [0.9, 0.95, 0.9, 0.8, 0.7, 0.3, 0.99]
        p_hom = [0.0, 0.9, 0.9, 0.1, 0.1, 0.0, 0.0]
        recs = assemble_records(keys, p_any, p_hom, 0.5, True, donor, allowed=lambda k: k[0] < 500)
        got = {r[:5] for r in recs}
        self.assertEqual(got, {(100, "A", "G", "het", False), (150, "C", "CT", "hom", False),    # donor, as published
                               (200, "T", "C", "hom", True),                                     # model call, P(hom) 0.9
                               (300, "G", "A", "het", True), (300, "G", "GA", "het", True)})     # top 2 of 3 -> 1/2
        # 400: below threshold; 500: outside the donor's confident region

    def test_kg_round(self):
        from answer_set.k2.assemble import assemble_records
        recs = assemble_records([(100, "A", "G")], [0.9], [0.2], 0.5, False, {})
        self.assertEqual([r[:5] for r in recs], [(100, "A", "G", "het", False)])        # 1000G: nothing SYNTHETIC


ROUND_1641 = os.path.join(SETTINGS.rounds_dir, "20261001T164100Z")


@unittest.skipUnless(os.environ.get("RUN_SLOW") and os.path.exists(os.path.join(ROUND_1641, "input.bam")),
                     "set RUN_SLOW=1 (needs round 2026-10-01 16:41, the answer_set/k2 data, Docker, internet; ~1 min)")
class K2EndToEndTest(unittest.TestCase):
    def test_giab_round_matches_released_truth(self):
        import tempfile
        from answer_set import K2Provider
        from answer_set.k2.util import read_alleles, ref_path
        from backtest.data import load
        r = load("20261001T164100Z", reference=str(ref_path("chr18")))
        with tempfile.TemporaryDirectory() as d:
            a = K2Provider(work_dir=d).build(r)
            ours = {k for k, _, _ in read_alleles(a.notes["truth_vcf"])}
            real = {k: "SYNTHETIC" in i for k, _, i in read_alleles(os.path.join(r.released_dir, "truth.vcf.gz"),
                                                                     keep_info=True)}
        self.assertEqual((a.regime, a.donor), (Regime.GIAB, "HG002"))
        self.assertNotIn("error", a.notes["report"]["readdiff"])                   # Google originals reached
        self.assertEqual(a.notes["warnings"], [])                                  # inputs as recorded, nothing degraded
        self.assertEqual({k for k, syn in real.items() if not syn} - ours, set())  # donor background copied exactly
        m = len(ours & real.keys())
        p, rc = m / len(ours), m / len(real)
        self.assertGreaterEqual(2 * p * rc / (p + rc), 0.998)                      # measured 0.9991


if __name__ == "__main__":
    unittest.main()
