"""main/answer_reuse.py: a finished k2 build is read again after a restart; a partial one is built again."""
import tempfile
import unittest
import unittest.mock
from pathlib import Path
from types import SimpleNamespace

from main.answer_reuse import DONE_FILES, reusing


class K2Provider:                                     # the name the wrapper checks; builds count calls
    def __init__(self, work_dir):
        self.work_dir, self.builds = work_dir, 0

    def build(self, round):
        self.builds += 1
        return "built"

    def update(self, a, cs):
        return a


class ReuseTest(unittest.TestCase):
    def test_finished_build_is_read_and_partial_build_is_rebuilt(self):
        with tempfile.TemporaryDirectory() as d:
            inner = K2Provider(d)
            b = reusing(inner)
            rnd = SimpleNamespace(tag="20261021T022600Z")
            out = Path(d) / rnd.tag / "v1"
            out.mkdir(parents=True)
            (out / "work").mkdir()                    # a killed build: only work/
            self.assertEqual(b.build(rnd), "built")
            self.assertEqual(inner.builds, 1)
            for f in DONE_FILES:
                (out / f).write_text("x")
            with unittest.mock.patch("answer_set.k2.provider.read_output", return_value="read") as ro:
                self.assertEqual(b.build(rnd), "read")
                ro.assert_called_once()
            self.assertEqual(inner.builds, 1)         # no second build


if __name__ == "__main__":
    unittest.main()
