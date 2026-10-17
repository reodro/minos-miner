"""answer_set.k2.kg_origin: the AWS mirror address of a person's CRAM, and which source is tried first."""
import unittest
from unittest import mock

from answer_set.k2 import kg_origin


class MirrorTest(unittest.TestCase):
    def tearDown(self):
        kg_origin.PREFER_MIRROR.discard("NA20887")
        kg_origin.PREFER_MIRROR.discard("HG00405")

    def test_main_set_mirror_and_order(self):
        ebi = "https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240025/NA20887.final.cram"
        with mock.patch.object(kg_origin, "cram_url", return_value=(ebi, kg_origin.INDEXES[0])):
            urls = kg_origin.cram_urls("NA20887")
            self.assertEqual(urls, [ebi, "https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data/"
                                         "ERR3240025/NA20887.final.cram"])
            kg_origin.PREFER_MIRROR.add("NA20887")             # after a failure the mirror goes first
            self.assertEqual(kg_origin.cram_urls("NA20887"), urls[::-1])

    def test_related_set_mirror(self):
        ebi = "https://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3988761/HG00405.final.cram"
        with mock.patch.object(kg_origin, "cram_url", return_value=(ebi, kg_origin.INDEXES[1])):
            self.assertEqual(kg_origin.cram_urls("HG00405")[1],
                             "https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/additional_698_related/"
                             "data/ERR3988761/HG00405.final.cram")

    def test_fetch_falls_back_to_the_mirror(self):
        tried = []

        def fake(url, *a, **k):
            tried.append(url.split("/")[2])
            if "ebi.ac.uk" in url:
                raise ConnectionRefusedError("Connection refused")
            return "local.bam"
        ebi = "https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240025/NA20887.final.cram"
        with mock.patch.object(kg_origin, "cram_url", return_value=(ebi, kg_origin.INDEXES[0])), \
                mock.patch.object(kg_origin, "_fetch_from", side_effect=fake), \
                mock.patch.object(kg_origin, "path", return_value=mock.MagicMock(exists=lambda: False)):
            self.assertEqual(kg_origin.fetch("NA20887", "chr16", 100, 200, budget_s=60), "local.bam")
        self.assertEqual(tried, ["ftp.sra.ebi.ac.uk", "1000genomes.s3.amazonaws.com"])
        self.assertIn("NA20887", kg_origin.PREFER_MIRROR)


if __name__ == "__main__":
    unittest.main()
