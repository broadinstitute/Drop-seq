# MIT License
#
# Copyright 2026 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import os
import shutil
import tempfile
import unittest

from dropseq.mmc import convert_mmc_csv

DATA_DIR = "tests/data/mmc"
MODEL_PREFIX = "HMBA_Human_WB_v0.5_"


class TestConvertMmcCsv(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(".tmp", "convert_mmc_csv.")
        self.csv = os.path.join(DATA_DIR, "mmc.csv")

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def test_matches_zamboni_output(self):
        """The output must be identical to the tsv written by convertPrefixMMCCsvClp for the same csv."""
        output = os.path.join(self.tmpDir, "out.tsv")
        self.assertEqual(convert_mmc_csv.main(["--input", self.csv, "--output", output,
                                               "--column-prefix", MODEL_PREFIX]), 0)
        with open(output) as actual, open(os.path.join(DATA_DIR, "mmc.prefixed.tsv")) as expected:
            self.assertEqual(actual.read(), expected.read())

    def test_no_prefix(self):
        df = convert_mmc_csv.convert_mmc_csv(self.csv)
        self.assertEqual(df.columns[0], "cell_id")
        self.assertEqual(df.columns[1], "neighborhood_label")

    def test_prefix_skips_cell_id(self):
        df = convert_mmc_csv.convert_mmc_csv(self.csv, MODEL_PREFIX)
        self.assertEqual(df.columns[0], "cell_id")
        self.assertTrue(all(c.startswith(MODEL_PREFIX) for c in df.columns[1:]))

    def test_read_taxonomy_hierarchy(self):
        self.assertEqual(convert_mmc_csv.read_taxonomy_hierarchy(self.csv),
                         ["neighborhood", "class", "subclass", "cluster"])

    def test_taxonomy_hierarchy_missing(self):
        path = os.path.join(self.tmpDir, "no_hierarchy.csv")
        with open(path, "w") as f:
            f.write("# metadata = x.json\ncell_id,a\nAAA,1\n")
        with self.assertRaises(ValueError):
            convert_mmc_csv.read_taxonomy_hierarchy(path)

    def _convert(self, *extra):
        column_file = os.path.join(self.tmpDir, "column.txt")
        rc = convert_mmc_csv.main(["--input", self.csv, "--output", os.path.join(self.tmpDir, "out.tsv"),
                                   "--column-prefix", MODEL_PREFIX, "--cell-type-column-output", column_file, *extra])
        return rc, column_file

    def test_cell_type_column_default_level(self):
        rc, column_file = self._convert()
        self.assertEqual(rc, 0)
        with open(column_file) as f:
            self.assertEqual(f.read(), MODEL_PREFIX + "neighborhood_name\n")

    def test_cell_type_column_level(self):
        rc, column_file = self._convert("--level", "class")
        self.assertEqual(rc, 0)
        with open(column_file) as f:
            self.assertEqual(f.read(), MODEL_PREFIX + "class_name\n")

    def test_cell_type_column_unknown_level(self):
        rc, column_file = self._convert("--level", "nonsense")
        self.assertEqual(rc, 1)
        self.assertFalse(os.path.exists(column_file))


if __name__ == "__main__":
    unittest.main()
