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

import pandas as pd

from dropseq.aggregation import donor_cell_type_report as report

CELL_TYPE = "model_class_name"

# d1 has 4 cells (one without a cell type), d2 has 2, and two cells have no donor
METADATA = pd.DataFrame({
    "donor": ["d2", "d1", "d1", "d1", "d1", "d2", None, None],
    CELL_TYPE: ["glia", "neuron", "neuron", "glia", None, "glia", "glia", "neuron"],
    "num_retained_transcripts": [10, 100, 300, 50, 7, 30, 1, 1]})


class TestDonorCellTypeReport(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(".tmp", "donor_cell_type_report.")

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def test_report(self):
        df = report.donor_cell_type_report(METADATA, CELL_TYPE)
        self.assertEqual(list(df.columns), report.OUTPUT_COLUMNS)
        self.assertEqual(df.values.tolist(), [
            ["d2", "glia", 2, 1.0, 20.0],
            ["d1", "neuron", 2, 0.5, 200.0],
            ["d1", "glia", 1, 0.25, 50.0]])

    def test_fraction_is_rounded(self):
        metadata = pd.DataFrame({"donor": ["d"] * 3, CELL_TYPE: ["a", "b", "b"],
                                 "num_retained_transcripts": [1, 2, 3]})
        df = report.donor_cell_type_report(metadata, CELL_TYPE)
        self.assertEqual(df["fraction_nuclei"].tolist(), [0.333, 0.667])

    def test_main(self):
        # numeric donors must be handled as labels, and an empty donor or NA is no donor
        metadata = METADATA.assign(donor=["12", "7", "7", "7", "7", "12", "", "NA"])
        input_file = os.path.join(self.tmpDir, "cmd.tsv")
        output_file = os.path.join(self.tmpDir, "report.tsv")
        metadata.to_csv(input_file, sep="\t", index=False)
        self.assertEqual(report.main(["-i", input_file, "-o", output_file, "--cell-type-column", CELL_TYPE]), 0)
        actual = pd.read_csv(output_file, sep="\t", dtype={"donor": str})
        self.assertEqual(actual["donor"].tolist(), ["12", "7", "7"])
        self.assertEqual(actual["num_nuclei"].tolist(), [2, 2, 1])

    def test_main_missing_column(self):
        input_file = os.path.join(self.tmpDir, "cmd.tsv")
        METADATA.to_csv(input_file, sep="\t", index=False)
        self.assertEqual(report.main(["-i", input_file, "-o", os.path.join(self.tmpDir, "r.tsv"),
                                      "--cell-type-column", "nope"]), 1)


if __name__ == "__main__":
    unittest.main()
