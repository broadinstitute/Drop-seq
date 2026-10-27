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

import numpy as np
import pandas as pd

from dropseq.aggregation.plots import qc_tables as qc

CT = "cell_type"

# lib1: d1 (3 cells), d2 (1 cell), unexpected d9 (1 cell); lib2: d1 (1 cell), a cell without a donor
METADATA = pd.DataFrame({
    "PREFIX": ["lib1", "lib1", "lib1", "lib1", "lib1", "lib2", "lib2"],
    "donor": ["d1", "d1", "d1", "d2", "d9", "d1", None],
    CT: ["spn", "spn", "glia", "glia", "spn", None, "spn"],
    "num_retained_transcripts": [1000, 3000, 500, 200, 800, 2000, 50]})


class TestQcTables(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(".tmp", "qc_tables.")

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def write(self, name, text):
        path = os.path.join(self.tmpDir, name)
        with open(path, "w") as f:
            f.write(text)
        return path

    # expected donors

    def test_expected_donors_same_contents_different_paths(self):
        a = self.write("a.txt", "d2\nd1\n")
        b = self.write("b.txt", "d1\nd2")  # no trailing newline
        c = self.write("c.txt", "d1\r\nd2\r\n\r\n")  # Windows line endings, as in real donor files
        self.assertEqual(qc.expected_donors_from_files([a, b, c]), ["d1", "d2"])

    def test_expected_donors_different_sets(self):
        a = self.write("a.txt", "d1\nd2\n")
        b = self.write("b.txt", "d1\nd3\n")
        self.assertIsNone(qc.expected_donors_from_files([a, b]))

    def test_expected_donors_no_files(self):
        self.assertIsNone(qc.expected_donors_from_files([]))
        self.assertIsNone(qc.expected_donors_from_files(None))

    # conditions

    def test_conditions_file(self):
        path = self.write("cond.txt", "PREFIX\tCONDITION\nlib2\tNeuN+\nlib1\tNeuN-\nunused\tx\n")
        result = qc.resolve_conditions(METADATA, "PREFIX", conditions_file=path)
        self.assertEqual(result.to_dict("list"), {"library": ["lib1", "lib2"], "condition": ["NeuN-", "NeuN+"]})

    def test_conditions_file_missing_library(self):
        path = self.write("cond.txt", "PREFIX\tCONDITION\nlib1\ta\n")
        with self.assertRaisesRegex(ValueError, "lib2"):
            qc.resolve_conditions(METADATA, "PREFIX", conditions_file=path)

    def test_conditions_column(self):
        metadata = METADATA.assign(sort=["hi"] * 5 + ["lo"] * 2)
        result = qc.resolve_conditions(metadata, "PREFIX", condition_column="sort")
        self.assertEqual(result.to_dict("list"), {"library": ["lib1", "lib2"], "condition": ["hi", "lo"]})

    def test_conditions_column_varies_within_library(self):
        metadata = METADATA.assign(sort=["hi", "lo", "hi", "hi", "hi", "lo", "lo"])
        with self.assertRaisesRegex(ValueError, "more than one value.*lib1"):
            qc.resolve_conditions(metadata, "PREFIX", condition_column="sort")

    # donor_qc

    def test_donor_qc_basic(self):
        result = qc.donor_qc_table(METADATA, ["d1", "d2", "d3"])
        self.assertEqual(list(result.columns), ["donor", "num_nuclei", "fraction_nuclei", "total_umis_m",
                                                "median_umis_per_nucleus", "expected"])
        self.assertEqual(list(result["donor"]), ["d1", "d2", "d3", "d9"])
        self.assertEqual(list(result["num_nuclei"]), [4, 1, 0, 1])
        # 6 cells have a donor
        self.assertEqual(list(result["fraction_nuclei"]), [0.667, 0.167, 0.0, 0.167])
        self.assertEqual(result["total_umis_m"].iloc[1], 0.0)  # 200 UMIs
        self.assertTrue(np.isnan(result["total_umis_m"].iloc[2]))
        self.assertEqual(result["median_umis_per_nucleus"].iloc[0], 1500)  # median(1000, 3000, 500, 2000)
        self.assertEqual(result["median_umis_per_nucleus"].iloc[1], 200)
        self.assertTrue(np.isnan(result["median_umis_per_nucleus"].iloc[2]))
        self.assertEqual(list(result["expected"]), [True, True, True, False])

    def test_donor_qc_core_cell_types(self):
        result = qc.donor_qc_table(METADATA, ["d1", "d2", "d3"], core_cell_types=["spn", "absent"],
                                   cell_type_column=CT)
        self.assertEqual(list(result.columns)[-3:], ["median_umis_per_core_cell_type",
                                                     "median_umis_per_nucleus_spn", "expected"])
        self.assertNotIn("median_umis_per_nucleus_absent", result.columns)
        spn = result["median_umis_per_nucleus_spn"]
        self.assertEqual(spn.iloc[0], 2000)  # d1: 1000, 3000
        self.assertTrue(spn.iloc[1:3].isna().all())  # d2 has only glia; d3 has no cells
        self.assertTrue(np.isnan(result["median_umis_per_core_cell_type"].iloc[2]))

    # library summary

    def test_library_cell_summary(self):
        result = qc.library_cell_summary(METADATA, CT)
        self.assertEqual(list(result["library"]), ["lib1", "lib2"])
        self.assertEqual(list(result["num_nuclei"]), [5, 2])
        self.assertEqual(list(result["median_umis_per_nucleus"]), [800, 1025])
        self.assertEqual(list(result["num_labelled_nuclei"]), [5, 1])
        self.assertEqual(list(result["median_umis_per_labelled_nucleus"]), [800, 50])

    # complete_donor_cell_type

    def test_complete_donor_cell_type(self):
        result = qc.complete_donor_cell_type(METADATA, CT, ["d1", "d2", "d3"])
        self.assertEqual(len(result), 4 * 2)  # d1, d2, d3, d9 x glia, spn
        row = result[(result["donor"] == "d3") & (result["cell_type"] == "spn")].iloc[0]
        self.assertEqual((row["num_nuclei"], row["fraction_nuclei"]), (0, 0))
        self.assertTrue(np.isnan(row["median_umis_per_nucleus"]))
        self.assertTrue(row["expected"])
        row = result[(result["donor"] == "d1") & (result["cell_type"] == "spn")].iloc[0]
        self.assertEqual((row["num_nuclei"], row["median_umis_per_nucleus"]), (2, 2000))
        self.assertEqual(set(result.loc[result["donor"] == "d9", "expected"]), {False})
        self.assertEqual(set(result.loc[result["donor"] == "d3", "expected"]), {True})

    # normalize / cell_type_table / donor_library_counts

    def normalized(self, conditions=None):
        return qc.normalize_cell_metadata(METADATA, cell_type_column=CT, conditions=conditions)

    def test_cell_type_table_by_library_has_zero_rows(self):
        metadata = METADATA.assign(**{CT: ["spn"] * 5 + [None, "glia"]})
        result = qc.cell_type_table(qc.normalize_cell_metadata(metadata, cell_type_column=CT), ["library"])
        self.assertEqual(len(result), 4)
        row = result[(result["library"] == "lib1") & (result["cell_type"] == "glia")].iloc[0]
        self.assertEqual((row["num_nuclei"], row["fraction_nuclei"]), (0, 0))
        self.assertTrue(np.isnan(row["median_umis_per_nucleus"]))
        row = result[(result["library"] == "lib2") & (result["cell_type"] == "spn")].iloc[0]
        self.assertEqual((row["num_nuclei"], row["fraction_nuclei"]), (0, 0))
        row = result[(result["library"] == "lib2") & (result["cell_type"] == "glia")].iloc[0]
        self.assertEqual((row["num_nuclei"], row["fraction_nuclei"]), (1, 0.5))


if __name__ == "__main__":
    unittest.main()
