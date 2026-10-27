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
import re
import shutil
import tempfile
import unittest

import pandas as pd

from dropseq.aggregation.plots import aggregation_qc

CT = "model_neighborhood_name"
DONORS = ["d1", "d2", "d3"]

PDFS = ["library_metrics", "filtered_library_metrics", "donor_summary", "donor_nuclei_counts", "cell_type_counts",
        "cell_type_median_UMIs", "cell_type_proportions", "donor_cell_type_distributions"]


def page_count(path):
    with open(path, "rb") as f:
        return len(re.findall(rb"/Type /Page\b", f.read()))


def make_cell_metadata():
    rows = []
    for library, condition in (("lib1", "hi"), ("lib2", "hi"), ("lib3", "lo")):
        for i in range(30):
            rows.append({"cell_barcode": f"{library}_{i}", "PREFIX": library, "sort": condition,
                         "donor": ["d1", "d2", "d9"][i % 3],   # d3 is expected but never seen; d9 is not expected
                         CT: ["spn", "glia", "spn", None][i % 4],
                         "num_retained_transcripts": 500 + 37 * i + 100 * len(library)})
    return pd.DataFrame(rows)


def make_library_metrics():
    return pd.DataFrame({"library": ["lib1", "lib2", "lib3"], "total_reads_M": [100.0, 120.0, 90.0],
                         "pct_hq_mapped_reads": [90.0, 91.0, 85.0], "pct_chimeric_UMIs": [3.0, 2.0, 4.0],
                         "n_STAMPs": [1000, 1100, 900], "median_cell_UMIs_cbrb": [2000, 2100, 1900],
                         "pct_confident_doublets": [5.0, 6.0, 4.0], "assignable_singlets": [800, 900, 700],
                         "reads_per_umi_donor": [3.1, 3.0, 2.9], "totalUMIs": [2.0e6, 2.5e6, 1.5e6],
                         "donor_equitability": [0.9, 0.8, 0.7], "scDblFinder_pct_doublets": [None, None, None]})


class TestAggregationQc(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(".tmp", "aggregation_qc.")
        self.out = os.path.join(self.tmpDir, "out")
        self.cell_metadata = self.write_table("cells.txt", make_cell_metadata())
        self.library_metrics = self.write_table("library_metrics.txt", make_library_metrics())
        # the donor file of each library, with no trailing newline
        self.donor_files = []
        for library in ("lib1", "lib2", "lib3"):
            self.donor_files.append(os.path.join(self.tmpDir, f"{library}.donors.txt"))
            with open(self.donor_files[-1], "w") as f:
                f.write("\n".join(DONORS))

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def write_table(self, name, table):
        path = os.path.join(self.tmpDir, name)
        table.to_csv(path, sep="\t", index=False)
        return path

    def run_qc(self, *extra, donor_files=None, cell_metadata=None, expected_status=0):
        donor_files = self.donor_files if donor_files is None else donor_files
        args = ["--cell-metadata", cell_metadata or self.cell_metadata, "--library-metrics", self.library_metrics,
                "--analysis-id", "AGG", "--output-dir", self.out, *(["--donor-file", *donor_files] if donor_files
                                                                    else []), *extra]
        self.assertEqual(aggregation_qc.main(args), expected_status)

    def produced(self):
        return sorted(os.listdir(self.out)) if os.path.isdir(self.out) else []

    def expected_files(self, pdfs, donor_qc):
        return sorted([f"AGG.{p}.pdf" for p in pdfs] + (["AGG.donor_qc.txt"] if donor_qc else []))

    def test_everything(self):
        self.run_qc("--cell-type-column", CT, "--conditions-file", self.write_table(
            "conditions.txt", pd.DataFrame({"PREFIX": ["lib1", "lib2", "lib3"], "CONDITION": ["hi", "hi", "lo"]})),
            "--core-cell-types", self.write_lines("core.txt", ["spn"]))
        self.assertEqual(self.produced(), self.expected_files(PDFS, True))
        for pdf in PDFS:
            self.assertGreater(page_count(os.path.join(self.out, f"AGG.{pdf}.pdf")), 0, pdf)

    def test_no_cell_type_no_conditions(self):
        self.run_qc()
        self.assertEqual(self.produced(), self.expected_files(
            ["library_metrics", "filtered_library_metrics", "donor_summary", "donor_nuclei_counts"], True))
        # 9 metrics are available; scDblFinder is all NA
        self.assertEqual(page_count(os.path.join(self.out, "AGG.library_metrics.pdf")), 9)
        self.assertEqual(page_count(os.path.join(self.out, "AGG.filtered_library_metrics.pdf")), 2)
        # 4 donor bar pages
        self.assertEqual(page_count(os.path.join(self.out, "AGG.donor_summary.pdf")), 4)
        # total, by library, donor composition (no condition pages)
        self.assertEqual(page_count(os.path.join(self.out, "AGG.donor_nuclei_counts.pdf")), 4)

    def test_cell_type_adds_pages(self):
        self.run_qc("--cell-type-column", CT)
        self.assertEqual(self.produced(), self.expected_files(PDFS, True))
        self.assertEqual(page_count(os.path.join(self.out, "AGG.library_metrics.pdf")), 11)
        # stacked, then box + facets for counts, proportions and median UMIs
        self.assertEqual(page_count(os.path.join(self.out, "AGG.donor_summary.pdf")), 4 + 1 + 3 * 2)
        # overall, per library, per donor box and facets
        self.assertEqual(page_count(os.path.join(self.out, "AGG.cell_type_counts.pdf")), 4)
        self.assertEqual(page_count(os.path.join(self.out, "AGG.donor_cell_type_distributions.pdf")), 2)

    def test_conditions_add_pages(self):
        self.run_qc("--cell-type-column", CT, "--condition-column", "sort")
        # + per condition, dodged box and free facets
        self.assertEqual(page_count(os.path.join(self.out, "AGG.cell_type_counts.pdf")), 7)
        self.assertEqual(page_count(os.path.join(self.out, "AGG.donor_cell_type_distributions.pdf")), 3)
        # + by condition stacked and dodged
        self.assertEqual(page_count(os.path.join(self.out, "AGG.donor_nuclei_counts.pdf")), 6)

    def test_donor_qc_table(self):
        self.run_qc("--cell-type-column", CT, "--core-cell-types", self.write_lines("core.txt", ["spn", "absent"]))
        table = pd.read_csv(os.path.join(self.out, "AGG.donor_qc.txt"), sep="\t", dtype={"donor": str})
        self.assertEqual(list(table.columns), ["donor", "num_nuclei", "fraction_nuclei", "total_umis_m",
                                               "median_umis_per_nucleus", "median_umis_per_core_cell_type",
                                               "median_umis_per_nucleus_spn", "expected"])
        self.assertEqual(list(table["donor"]), ["d1", "d2", "d3", "d9"])
        self.assertEqual(list(table["num_nuclei"]), [30, 30, 0, 30])
        self.assertEqual(list(table["expected"]), [True, True, True, False])
        self.assertTrue(table.loc[2, "median_umis_per_nucleus"] != table.loc[2, "median_umis_per_nucleus"])

    def write_lines(self, name, lines):
        path = os.path.join(self.tmpDir, name)
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
        return path

    def test_different_donor_files_skip_donor_outputs(self):
        other = self.write_lines("other.txt", ["d1", "d2"])
        self.run_qc("--cell-type-column", CT,
                    donor_files=[other, other, self.write_lines("donors3.txt", ["d1", "d2", "d3"])])
        self.assertEqual(self.produced(), self.expected_files(
            ["library_metrics", "filtered_library_metrics", "cell_type_counts", "cell_type_median_UMIs",
             "cell_type_proportions"], False))

    def test_no_donor_file_skips_donor_outputs(self):
        self.run_qc(donor_files=[])
        self.assertEqual(self.produced(), self.expected_files(["library_metrics", "filtered_library_metrics"], False))

    def test_no_donor_column_skips_donor_outputs(self):
        cell_metadata = self.write_table("nodonor.txt", make_cell_metadata().drop(columns="donor"))
        self.run_qc(cell_metadata=cell_metadata)
        self.assertEqual(self.produced(), self.expected_files(["library_metrics", "filtered_library_metrics"], False))

    def test_condition_varies_within_library_is_invalid(self):
        cells = make_cell_metadata()
        cells.loc[0, "sort"] = "lo"
        self.run_qc("--condition-column", "sort", cell_metadata=self.write_table("varying.txt", cells),
                    expected_status=1)

    def test_numeric_looking_donors(self):
        cells = make_cell_metadata().assign(donor=lambda d: d["donor"].map({"d1": "007", "d2": "8", "d9": "9"}))
        donor_file = self.write_lines("numeric.txt", ["007", "8", "9"])
        self.run_qc(cell_metadata=self.write_table("numeric_cells.txt", cells),
                    donor_files=[donor_file])
        table = pd.read_csv(os.path.join(self.out, "AGG.donor_qc.txt"), sep="\t", dtype={"donor": str})
        self.assertEqual(list(table["donor"]), ["007", "8", "9"])
        self.assertEqual(set(table["expected"]), {True})


if __name__ == "__main__":
    unittest.main()
