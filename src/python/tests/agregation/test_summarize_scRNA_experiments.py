#!/usr/bin/env python3
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
import contextlib
import csv
import gzip
import io
import os
import shutil
import tempfile
import unittest

import yaml

import dropseq.aggregation.locate_scRNA_artifacts as locator
import dropseq.aggregation.summarize_scRNA_experiments as summarizer
from dropseq.aggregation.locate_scRNA_artifacts import NA
from dropseq.util.storage import LocalStore

UEI = "exp"
ALIGNMENT_REL = "GRCh38_ensembl_v43"
CBRB_REL = ALIGNMENT_REL + "/cbrb/auto"
CELL_SELECTION_REL = CBRB_REL + "/cell_selection/auto"
STD_REL = CELL_SELECTION_REL + "/standard_analysis"
VILLAGE_REL = STD_REL + "/village/donors"

PICARD_HEADER = "## htsjdk.samtools.metrics.StringHeader\n# SomeProgram --INPUT x\n\n"

# Cells A, B are selected; C, D are not.
CELL_FEATURES = """cell_barcode\tnum_transcripts\tnum_reads\tpct_ribosomal\tpct_coding\tpct_intronic\tpct_intergenic\tpct_utr\tpct_genic\tpct_mt
A\t100\t300000\t0\t0.2\t0.6\t0.1\t0.1\t0.9\t0.01
B\t200\t100000\t0\t0.1\t0.7\t0.1\t0.1\t0.9\t0.03
C\t20\t50000\t0\t0.3\t0.5\t0.1\t0.1\t0.9\t0.1
D\t10\t50000\t0\t0.3\t0.5\t0.1\t0.1\t0.9\t0.2
"""

FILES = {
    "properties.yaml": yaml.safe_dump({"library": UEI, "experimentDate": "2024-03-27", "stage": "beginning"}),
    "pipeline_info/execution_trace_2026-10-01_17-08-40.txt":
        "task_id\thash\tname\tstatus\texit\n"
        "1\tab/1\tALIGN (exp)\tFAILED\t1\n"
        "2\tab/2\tALIGN (exp)\tCOMPLETED\t0\n"
        "3\tab/3\tCBRB (exp)\tCACHED\t0\n",
    ALIGNMENT_REL + "/{uei}.ReadQualityMetrics.txt": PICARD_HEADER +
        "aggregate\ttotalReads\tmappedReads\thqMappedReads\thqMappedReadsNoPCRDupes\n"
        "all\t2000000\t1900000\t1800000\t1800000\n",
    ALIGNMENT_REL + "/{uei}.fracIntronicExonic.txt": PICARD_HEADER +
        "PCT_RIBOSOMAL_BASES\tPCT_CODING_BASES\tPCT_UTR_BASES\tPCT_INTRONIC_BASES\tPCT_INTERGENIC_BASES\tSAMPLE\n"
        "0.001\t0.15\t0.09\t0.66\t0.1\t\n",
    ALIGNMENT_REL + "/{uei}.chimeric_read_metrics": PICARD_HEADER +
        "NUM_UMIS\tNUM_MARKED_UMIS\tNUM_T_RICH_UMIS\tNUM_REUSED_UMIS\n100000\t25\t0\t25\n",
    ALIGNMENT_REL + "/{uei}.digital_expression_summary.txt": PICARD_HEADER +
        "CELL_BARCODE\tNUM_GENIC_READS\tNUM_TRANSCRIPTS\tNUM_GENES\n"
        "A\t150\t100\t50\nB\t300\t200\t90\nC\t30\t20\t10\nD\t15\t10\t5\n",
    ALIGNMENT_REL + "/{uei}.numReads_perCell.txt.gz": "#INPUT=x.bam\tTAG=XC\n300\tA\n800\tB\n50\tC\n20\tD\n",
    ALIGNMENT_REL + "/{uei}.cell_features.txt": CELL_FEATURES,
    ALIGNMENT_REL + "/{uei}.alignment_summary.pdf": "",
    CBRB_REL + "/{uei}.cbrb.cell_features.txt":
        "cell_barcode\tnum_transcripts\tnum_retained_transcripts\nA\t100\t90\nB\t200\t150\nC\t20\t0\n",
    CBRB_REL + "/{uei}.cbrb_tearsheet.pdf": "",
    CBRB_REL + "/{uei}_report.html": "",
    CELL_SELECTION_REL + "/{uei}.selectedCellBarcodes.txt": "# CallSTAMPs(outCellFile='x')\nA\nB\n",
    CELL_SELECTION_REL + "/{uei}.cell_selection_assignments_summary.txt":
        "version\tn_STAMPs\tUMIs_ambient\tUMIs_STAMPs_median\tthresholdMessage\tmedianReadsPerUMI\n"
        "v5.0_CBRB\t2\t61.57\t150\t\t1.4989\n",
    CELL_SELECTION_REL + "/{uei}.cell_selection_assignments.pdf": "",
    STD_REL + "/{uei}.selected.digital_expression.txt.gz": "",
    STD_REL + "/{uei}.selected.digital_expression_summary.txt": PICARD_HEADER +
        "CELL_BARCODE\tNUM_GENIC_READS\tNUM_TRANSCRIPTS\tNUM_GENES\nA\t150\t90\t50\nB\t300\t150\t90\n",
    STD_REL + "/{uei}.gmg.digital_expression.txt.gz": "",
    STD_REL + "/{uei}.cmd.tsv": "cell_barcode\tnum_transcripts\tdoublet\tscDblFinder_score\n"
                                "A\t100\tsinglet\t0.1\nB\t200\tdoublet\t0.9\n",
    VILLAGE_REL + "/{uei}.donors.digital_expression.txt.gz": "",
    VILLAGE_REL + "/{uei}.dropulation_summary_stats.txt":
        "expName\ttotal_cells\tpct_all_doublets\tpct_confident_doublets\tsinglets\tassignable_singlets\t"
        "cell_equitability\tdiversity\tequitability\ttotalUMIs\treads_per_umi\n"
        "exp\t2\t3.52\t0.32\t1\t1\t0.8\t2.55\t0.9\t300\t1.55\n",
    VILLAGE_REL + "/{uei}.donors.metacells.txt": "",
    VILLAGE_REL + "/{uei}.dropulation_report.pdf": "",
}


class TestSummarizeScRNAExperiments(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(".tmp", "summarize_scRNA_experiments.")
        self.root = os.path.join(self.tmpDir, UEI + "_output")
        for name, content in FILES.items():
            path = os.path.join(self.root, name.format(uei=UEI))
            os.makedirs(os.path.dirname(path), exist_ok=True)
            opener = gzip.open if path.endswith(".gz") else open
            with opener(path, "wt") as f:
                f.write(content)
        self.store = LocalStore()

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def path(self, name):
        return os.path.join(self.root, name.format(uei=UEI))

    def village_dge(self):
        return self.path(VILLAGE_REL + "/{uei}.donors.digital_expression.txt.gz")

    def selected_dge(self):
        return self.path(STD_REL + "/{uei}.selected.digital_expression.txt.gz")

    def record(self, dge=None):
        return locator.locate_scRNA_artifacts(dge or self.village_dge(), self.store)

    def test_summary(self):
        summary = summarizer.summarize_experiment(self.record(), self.store)
        self.assertEqual(list(summary.keys()), summarizer.SUMMARY_COLUMNS)
        expected = {
            "library": UEI, "exp_date": "2024-03-27", "reference": ALIGNMENT_REL,
            "standard_analysis_dir": "standard_analysis", "run_status": "COMPLETED", "run_date": "2026-10-01",
            "total_reads_M": 2.0, "pct_mapped_reads": 95.0, "pct_hq_mapped_reads": 90.0,
            "pct_genic": 90.0, "pct_exonic": 24.0, "pct_intronic": 66.0, "pct_intergenic": 10.0,
            "pct_ribosomal": 0.1,
            "n_STAMPs": 2, "median_cell_UMIs": 150, "median_ambient_UMIs": 15, "pct_ambient": 10.0,
            # A: 300/100, B: 800/200
            "reads_per_umi": 3.5,
            "pct_chimeric_UMIs": 0.025, "median_cell_UMIs_cbrb": 120,
            # (90 + 150) / (100 + 200)
            "cbrb_retained_umi_pct": 80.0,
            "pct_all_doublets": 3.52, "pct_confident_doublets": 0.32, "singlets": 1, "assignable_singlets": 1,
            "reads_per_umi_donor": 1.55, "totalUMIs": 300, "donor_equitability": 0.9, "cell_equitability": 0.8,
            "scDblFinder_pct_doublets": 50.0,
            "AlignmentPDF": self.path(ALIGNMENT_REL + "/{uei}.alignment_summary.pdf"),
            "CBRB_tear_sheet_path": self.path(CBRB_REL + "/{uei}.cbrb_tearsheet.pdf"),
            "CBRB_html_path": self.path(CBRB_REL + "/{uei}_report.html"),
            "SelectionPlot": self.path(CELL_SELECTION_REL + "/{uei}.cell_selection_assignments.pdf"),
            "DropulationReport": self.path(VILLAGE_REL + "/{uei}.dropulation_report.pdf"),
            "Metacell": self.path(VILLAGE_REL + "/{uei}.donors.metacells.txt"),
            "gmg_dge_path": self.path(STD_REL + "/{uei}.gmg.digital_expression.txt.gz"),
        }
        for column, value in expected.items():
            self.assertEqual(summary[column], value, column)
        # read-weighted: (0.01*300000 + 0.03*100000 + 0.1*50000 + 0.2*50000) / 500000
        approximately = {
            "pct_mt": 4.2, "nonselected_cbcs_pct_mt": 15.0, "nonselected_cbcs_total_reads_M": 0.1,
            "selected_cells_total_reads_M": 0.4,
            "selected_cells_pct_genic": 90.0, "selected_cells_pct_exonic": 25.0, "selected_cells_pct_intronic": 65.0,
            "selected_cells_pct_intergenic": 10.0, "selected_cells_pct_ribosomal": 0.0, "selected_cells_pct_mt": 2.0,
        }
        for column, value in approximately.items():
            self.assertAlmostEqual(summary[column], value, msg=column)

    def test_selected_dge_input(self):
        self.assertEqual(summarizer.summarize_experiment(self.record(self.selected_dge()), self.store),
                         summarizer.summarize_experiment(self.record(), self.store))

    def test_no_village(self):
        shutil.rmtree(os.path.dirname(os.path.dirname(self.village_dge())))
        summary = summarizer.summarize_experiment(self.record(self.selected_dge()), self.store)
        for column in list(summarizer.DROPULATION_COLUMNS.values()) + ["DropulationReport", "Metacell"]:
            self.assertIsNone(summary[column], column)
        self.assertEqual(summary["scDblFinder_pct_doublets"], 50.0)

    def test_missing_inputs_are_na(self):
        os.remove(self.path(ALIGNMENT_REL + "/{uei}.cell_features.txt"))
        os.remove(self.path(ALIGNMENT_REL + "/{uei}.numReads_perCell.txt.gz"))
        shutil.rmtree(self.path("pipeline_info"))
        summary = summarizer.summarize_experiment(self.record(), self.store)
        for column in ["run_status", "run_date", "pct_mt", "selected_cells_pct_mt", "nonselected_cbcs_pct_mt",
                       "reads_per_umi"]:
            self.assertIsNone(summary[column], column)
        # the other selection stats don't need the missing files
        self.assertEqual(summary["median_cell_UMIs"], 150)

    def test_missing_selected_barcodes(self):
        os.remove(self.path(CELL_SELECTION_REL + "/{uei}.selectedCellBarcodes.txt"))
        summary = summarizer.summarize_experiment(self.record(), self.store)
        for column in ["n_STAMPs", "median_cell_UMIs", "median_cell_UMIs_cbrb", "cbrb_retained_umi_pct",
                       "selected_cells_pct_mt"]:
            self.assertIsNone(summary[column], column)
        self.assertEqual(summary["total_reads_M"], 2.0)

    def write_trace(self, rows):
        with open(self.path("pipeline_info/execution_trace_2026-10-01_17-08-40.txt"), "w") as f:
            f.write("task_id\tname\tstatus\n")
            for i, (name, status) in enumerate(rows):
                f.write(f"{i}\t{name}\t{status}\n")

    def test_run_status_failed(self):
        self.write_trace([("ALIGN", "COMPLETED"), ("CBRB", "FAILED"), ("CBRB", "ABORTED")])
        self.assertEqual(summarizer.run_status(self.record(), self.store),
                         {"run_status": "FAILED", "run_date": "2026-10-01"})

    def test_run_status_retried(self):
        self.write_trace([("CBRB", "FAILED"), ("CBRB", "COMPLETED"), ("ALIGN", "CACHED")])
        self.assertEqual(summarizer.run_status(self.record(), self.store)["run_status"], "COMPLETED")

    def test_read_table_stops_at_histogram(self):
        with open(self.path(ALIGNMENT_REL + "/{uei}.chimeric_read_metrics"), "a") as f:
            f.write("\n## HISTOGRAM\nbin\tcount\n1\t2\n")
        table = summarizer.read_table(self.path(ALIGNMENT_REL + "/{uei}.chimeric_read_metrics"), self.store)
        self.assertEqual(list(table.columns), ["NUM_UMIS", "NUM_MARKED_UMIS", "NUM_T_RICH_UMIS", "NUM_REUSED_UMIS"])
        self.assertEqual(len(table), 1)

    def test_tearsheet(self):
        rows = summarizer.standard_analysis_tearsheet(self.record(), self.store)
        self.assertEqual(rows, [
            ("PF reads", "2000000"),
            ("HQ mapped reads", "1800000"),
            ("% HQ mapped reads", "90.0%"),
            ("% chimeric UMIs", "0.03%"),
            ("median Reads/UMI(STAMPS)", "1.5"),
            ("STAMPS", "2"),
            # median of the selected cells after CBRB (90, 150), not UMIs_STAMPs_median from the summary file
            ("median UMIs(STAMPS)", "120"),
            ("ambient peak(UMIs)", "61"),
        ])

    def test_tearsheet_median_umis_missing_input(self):
        os.remove(self.path(STD_REL + "/{uei}.selected.digital_expression_summary.txt"))
        rows = dict(summarizer.standard_analysis_tearsheet(self.record(), self.store))
        self.assertEqual(rows["median UMIs(STAMPS)"], "NA")

    def test_tearsheet_optional_rows(self):
        with open(self.path(CELL_SELECTION_REL + "/{uei}.cell_selection_assignments_summary.txt"), "w") as f:
            f.write("version\tn_STAMPs\tUMIs_ambient\tUMIs_STAMPs_median\tthresholdMessage\tmedianReadsPerUMI\n"
                    "v5\t2\t\t150\tmanual threshold\t2\n")
        os.remove(self.path(ALIGNMENT_REL + "/{uei}.chimeric_read_metrics"))
        rows = dict(summarizer.standard_analysis_tearsheet(self.record(), self.store))
        self.assertEqual(rows["M/O"], "manual threshold")
        self.assertEqual(rows["% chimeric UMIs"], "NA")
        self.assertEqual(rows["median Reads/UMI(STAMPS)"], "2")
        self.assertNotIn("ambient peak(UMIs)", rows)

    def read_summary(self, path):
        with open(path) as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def test_write_summary(self):
        no_village = self.record()
        no_village["dropulation_summary_stats"] = NA
        no_village["uei"] = "other"
        summary = summarizer.summarize_datasets([self.record(), no_village], self.store)
        out = io.StringIO()
        summarizer.write_summary(summary, out)
        rows = list(csv.DictReader(io.StringIO(out.getvalue()), delimiter="\t"))
        self.assertEqual(list(rows[0].keys()), summarizer.SUMMARY_COLUMNS)
        self.assertEqual([r["library"] for r in rows], [UEI, "other"])
        self.assertEqual(rows[0]["n_STAMPs"], "2")
        self.assertEqual(rows[0]["total_reads_M"], "2")
        self.assertEqual(rows[0]["pct_ribosomal"], "0.1")
        self.assertEqual(rows[0]["singlets"], "1")
        self.assertEqual(rows[1]["singlets"], "NA")

        out = io.StringIO()
        no_village["dropulation_report"] = NA
        no_village["meta_cell"] = NA
        summarizer.write_summary(summarizer.summarize_datasets([no_village], self.store), out,
                                 remove_empty_columns=True)
        header = out.getvalue().splitlines()[0].split("\t")
        self.assertNotIn("singlets", header)
        self.assertNotIn("Metacell", header)
        self.assertIn("n_STAMPs", header)

    def test_main_dge(self):
        summary = os.path.join(self.tmpDir, "summary.tsv")
        tearsheets = os.path.join(self.tmpDir, "tearsheets")
        self.assertEqual(summarizer.main(["--dge", self.village_dge(), "--summary", summary,
                                          "--tearsheet-dir", tearsheets]), 0)
        rows = self.read_summary(summary)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["scDblFinder_pct_doublets"], "50")
        self.assertEqual(os.listdir(tearsheets), [f"{UEI}.standard_analysis_tearsheet.txt"])
        with open(os.path.join(tearsheets, f"{UEI}.standard_analysis_tearsheet.txt")) as f:
            lines = f.read().splitlines()
        self.assertEqual(lines[0], "label\tvalue")
        self.assertEqual(lines[1], "PF reads\t2000000")

    def test_main_tearsheet_only(self):
        tearsheets = os.path.join(self.tmpDir, "tearsheets")
        self.assertEqual(summarizer.main(["--dge", self.village_dge(), "--tearsheet-dir", tearsheets]), 0)
        self.assertEqual(os.listdir(tearsheets), [f"{UEI}.standard_analysis_tearsheet.txt"])
        self.assertEqual(sorted(os.listdir(self.tmpDir)), sorted([UEI + "_output", "tearsheets"]))

    def test_main_summary_only(self):
        summary = os.path.join(self.tmpDir, "summary.tsv")
        self.assertEqual(summarizer.main(["--dge", self.village_dge(), "--summary", summary]), 0)
        self.assertEqual(len(self.read_summary(summary)), 1)
        self.assertEqual(sorted(os.listdir(self.tmpDir)), sorted([UEI + "_output", "summary.tsv"]))

    def test_main_requires_an_output(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                summarizer.main(["--dge", self.village_dge()])

    def test_main_artifacts(self):
        artifacts = os.path.join(self.tmpDir, "artifacts.yaml")
        self.assertEqual(locator.main(["--dge", self.village_dge(), "--output", artifacts]), 0)
        from_artifacts = os.path.join(self.tmpDir, "from_artifacts.tsv")
        from_dge = os.path.join(self.tmpDir, "from_dge.tsv")
        self.assertEqual(summarizer.main(["--artifacts", artifacts, "--summary", from_artifacts]), 0)
        self.assertEqual(summarizer.main(["--dge", self.village_dge(), "--summary", from_dge]), 0)
        self.assertEqual(self.read_summary(from_artifacts), self.read_summary(from_dge))

    def test_duplicate_uei_tearsheets(self):
        record = self.record()
        output_dir = os.path.join(self.tmpDir, "tearsheets")
        with self.assertRaises(ValueError):
            summarizer.write_tearsheets([record, record], output_dir, self.store)
        self.assertFalse(os.path.exists(output_dir))


if __name__ == '__main__':
    unittest.main()
