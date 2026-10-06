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
import io
import os
import shutil
import tempfile
import unittest

import yaml

import dropseq.aggregation.locate_scRNA_artifacts as locator
from dropseq.aggregation.locate_scRNA_artifacts import NA
from dropseq.util.storage import LocalStore

ALIGNMENT_REL = "GRCh38_ensembl_v43"
CBRB_REL = ALIGNMENT_REL + "/cbrb/auto"
CELL_SELECTION_REL = CBRB_REL + "/cell_selection/auto"
STD_REL = CELL_SELECTION_REL + "/standard_analysis"
VILLAGE_NAME = "donors.donor_assignment_filtered.2024-03-27_SMRI_V11"

# Mirrors the layout of the Nextflow output in
# gs://mccarroll_scrnaseq_central_temp/nextflow/pipeline-testing/jn_SMRI_V11_rxn1_output/
ROOT_FILES = ["{uei}.corrected_barcode_metrics", "{uei}.0.unmapped.bam"]
PIPELINE_INFO_FILES = ["execution_trace_2026-09-30_10-00-00.txt", "execution_trace_2026-10-01_17-08-40.txt",
                       "execution_report_2026-10-01_17-08-40.html", "params_2026-10-01_17-09-03.json"]
ALIGNMENT_FILES = ["properties.yaml", "{uei}.fracIntronicExonic.txt", "{uei}.ReadQualityMetrics.txt",
                   "{uei}.numReads_perCell.txt.gz", "{uei}.cell_features.txt", "{uei}.digital_expression.txt.gz",
                   "{uei}.digital_expression_summary.txt", "{uei}.chimeric_read_metrics",
                   "{uei}.alignment_summary.pdf"]
CBRB_FILES = ["properties.yaml", "{uei}.cbrb.cell_features.txt", "{uei}_report.html", "{uei}.cbrb_tearsheet.pdf",
              "{uei}.svm_cbrb_parameter_estimation.pdf", "{uei}.cbrb.digital_expression.txt.gz"]
CELL_SELECTION_FILES = ["properties.yaml", "{uei}.selectedCellBarcodes.txt",
                        "{uei}.cell_selection_assignments_summary.txt", "{uei}.cell_selection_assignments.pdf"]
STD_FILES = ["properties.yaml", "{uei}.cmd.tsv", "{uei}.selected.digital_expression.txt.gz",
             "{uei}.selected.digital_expression_summary.txt", "{uei}.gmg.digital_expression.txt.gz",
             "{uei}.gmg.digital_expression_summary.txt", "{uei}.metagene.digital_expression.txt.gz",
             "{uei}.metagene.digital_expression_summary.txt"]
MMC_FILES = ["properties.yaml", "{uei}.csv", "{uei}.json", "{uei}.cell_type_counts.tsv"]
VILLAGE_FILES = ["properties.yaml", "{uei}.cmd.tsv", "{uei}.dropulation_summary_stats.txt",
                 "{uei}.donors.digital_expression.txt.gz", "{uei}.donors.digital_expression_summary.txt",
                 "{uei}.donor_cell_map.txt", "{uei}.donor_assignments.txt",
                 "{uei}.donors.metacells.txt", "{uei}.dropulation_report.pdf", "{uei}.dropulation_tearsheet.pdf"]


class TestLocateScRNAArtifacts(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(".tmp", "locate_scRNA_artifacts.")

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def touch(self, directory, names, uei):
        os.makedirs(directory, exist_ok=True)
        for name in names:
            open(os.path.join(directory, name.format(uei=uei)), "w").close()

    def reference_dir(self):
        return os.path.join(self.tmpDir, "reference", "GRCh38")

    def make_experiment(self, uei, mmc_references=("HMBA_Human_WB_v0.5",), villages=(VILLAGE_NAME,),
                        library=None, reference_files=("GRCh38.fasta.gz", "GRCh38.reduced.gtf")):
        """
        Build a Nextflow output tree for one experiment and return the selected DGE path.  The reference of the
        alignment is in a directory shared by the experiments that has reference_files, or is not recorded if
        reference_files is None.
        """
        root = os.path.join(self.tmpDir, uei + "_output")
        os.makedirs(root)
        with open(os.path.join(root, "properties.yaml"), "w") as f:
            yaml.safe_dump({"library": library or uei, "stage": "beginning"}, f)
        self.touch(root, ROOT_FILES, uei)
        self.touch(os.path.join(root, "pipeline_info"), PIPELINE_INFO_FILES, uei)
        alignment = os.path.join(root, ALIGNMENT_REL)
        self.touch(alignment, ALIGNMENT_FILES, uei)
        alignment_properties = {"stage": "alignment"}
        if reference_files is not None:
            self.touch(self.reference_dir(), reference_files, uei)
            alignment_properties["reference"] = os.path.join(self.reference_dir(), "GRCh38.fasta.gz")
        with open(os.path.join(alignment, "properties.yaml"), "w") as f:
            yaml.safe_dump(alignment_properties, f)
        self.touch(os.path.join(root, CBRB_REL), CBRB_FILES, uei)
        self.touch(os.path.join(root, CELL_SELECTION_REL), CELL_SELECTION_FILES, uei)
        std = os.path.join(root, STD_REL)
        self.touch(std, STD_FILES, uei)
        for reference in mmc_references:
            self.touch(os.path.join(std, "mmc", reference), MMC_FILES, uei)
        for village in villages:
            self.touch(os.path.join(std, "village", village), VILLAGE_FILES, uei)
        return os.path.join(std, f"{uei}.selected.digital_expression.txt.gz")

    def test_derive_stage_dirs(self):
        root = "gs://bucket/nextflow/exp_output"
        dge = f"{root}/{STD_REL}/exp.selected.digital_expression.txt.gz"
        self.assertEqual(locator.derive_stage_dirs(dge), {
            "root": root,
            "alignment": f"{root}/{ALIGNMENT_REL}",
            "cbrb": f"{root}/{CBRB_REL}",
            "cell_selection": f"{root}/{CELL_SELECTION_REL}",
            "std_analysis": f"{root}/{STD_REL}",
        })

    def test_derive_stage_dirs_village(self):
        root = "gs://bucket/nextflow/exp_output"
        village = f"{root}/{STD_REL}/village/{VILLAGE_NAME}"
        stage_dirs = locator.derive_stage_dirs(f"{village}/exp.donors.digital_expression.txt.gz")
        self.assertEqual(stage_dirs["std_analysis"], f"{root}/{STD_REL}")
        self.assertEqual(stage_dirs["village"], village)
        self.assertEqual(stage_dirs["root"], root)

    def test_village_dge(self):
        uei = "SMRI_V11_rxn1"
        std = os.path.dirname(self.make_experiment(uei))
        village = os.path.join(std, "village", VILLAGE_NAME)
        dge = os.path.join(village, f"{uei}.donors.digital_expression.txt.gz")
        record = locator.locate_scRNA_artifacts(dge)
        self.assertEqual(list(record.keys()), locator.CANONICAL_FIELDS)
        self.assertEqual(record["uei"], uei)
        self.assertEqual(record["standard_analysis_dir"], std)
        self.assertEqual(record["dge_donors"], dge)
        self.assertEqual(record["dge_summary_donors"],
                         os.path.join(village, f"{uei}.donors.digital_expression_summary.txt"))
        self.assertEqual(record["dge_selected_cells"], os.path.join(std, f"{uei}.selected.digital_expression.txt.gz"))
        self.assertEqual(record["dge_summary_selected_cells"],
                         os.path.join(std, f"{uei}.selected.digital_expression_summary.txt"))
        self.assertEqual(record["cell_metadata"], os.path.join(village, f"{uei}.cmd.tsv"))
        self.assertEqual(record["village_properties"], os.path.join(village, "properties.yaml"))
        self.assertEqual(record["mmc"]["HMBA_Human_WB_v0.5"]["mmc_annotations"],
                         os.path.join(std, "mmc", "HMBA_Human_WB_v0.5", f"{uei}.csv"))

    def test_derive_stage_dirs_bad_layout(self):
        with self.assertRaises(ValueError):
            locator.derive_stage_dirs(f"gs://bucket/exp_output/{ALIGNMENT_REL}/exp.digital_expression.txt.gz")
        with self.assertRaises(ValueError):
            locator.derive_stage_dirs(
                f"gs://bucket/exp_output/{ALIGNMENT_REL}/cbrb/auto/foo/auto/standard_analysis/"
                f"exp.selected.digital_expression.txt.gz")

    def test_single_dge(self):
        uei = "SMRI_V11_rxn1"
        dge = self.make_experiment(uei)
        record = locator.locate_scRNA_artifacts(dge)
        self.assertEqual(list(record.keys()), locator.CANONICAL_FIELDS)
        root = os.path.join(self.tmpDir, uei + "_output")
        std = os.path.join(root, STD_REL)
        village = os.path.join(std, "village", VILLAGE_NAME)
        self.assertEqual(record["uei"], uei)
        self.assertEqual(record["root_properties"], os.path.join(root, "properties.yaml"))
        self.assertEqual(record["alignment_dir"], os.path.join(root, ALIGNMENT_REL))
        self.assertEqual(record["dge_unfiltered"],
                         os.path.join(root, ALIGNMENT_REL, f"{uei}.digital_expression.txt.gz"))
        self.assertEqual(record["dge_summary_unfiltered"],
                         os.path.join(root, ALIGNMENT_REL, f"{uei}.digital_expression_summary.txt"))
        self.assertEqual(record["cbrb_properties"], os.path.join(root, CBRB_REL, "properties.yaml"))
        self.assertEqual(record["selected_cell_barcodes"],
                         os.path.join(root, CELL_SELECTION_REL, f"{uei}.selectedCellBarcodes.txt"))
        self.assertEqual(record["standard_analysis_dir"], std)
        self.assertEqual(record["dge_selected_cells"], dge)
        self.assertEqual(record["dge_summary_selected_cells"],
                         os.path.join(std, f"{uei}.selected.digital_expression_summary.txt"))
        self.assertEqual(record["cell_metadata"], os.path.join(village, f"{uei}.cmd.tsv"))
        self.assertEqual(record["selected_cell_metadata"], os.path.join(std, f"{uei}.cmd.tsv"))
        self.assertEqual(record["dge_gmg"], os.path.join(std, f"{uei}.gmg.digital_expression.txt.gz"))
        self.assertEqual(record["dge_summary_gmg"], os.path.join(std, f"{uei}.gmg.digital_expression_summary.txt"))
        self.assertEqual(record["cell_features"], os.path.join(root, ALIGNMENT_REL, f"{uei}.cell_features.txt"))
        self.assertEqual(record["cbrb_cell_features"], os.path.join(root, CBRB_REL, f"{uei}.cbrb.cell_features.txt"))
        # the most recent of the two traces
        self.assertEqual(record["pipeline_trace"],
                         os.path.join(root, "pipeline_info", "execution_trace_2026-10-01_17-08-40.txt"))
        self.assertEqual(record["dge_donors"], os.path.join(village, f"{uei}.donors.digital_expression.txt.gz"))
        self.assertEqual(record["dge_summary_donors"],
                         os.path.join(village, f"{uei}.donors.digital_expression_summary.txt"))
        self.assertEqual(record["dropulation_tearsheet"], os.path.join(village, f"{uei}.dropulation_tearsheet.pdf"))
        self.assertEqual(record["mmc"]["HMBA_Human_WB_v0.5"]["mmc_annotations"],
                         os.path.join(std, "mmc", "HMBA_Human_WB_v0.5", f"{uei}.csv"))
        # donor is only set by the user, so it is NA when none is given
        self.assertEqual(record["donor"], NA)
        for field in locator.CANONICAL_FIELDS:
            if field != "donor":
                self.assertNotEqual(record[field], NA, field)

    def test_other_dge_type(self):
        uei = "exp"
        selected = self.make_experiment(uei)
        dge = selected.replace(".selected.", ".metagene.")
        record = locator.locate_scRNA_artifacts(dge)
        fields = list(record.keys())
        extra = ("dge_metagene", "dge_summary_metagene")
        # the extra pair immediately follows the gmg pair
        i = fields.index("dge_summary_gmg")
        self.assertEqual(fields[i + 1:i + 3], list(extra))
        self.assertEqual([f for f in fields if f not in extra], locator.CANONICAL_FIELDS)
        self.assertEqual(record["dge_metagene"], dge)
        self.assertEqual(record["dge_summary_metagene"],
                         os.path.join(os.path.dirname(dge), f"{uei}.metagene.digital_expression_summary.txt"))
        self.assertEqual(record["dge_selected_cells"], selected)

    def test_gmg_dge(self):
        selected = self.make_experiment("exp")
        dge = selected.replace(".selected.", ".gmg.")
        record = locator.locate_scRNA_artifacts(dge)
        self.assertEqual(list(record.keys()), locator.CANONICAL_FIELDS)
        self.assertEqual(record["dge_gmg"], dge)
        self.assertEqual(record["dge_selected_cells"], selected)

    def test_no_pipeline_trace(self):
        uei = "exp"
        dge = self.make_experiment(uei)
        pipeline_info = os.path.join(self.tmpDir, uei + "_output", "pipeline_info")
        for name in os.listdir(pipeline_info):
            if name.startswith("execution_trace_"):
                os.remove(os.path.join(pipeline_info, name))
        self.assertEqual(locator.locate_scRNA_artifacts(dge)["pipeline_trace"], NA)
        shutil.rmtree(pipeline_info)
        self.assertEqual(locator.locate_scRNA_artifacts(dge)["pipeline_trace"], NA)

    def test_open_text(self):
        import gzip
        path = os.path.join(self.tmpDir, "a.txt.gz")
        with gzip.open(path, "wt") as f:
            f.write("x\ty\n")
        with LocalStore().open_text(path) as f:
            self.assertEqual(f.read(), "x\ty\n")

    def test_dge_not_named_for_uei_raises(self):
        dge = self.make_experiment("exp")
        other = os.path.join(os.path.dirname(dge), "other.selected.digital_expression.txt.gz")
        open(other, "w").close()
        with self.assertRaises(ValueError):
            locator.locate_scRNA_artifacts(other)

    def test_missing_scalar_is_na(self):
        uei = "exp"
        dge = self.make_experiment(uei)
        os.remove(os.path.join(self.tmpDir, uei + "_output", ALIGNMENT_REL, f"{uei}.chimeric_read_metrics"))
        os.remove(os.path.join(os.path.dirname(dge), f"{uei}.selected.digital_expression_summary.txt"))
        record = locator.locate_scRNA_artifacts(dge)
        self.assertEqual(record["chimeric_metrics"], NA)
        self.assertEqual(record["dge_summary_selected_cells"], NA)

    def test_missing_dge_raises(self):
        dge = self.make_experiment("exp")
        os.remove(dge)
        with self.assertRaises(FileNotFoundError):
            locator.locate_scRNA_artifacts(dge)

    def test_missing_library_raises(self):
        dge = self.make_experiment("exp")
        with open(os.path.join(self.tmpDir, "exp_output", "properties.yaml"), "w") as f:
            yaml.safe_dump({"stage": "beginning"}, f)
        with self.assertRaises(ValueError):
            locator.locate_scRNA_artifacts(dge)

    def test_no_village(self):
        uei = "exp"
        dge = self.make_experiment(uei, villages=())
        record = locator.locate_scRNA_artifacts(dge)
        for field, stage, _ in locator.STAGE_FILES:
            if stage == "village":
                self.assertEqual(record[field], NA, field)
        self.assertEqual(record["cell_metadata"], os.path.join(os.path.dirname(dge), f"{uei}.cmd.tsv"))
        self.assertEqual(record["selected_cell_metadata"], os.path.join(os.path.dirname(dge), f"{uei}.cmd.tsv"))

    def test_village_without_cell_metadata_falls_back(self):
        uei = "exp"
        dge = self.make_experiment(uei)
        os.remove(os.path.join(os.path.dirname(dge), "village", VILLAGE_NAME, f"{uei}.cmd.tsv"))
        record = locator.locate_scRNA_artifacts(dge)
        self.assertEqual(record["cell_metadata"], os.path.join(os.path.dirname(dge), f"{uei}.cmd.tsv"))

    def test_multiple_villages_raises(self):
        dge = self.make_experiment("exp", villages=("v1", "v2"))
        with self.assertRaises(ValueError):
            locator.locate_scRNA_artifacts(dge)

    def test_multiple_mmc(self):
        uei = "exp"
        dge = self.make_experiment(uei, mmc_references=("refB", "refA"))
        os.remove(os.path.join(os.path.dirname(dge), "mmc", "refB", f"{uei}.csv"))
        record = locator.locate_scRNA_artifacts(dge)
        mmc = os.path.join(os.path.dirname(dge), "mmc")
        self.assertEqual(list(record["mmc"].keys()), ["refA", "refB"])
        self.assertEqual(record["mmc"]["refA"], {
            "mmc_properties": os.path.join(mmc, "refA", "properties.yaml"),
            "mmc_annotations": os.path.join(mmc, "refA", f"{uei}.csv"),
            "mmc_json": os.path.join(mmc, "refA", f"{uei}.json"),
            "mmc_cell_type_counts": os.path.join(mmc, "refA", f"{uei}.cell_type_counts.tsv"),
        })
        self.assertEqual(record["mmc"]["refB"]["mmc_properties"], os.path.join(mmc, "refB", "properties.yaml"))
        self.assertEqual(record["mmc"]["refB"]["mmc_annotations"], NA)

    def test_reference_and_reduced_gtf(self):
        record = locator.locate_scRNA_artifacts(self.make_experiment("exp"))
        self.assertEqual(record["reference"], os.path.join(self.reference_dir(), "GRCh38.fasta.gz"))
        self.assertEqual(record["reduced_gtf"], os.path.join(self.reference_dir(), "GRCh38.reduced.gtf"))

    def test_reference_without_reduced_gtf(self):
        record = locator.locate_scRNA_artifacts(self.make_experiment("exp", reference_files=("GRCh38.fasta.gz",)))
        self.assertEqual(record["reference"], os.path.join(self.reference_dir(), "GRCh38.fasta.gz"))
        self.assertEqual(record["reduced_gtf"], NA)

    def test_more_than_one_reduced_gtf(self):
        dge = self.make_experiment("exp", reference_files=("GRCh38.fasta.gz", "a.reduced.gtf", "b.reduced.gtf"))
        with self.assertRaises(ValueError):
            locator.locate_scRNA_artifacts(dge)

    def test_no_reference_in_alignment_properties(self):
        record = locator.locate_scRNA_artifacts(self.make_experiment("exp", reference_files=None))
        self.assertEqual((record["reference"], record["reduced_gtf"]), (NA, NA))

    def test_no_mmc(self):
        record = locator.locate_scRNA_artifacts(self.make_experiment("exp", mmc_references=()))
        self.assertEqual(record["mmc"], NA)

    def test_cbrb_files(self):
        uei = "exp"
        record = locator.locate_scRNA_artifacts(self.make_experiment(uei))
        cbrb = os.path.join(self.tmpDir, uei + "_output", CBRB_REL)
        self.assertEqual(record["cbrb_report"], os.path.join(cbrb, f"{uei}_report.html"))
        self.assertEqual(record["cbrb_tearsheet"], os.path.join(cbrb, f"{uei}.cbrb_tearsheet.pdf"))
        self.assertEqual(record["cbrb_parameter_estimation_pdf"],
                         os.path.join(cbrb, f"{uei}.svm_cbrb_parameter_estimation.pdf"))

    def test_dge_manifest_order(self):
        dge_b = self.make_experiment("B")
        dge_a = self.make_experiment("A")
        dges = locator.read_dge_manifest(io.StringIO(yaml.safe_dump({"dges": [{"dge": dge_b}, {"dge": dge_a}]})))
        self.assertEqual(dges, [dge_b, dge_a])
        result = locator.locate_datasets(dges)
        self.assertEqual([d["uei"] for d in result["datasets"]], ["B", "A"])
        self.assertEqual([d["dge_selected_cells"] for d in result["datasets"]], [dge_b, dge_a])

    def assertRoundTrip(self, result):
        """
        Serialize result, parse it back, and check the data structure (values, field order, list values) is
        recovered exactly.
        """
        out = io.StringIO()
        locator.write_artifact_manifest(result, out)
        text = out.getvalue()
        self.assertTrue(text.startswith("datasets:\n- uei: "), text[:40])
        self.assertNotIn("'", text)
        self.assertNotIn('"', text)
        datasets = locator.load_artifact_manifest(io.StringIO(text))
        self.assertEqual(datasets, result["datasets"])
        for parsed, original in zip(datasets, result["datasets"]):
            self.assertKeyOrderEqual(parsed, original)
        return text

    def assertKeyOrderEqual(self, parsed, original):
        self.assertEqual(list(parsed.keys()), list(original.keys()))
        for key, value in original.items():
            if isinstance(value, dict):
                self.assertKeyOrderEqual(parsed[key], value)

    def test_yaml_round_trip(self):
        std = os.path.dirname(self.make_experiment("A", mmc_references=("refA", "refB")))
        result = locator.locate_datasets([
            # village DGE, two MMC references
            os.path.join(std, "village", VILLAGE_NAME, "A.donors.digital_expression.txt.gz"),
            # selected DGE, no village or MMC
            self.make_experiment("B", mmc_references=(), villages=()),
            # non-standard DGE type
            os.path.join(std, "A.metagene.digital_expression.txt.gz"),
        ])
        text = self.assertRoundTrip(result)
        self.assertIn("  mmc: NA\n", text)
        self.assertIn("  mmc:\n    refA:\n      mmc_properties: ", text)
        self.assertIn("dge_donors: NA\n", text)
        # each DGE line is immediately followed by its summary
        lines = text.splitlines()
        for i, line in enumerate(lines):
            if line.strip().startswith("dge_") and not line.strip().startswith("dge_summary_"):
                dge_type = line.split(":")[0].strip()[len("dge_"):]
                self.assertTrue(lines[i + 1].strip().startswith(f"dge_summary_{dge_type}:"), lines[i:i + 2])
        self.assertEqual(text.count("- uei: "), 3)

    def test_output_dir(self):
        result = locator.locate_datasets([self.make_experiment("A"), self.make_experiment("B")])
        output_dir = os.path.join(self.tmpDir, "manifests")
        locator.write_artifact_manifest_dir(result, output_dir)
        self.assertEqual(sorted(os.listdir(output_dir)), ["A.yaml", "B.yaml"])
        self.assertEqual(locator.load_artifact_manifest(os.path.join(output_dir, "B.yaml")), [result["datasets"][1]])

    def test_output_dir_duplicate_uei(self):
        dge = self.make_experiment("A")
        result = locator.locate_datasets([dge, dge])
        output_dir = os.path.join(self.tmpDir, "manifests")
        with self.assertRaises(ValueError):
            locator.write_artifact_manifest_dir(result, output_dir)
        self.assertFalse(os.path.exists(output_dir))

    def test_main_dge(self):
        dge = self.make_experiment("A")
        output = os.path.join(self.tmpDir, "out.yaml")
        self.assertEqual(locator.main(["--dge", dge, "--output", output]), 0)
        datasets = locator.load_artifact_manifest(output)
        self.assertEqual(len(datasets), 1)
        self.assertEqual(datasets[0]["dge_selected_cells"], dge)
        self.assertEqual(datasets, locator.locate_datasets([dge])["datasets"])

    def test_main_manifest_output_dir(self):
        manifest = os.path.join(self.tmpDir, "dges.yaml")
        with open(manifest, "w") as f:
            yaml.safe_dump({"dges": [{"dge": self.make_experiment("A")}, {"dge": self.make_experiment("B")}]}, f)
        output_dir = os.path.join(self.tmpDir, "manifests")
        self.assertEqual(locator.main(["--manifest", manifest, "--output-dir", output_dir]), 0)
        self.assertEqual(sorted(os.listdir(output_dir)), ["A.yaml", "B.yaml"])

    def test_user_dge_and_donor_default_to_na(self):
        dge = self.make_experiment("A")
        record = locator.locate_scRNA_artifacts(dge)
        self.assertEqual(record["user_dge"], dge)
        self.assertEqual(record["donor"], NA)

    def test_user_dge_is_the_input_dge_of_any_type(self):
        selected = self.make_experiment("A")
        village = os.path.join(os.path.dirname(selected), "village", VILLAGE_NAME, "A.donors.digital_expression.txt.gz")
        self.assertEqual(locator.locate_scRNA_artifacts(village)["user_dge"], village)
        self.assertEqual(locator.locate_scRNA_artifacts(selected)["user_dge"], selected)

    def test_donor(self):
        record = locator.locate_scRNA_artifacts(self.make_experiment("A"), donor="N1")
        self.assertEqual(record["donor"], "N1")
        # a donor that YAML parses as a number is still a string
        self.assertEqual(locator.locate_scRNA_artifacts(self.make_experiment("B"), donor=12)["donor"], "12")

    def test_read_dge_entries_applies_defaults(self):
        manifest = {"dgeDefaults": {"donor": "D0", "filters": {"x": {"min": 1}}},
                    "dges": [{"dge": "a"}, {"dge": "b", "donor": "D1"}]}
        entries = locator.read_dge_entries(io.StringIO(yaml.safe_dump(manifest)))
        self.assertEqual([e["donor"] for e in entries], ["D0", "D1"])
        self.assertEqual(locator.read_dge_manifest(io.StringIO(yaml.safe_dump(manifest))), ["a", "b"])
        entries[0]["filters"]["x"]["min"] = 2
        self.assertEqual(entries[1]["filters"]["x"]["min"], 1)

    def test_read_dge_entries_dge_in_defaults(self):
        with self.assertRaises(ValueError):
            locator.read_dge_entries(io.StringIO(yaml.safe_dump({"dgeDefaults": {"dge": "a"}, "dges": [{"dge": "b"}]})))

    def test_locate_datasets_donors(self):
        dges = [self.make_experiment("A"), self.make_experiment("B")]
        result = locator.locate_datasets(dges, donors=["N1", None])
        self.assertEqual([d["donor"] for d in result["datasets"]], ["N1", NA])
        with self.assertRaises(ValueError):
            locator.locate_datasets(dges, donors=["N1"])

    def test_main_manifest_donors(self):
        manifest = os.path.join(self.tmpDir, "dges.yaml")
        with open(manifest, "w") as f:
            yaml.safe_dump({"dgeDefaults": {"donor": "D0"},
                            "dges": [{"dge": self.make_experiment("A")}, {"dge": self.make_experiment("B"), "donor": 7}]}, f)
        output = os.path.join(self.tmpDir, "out.yaml")
        self.assertEqual(locator.main(["--manifest", manifest, "--output", output]), 0)
        self.assertEqual([d["donor"] for d in locator.load_artifact_manifest(output)], ["D0", "7"])

    def test_main_dge_donor(self):
        output = os.path.join(self.tmpDir, "out.yaml")
        self.assertEqual(locator.main(["--dge", self.make_experiment("A"), "--donor", "N1", "--output", output]), 0)
        self.assertEqual(locator.load_artifact_manifest(output)[0]["donor"], "N1")

    def test_main_donor_requires_dge(self):
        manifest = os.path.join(self.tmpDir, "dges.yaml")
        with open(manifest, "w") as f:
            yaml.safe_dump({"dges": [{"dge": self.make_experiment("A")}]}, f)
        with self.assertRaises(SystemExit):
            locator.main(["--manifest", manifest, "--donor", "N1"])

    def test_main_duplicate_uei(self):
        dge = self.make_experiment("A")
        village = os.path.join(os.path.dirname(dge), "village", VILLAGE_NAME, "A.donors.digital_expression.txt.gz")
        manifest = os.path.join(self.tmpDir, "dges.yaml")
        with open(manifest, "w") as f:
            yaml.safe_dump({"dges": [{"dge": dge}, {"dge": village}]}, f)
        output = os.path.join(self.tmpDir, "out.yaml")
        with self.assertRaises(ValueError):
            locator.main(["--manifest", manifest, "--output", output])
        # nothing was written to the output
        self.assertEqual(os.path.getsize(output), 0)


if __name__ == '__main__':
    unittest.main()
