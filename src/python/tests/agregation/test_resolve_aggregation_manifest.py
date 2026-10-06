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

import yaml

from dropseq.aggregation import resolve_aggregation_manifest as ram

ROOT = "gs://bucket/exp"
REDUCED_GTF = "gs://ref/GRCh38.reduced.gtf"


def dge_url(uei, kind="donors"):
    return f"{ROOT}/{uei}/GRCh38/cbrb/auto/cell_selection/auto/standard_analysis/{uei}.{kind}.digital_expression.txt.gz"


class TestResolveAggregationManifest(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(".tmp", "resolve_aggregation_manifest.")
        self.artifactsDir = os.path.join(self.tmpDir, "artifacts")
        os.makedirs(self.artifactsDir)

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def write_artifacts(self, uei, reference="GRCh38", donors=True, user_kind="donors", reference_url=None,
                        donor=None, mmc=("HMBA", "other"), reduced_gtf=REDUCED_GTF):
        dataset = {"uei": uei, "user_dge": dge_url(uei, user_kind), "alignment_dir": f"{ROOT}/{uei}/{reference}",
                   "reference": reference_url or "NA", "reduced_gtf": reduced_gtf,
                   "cell_metadata": f"{ROOT}/{uei}/{uei}.cmd.tsv", "donor": donor or "NA",
                   "mmc": {m: {"mmc_annotations": f"{ROOT}/{uei}/mmc/{m}/{uei}.csv"} for m in mmc} or "NA",
                   "dge_unfiltered": f"{ROOT}/{uei}/{uei}.digital_expression.txt.gz",
                   "dge_summary_unfiltered": f"{ROOT}/{uei}/{uei}.digital_expression_summary.txt",
                   "dge_selected_cells": dge_url(uei, "selected"),
                   "dge_donors": dge_url(uei, "donors") if donors else "NA"}
        with open(os.path.join(self.artifactsDir, uei + ".yaml"), "w") as f:
            yaml.safe_dump({"datasets": [dataset]}, f)

    def test_defaults_are_projected_unless_set(self):
        manifest = {"dgeDefaults": {"donor": "D0", "filters": {"doublet": {"exclude": "doublet"}}},
                    "dges": [{"dge": "a"}, {"dge": "b", "donor": "D1"}]}
        dges = ram.resolve_dges(manifest)
        self.assertEqual([d["donor"] for d in dges], ["D0", "D1"])
        self.assertTrue(all("doublet" in d["filters"] for d in dges))
        # the defaults must not be shared between dges
        dges[0]["filters"]["x"] = {}
        self.assertNotIn("x", dges[1]["filters"])

    def test_single_dge_dictionary(self):
        self.assertEqual(len(ram.resolve_dges({"dges": {"dge": "a"}})), 1)

    def test_no_dges(self):
        with self.assertRaises(ValueError):
            ram.resolve_dges({"dges": []})

    def test_unknown_keys(self):
        with self.assertRaises(ValueError):
            ram.resolve_dges({"dges": [{"dge": "a", "scPred": "x"}]})
        with self.assertRaises(ValueError):
            ram.resolve_dges({"dgeDefaults": {"dge": "a"}, "dges": [{"dge": "b"}]})

    def test_missing_dge(self):
        with self.assertRaises(ValueError):
            ram.resolve_dges({"dges": [{"libraryId": "x"}]})

    def test_args_empty(self):
        self.assertEqual(ram.join_and_filter_args({"dge": "a"}), [])

    def test_args_donor_is_not_an_argument(self):
        # the donor is recorded in the artifacts by locate_scRNA_artifacts
        self.assertEqual(ram.join_and_filter_args({"dge": "a", "donor": 12}), [])

    def test_args_filters(self):
        dge = {"dge": "a", "filters": {
            "doublet": {"exclude": "doublet"},
            "predClass": {"include": ["a", "b"], "includeFile": "in.txt", "excludeFile": "out.txt"},
            "pct_mt": {"min": 0.01, "max": 0.1}}}
        self.assertEqual(ram.join_and_filter_args(dge), [
            "--exclude", "doublet", "doublet",
            "--include-file", "predClass", "in.txt", "--exclude-file", "predClass", "out.txt",
            "--include", "predClass", "a", "b",
            "--min", "pct_mt", "0.01", "--max", "pct_mt", "0.1"])

    def test_args_filters_as_list(self):
        dge = {"dge": "a", "filters": [{"x": {"min": 1}}, {"y": {"max": 2}}]}
        self.assertEqual(ram.join_and_filter_args(dge), ["--min", "x", "1", "--max", "y", "2"])
        with self.assertRaises(ValueError):
            ram.join_and_filter_args({"dge": "a", "filters": [{"x": {"min": 1}}, {"x": {"max": 2}}]})

    def test_args_unknown_filter_key(self):
        with self.assertRaises(ValueError):
            ram.join_and_filter_args({"dge": "a", "filters": {"x": {"equals": 1}}})

    def test_args_joins(self):
        dge = {"dge": "a", "joins": [{"joinFile": "j1.tsv", "leftColumn": "l1", "joinColumn": "c1"},
                                     {"joinFile": "j2.tsv", "leftColumn": "l2", "joinColumn": "c2"}]}
        self.assertEqual(ram.join_and_filter_args(dge),
                         ["--join", "j1.tsv", "l1", "c1", "--join", "j2.tsv", "l2", "c2"])
        # a single join does not need to be in a list
        self.assertEqual(ram.join_and_filter_args({"dge": "a", "joins": dge["joins"][0]}),
                         ["--join", "j1.tsv", "l1", "c1"])

    def test_args_join_missing_key(self):
        with self.assertRaises(ValueError):
            ram.join_and_filter_args({"dge": "a", "joins": [{"joinFile": "j.tsv", "leftColumn": "l"}]})

    def test_resolve_matches_dges_by_url_not_position(self):
        for uei in ("lib1", "lib2"):
            self.write_artifacts(uei)
        manifest = {"dges": [{"dge": dge_url("lib2")}, {"dge": dge_url("lib1")}]}
        resolved = ram.resolve(manifest, self.artifactsDir, "agg")
        self.assertEqual([r[0] for r in resolved], ["lib2", "lib1"])
        self.assertEqual(os.path.basename(resolved[0][1]), "lib2.yaml")

    def test_resolve_matches_user_dge_not_other_dges_of_the_library(self):
        self.write_artifacts("lib1", user_kind="selected")
        # the donors DGE belongs to lib1 but is not the DGE the user gave
        with self.assertRaises(ValueError):
            ram.resolve({"dges": [{"dge": dge_url("lib1", "donors")}]}, self.artifactsDir, "agg")
        resolved = ram.resolve({"dges": [{"dge": dge_url("lib1", "selected"), "filters": {"x": {"min": 1}}}]},
                               self.artifactsDir, "agg")
        self.assertEqual(resolved[0][0], "lib1")
        self.assertEqual(resolved[0][2], ["--min", "x", "1"])

    def test_library_id_is_not_supported(self):
        with self.assertRaises(ValueError):
            ram.resolve_dges({"dges": [{"dge": "a", "libraryId": "renamed"}]})

    def test_artifacts_without_user_dge(self):
        with open(os.path.join(self.artifactsDir, "old.yaml"), "w") as f:
            yaml.safe_dump({"datasets": [{"uei": "old", "alignment_dir": f"{ROOT}/old/GRCh38"}]}, f)
        with self.assertRaises(ValueError):
            ram.load_artifacts_by_dge(self.artifactsDir)

    def test_resolve_library_without_donors(self):
        self.write_artifacts("lib1", donors=False, user_kind="selected")
        resolved = ram.resolve({"dges": [{"dge": dge_url("lib1", "selected")}]}, self.artifactsDir, "agg")
        self.assertEqual(resolved[0][0], "lib1")

    def test_resolve_dge_not_found(self):
        self.write_artifacts("lib1")
        with self.assertRaises(ValueError):
            ram.resolve({"dges": [{"dge": dge_url("other")}]}, self.artifactsDir, "agg")

    def test_resolve_library_id_is_analysis_id(self):
        self.write_artifacts("lib1")
        with self.assertRaises(ValueError):
            ram.resolve({"dges": [{"dge": dge_url("lib1")}]}, self.artifactsDir, "lib1")

    def test_resolve_different_references(self):
        self.write_artifacts("lib1", reference="GRCh38")
        self.write_artifacts("lib2", reference="mm10")
        manifest = {"dges": [{"dge": dge_url("lib1")}, {"dge": dge_url("lib2")}]}
        with self.assertRaises(ValueError):
            ram.resolve(manifest, self.artifactsDir, "agg")
        self.assertEqual(len(ram.resolve(manifest, self.artifactsDir, "agg", reference="GRCh38")), 2)

    def test_resolve_references_from_reference_url(self):
        # the same alignment directory name, but different reference files
        self.write_artifacts("lib1", reference_url="gs://ref/a/GRCh38.fasta.gz")
        self.write_artifacts("lib2", reference_url="gs://ref/b/GRCh38.fasta.gz")
        manifest = {"dges": [{"dge": dge_url("lib1")}, {"dge": dge_url("lib2")}]}
        with self.assertRaises(ValueError):
            ram.resolve(manifest, self.artifactsDir, "agg")
        self.write_artifacts("lib2", reference_url="gs://ref/a/GRCh38.fasta.gz")
        self.assertEqual(len(ram.resolve(manifest, self.artifactsDir, "agg")), 2)

    def test_resolve_library_fields(self):
        self.write_artifacts("lib1")
        library = ram.resolve({"dges": [{"dge": dge_url("lib1")}]}, self.artifactsDir, "agg")[0]
        self.assertEqual(library.dge, dge_url("lib1", "donors"))
        self.assertEqual(library.cell_metadata, f"{ROOT}/lib1/lib1.cmd.tsv")
        self.assertEqual((library.mmc_model, library.mmc_annotations), ("HMBA", f"{ROOT}/lib1/mmc/HMBA/lib1.csv"))
        self.assertEqual(library.reduced_gtf, REDUCED_GTF)
        self.assertEqual(library.args, [])

    def test_resolve_dge_without_donors_is_selected_cells(self):
        self.write_artifacts("lib1", donors=False, user_kind="selected")
        library = ram.resolve({"dges": [{"dge": dge_url("lib1", "selected")}]}, self.artifactsDir, "agg")[0]
        self.assertEqual(library.dge, dge_url("lib1", "selected"))

    def test_resolve_mmc_model(self):
        self.write_artifacts("lib1")
        manifest = {"dges": [{"dge": dge_url("lib1")}]}
        library = ram.resolve(manifest, self.artifactsDir, "agg", mmc_model="other")[0]
        self.assertEqual((library.mmc_model, library.mmc_annotations), ("other", f"{ROOT}/lib1/mmc/other/lib1.csv"))
        # a model the library does not have
        library = ram.resolve(manifest, self.artifactsDir, "agg", mmc_model="missing")[0]
        self.assertEqual((library.mmc_model, library.mmc_annotations), ("missing", "NA"))

    def test_resolve_no_mmc(self):
        self.write_artifacts("lib1", mmc=())
        library = ram.resolve({"dges": [{"dge": dge_url("lib1")}]}, self.artifactsDir, "agg")[0]
        self.assertEqual((library.mmc_model, library.mmc_annotations), ("NA", "NA"))

    def test_resolve_manifest_donor_is_set_first(self):
        self.write_artifacts("lib1", donor="N1")
        manifest = {"dges": [{"dge": dge_url("lib1"), "filters": {"x": {"min": 1}}}]}
        self.assertEqual(ram.resolve(manifest, self.artifactsDir, "agg")[0].args,
                         ["--set", "donor", "N1", "--min", "x", "1"])

    def test_resolve_requires_cell_metadata(self):
        self.write_artifacts("lib1")
        path = os.path.join(self.artifactsDir, "lib1.yaml")
        with open(path) as f:
            artifacts = yaml.safe_load(f)
        artifacts["datasets"][0]["cell_metadata"] = "NA"
        with open(path, "w") as f:
            yaml.safe_dump(artifacts, f)
        with self.assertRaises(ValueError):
            ram.resolve({"dges": [{"dge": dge_url("lib1")}]}, self.artifactsDir, "agg")

    def test_resolve_reduced_gtf(self):
        self.write_artifacts("lib1", reduced_gtf="NA")
        manifest = {"dges": [{"dge": dge_url("lib1")}]}
        with self.assertRaises(ValueError):
            ram.resolve(manifest, self.artifactsDir, "agg")
        self.write_artifacts("lib1")
        self.write_artifacts("lib2", reduced_gtf="gs://ref/other.reduced.gtf")
        manifest["dges"].append({"dge": dge_url("lib2")})
        with self.assertRaises(ValueError):
            ram.resolve(manifest, self.artifactsDir, "agg")

    def test_main(self):
        for uei in ("lib1", "lib2"):
            self.write_artifacts(uei)
        manifestFile = os.path.join(self.tmpDir, "manifest.yaml")
        with open(manifestFile, "w") as f:
            yaml.safe_dump({"dgeDefaults": {"filters": {"doublet": {"exclude": "doublet"}}},
                            "dges": [{"dge": dge_url("lib1")}, {"dge": dge_url("lib2")}]}, f)
        outDir = os.path.join(self.tmpDir, "out")
        self.assertEqual(ram.main(["--manifest", manifestFile, "--artifacts-dir", self.artifactsDir,
                                   "--analysis-id", "agg", "--output-dir", outDir]), 0)
        with open(os.path.join(outDir, ram.LIBRARIES_FILE)) as f:
            rows = [line.rstrip("\n").split("\t") for line in f]
        self.assertEqual(rows[0], ram.LIBRARIES_COLUMNS)
        rows = [dict(zip(rows[0], row)) for row in rows[1:]]
        self.assertEqual([r["library_id"] for r in rows], ["lib1", "lib2"])
        self.assertEqual(rows[0]["mmc_model"], "HMBA")
        self.assertEqual(rows[0]["reduced_gtf"], REDUCED_GTF)
        with open(rows[0]["args_file"]) as f:
            self.assertEqual(f.read().split("\n")[:-1], ["--exclude", "doublet", "doublet"])

    def test_main_returns_error(self):
        self.write_artifacts("lib1")
        manifestFile = os.path.join(self.tmpDir, "manifest.yaml")
        with open(manifestFile, "w") as f:
            yaml.safe_dump({"dges": [{"dge": dge_url("lib1")}]}, f)
        self.assertEqual(ram.main(["--manifest", manifestFile, "--artifacts-dir", self.artifactsDir,
                                   "--analysis-id", "lib1", "--output-dir", os.path.join(self.tmpDir, "o")]), 1)


if __name__ == "__main__":
    unittest.main()
