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
import gzip
import os
import shutil
import tempfile
import unittest

import anndata as ad
import pandas as pd

from dropseq.aggregation import build_library_h5ad as bl

DGE = """#DGE\tVERSION:1.1\tEXPRESSION_FORMAT:raw
#COMMAND\tCL:Something
GENE\tAAAA\tCCCC\tGGGG
GA\t1\t0\t3
GB\t0\t0\t0
GC\t2\t5\t0
"""
GTF = ("chr\tstart\tend\tstrand\tgene_name\tgene_id\ttranscript_name\ttranscript_id\ttranscriptType\tannotationType\n"
       "chr1\t1\t2\t+\tGA\tENSG1\tNA\tNA\tNA\tgene\nchr1\t1\t2\t+\tGA\tENSG1\tT1\tT1\tNA\texon\n"
       "chr1\t3\t4\t+\tGB\tENSG2\tNA\tNA\tNA\tgene\nchr1\t5\t6\t+\tGC\tENSG3\tNA\tNA\tNA\tgene\n")
# GGGG first, CCCC absent (filtered out); a missing cell type and a missing score
METADATA = ("PREFIX\tcell_barcode\tdonor\tcell_type\tscore\n"
            "L1\tGGGG\t1\tglia\t0.5\nL1\tAAAA\t2\tNA\tNA\n")


class TestBuildLibraryH5ad(unittest.TestCase):
    def setUp(self):
        self.tmpDir = tempfile.mkdtemp(".tmp", "build_library_h5ad.")
        self.dge = os.path.join(self.tmpDir, "dge.txt.gz")
        with gzip.open(self.dge, "wt") as f:
            f.write(DGE)
        self.metadata = self.write("cmd.txt", METADATA)

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def write(self, name, text):
        path = os.path.join(self.tmpDir, name)
        with open(path, "w") as f:
            f.write(text)
        return path

    def test_round_trip(self):
        out = os.path.join(self.tmpDir, "out.h5ad")
        self.assertEqual(bl.main(["--dge", self.dge, "--cell-metadata", self.metadata,
                                  "--output", out]), 0)
        adata = ad.read_h5ad(out)
        self.assertEqual(list(adata.obs_names), ["L1_GGGG", "L1_AAAA"])
        # the bare barcode is preserved, and the index is <PREFIX>_<barcode>
        self.assertEqual(list(adata.obs["cell_barcode"].astype(str)), ["GGGG", "AAAA"])
        self.assertEqual(list(adata.obs["PREFIX"].astype(str)), ["L1", "L1"])
        self.assertEqual(list(adata.var_names), ["GA", "GB", "GC"])
        self.assertEqual(list(adata.var.columns), ["gene_symbol"])
        self.assertEqual(list(adata.var["gene_symbol"]), ["GA", "GB", "GC"])
        self.assertEqual(adata.X.toarray().tolist(), [[3, 0, 0], [1, 0, 2]])
        # donors are labels, missing values are encoded rather than lost or written as NaN strings
        self.assertEqual(list(adata.obs["donor"].astype(str)), ["1", "2"])
        self.assertEqual(list(adata.obs["cell_type"].astype(str)), ["glia", "NA"])
        self.assertTrue(pd.isna(adata.obs["score"].iloc[1]))

    def test_gene_ids_from_reduced_gtf(self):
        gtf = self.write("ref.gtf", GTF)
        out = os.path.join(self.tmpDir, "out.h5ad")
        bl.main(["--dge", self.dge, "--cell-metadata", self.metadata, "--reduced-gtf", gtf, "--output", out])
        adata = ad.read_h5ad(out)
        self.assertEqual(list(adata.var_names), ["ENSG1", "ENSG2", "ENSG3"])
        self.assertEqual(adata.var.index.name, "gene_ids")
        self.assertEqual(list(adata.var.columns), ["gene_symbol"])
        self.assertEqual(list(adata.var["gene_symbol"]), ["GA", "GB", "GC"])
        self.assertEqual(adata.X.toarray().tolist(), [[3, 0, 0], [1, 0, 2]])

    def test_gene_not_in_reduced_gtf(self):
        gtf = self.write("short.gtf", "\n".join(GTF.splitlines()[:3]) + "\n")
        metadata = pd.DataFrame({"PREFIX": ["L1"], "cell_barcode": ["AAAA"]})
        with self.assertRaisesRegex(ValueError, "not in"):
            bl.build_library_h5ad(self.dge, metadata, gtf)

    def test_cell_not_in_dge(self):
        metadata = pd.DataFrame({"PREFIX": ["L1"], "cell_barcode": ["TTTT"]})
        with self.assertRaisesRegex(ValueError, "not in"):
            bl.build_library_h5ad(self.dge, metadata)


if __name__ == "__main__":
    unittest.main()
