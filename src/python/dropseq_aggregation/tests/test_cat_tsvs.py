#!/usr/bin/env python3
# MIT License
# 
# Copyright 2024 Broad Institute
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
import collections
import os
import shutil
import tempfile
import unittest
import pandas as pd
import dropseq_aggregation.cat_tsvs

OptionsTuple = collections.namedtuple("OptionsTuple", ["output", "input", "index_col"],
                                      defaults=[None, None])

class TestCatTsvs(unittest.TestCase):
    def setUp(self):
        self.testDataDir = "../../../testdata/python/dropseq_aggregation/cat_tsvs"
        self.tmpDir = tempfile.mkdtemp(".tmp", "cat_tsvs.")
        self.outputFile = os.path.join(self.tmpDir, "output.tsv")
        self.options = OptionsTuple(open(self.outputFile, "w"))
        self.inputs = [os.path.join(self.testDataDir, f"rxn{i+1}.joined_filtered_cell_metadata.tsv") for i in range(2)]
        self.index_cols = ["PREFIX", "CELL_BARCODE"]

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def test_basic(self):
        options = self.options._replace(index_col=self.index_cols, input=[open(f) for f in self.inputs])
        self.assertEqual(dropseq_aggregation.cat_tsvs.run(options), 0)
        outDf = pd.read_csv(self.outputFile, sep="\t")
        inDfs = [pd.read_csv(f, sep="\t") for f in self.inputs]
        self.assertEqual(len(outDf), sum(len(df) for df in inDfs))
        self.assertEqual(set(outDf.columns), set.union(*[set(df.columns) for df in inDfs]))

    def test_duplicate_keys(self):
        options = self.options._replace(index_col=self.index_cols, input=[open(self.inputs[0]), open(self.inputs[0])])
        self.assertNotEqual(dropseq_aggregation.cat_tsvs.run(options), 0)

    def test_fewer_columns(self):
        dfToClip = pd.read_csv(self.inputs[0], sep="\t", index_col=self.index_cols)
        dfToClip = dfToClip.drop(columns="doublet", axis=1)
        clippedFile = os.path.join(self.tmpDir, "clipped.tsv")
        dfToClip.to_csv(clippedFile, sep="\t")
        inputs = [clippedFile, self.inputs[1]]
        options = self.options._replace(index_col=self.index_cols, input=[open(f) for f in inputs])
        self.assertEqual(dropseq_aggregation.cat_tsvs.run(options), 0)
        outDf = pd.read_csv(self.outputFile, sep="\t")
        inDfs = [pd.read_csv(f, sep="\t") for f in inputs]
        self.assertEqual(len(outDf), sum(len(df) for df in inDfs))
        self.assertEqual(set(outDf.columns), set.union(*[set(df.columns) for df in inDfs]))

    def test_conflicting_column_type(self):
        inputs = [self.inputs[0], os.path.join(self. testDataDir, "conflicting_column_type.joined_filtered_cell_metadata.tsv")]
        options = self.options._replace(index_col=self.index_cols, input=[open(f) for f in inputs])
        self.assertRaisesRegex(Exception, 'Column types disagree', dropseq_aggregation.cat_tsvs.run, options)
