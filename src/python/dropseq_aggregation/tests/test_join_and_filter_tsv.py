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
import logging

import pandas.testing

import dropseq_aggregation.join_and_filter_tsv

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
OptionsTuple = collections.namedtuple("OptionsTuple", ["output", "input", "join", "set", "min",
                                                       "max", "include_file", "exclude_file", "include", "exclude"],
                                      defaults=(None, [], [], [], [], [], [], [], []))

class TestJoinAndFilterTSV(unittest.TestCase):
    def setUp(self):
        self.testDataDir = "../../../testdata/python/dropseq_aggregation/join_and_filter_tsv"
        self.tmpDir = tempfile.mkdtemp(".tmp", "join_and_filter_tsv.")
        self.outputFile = os.path.join(self.tmpDir, "output.tsv")
        self.options = OptionsTuple(open(self.outputFile, "w"))

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def test_basic(self):
        primary = os.path.join(self.testDataDir, "sample1.100.cell_metadata.txt")
        secondary = os.path.join(self.testDataDir, "sample1.100.scPred.txt")
        options = self.options._replace(input=open(primary),
                                        join=[(secondary, "CELL_BARCODE", "CELL_BARCODE")])
        self.assertEqual(dropseq_aggregation.join_and_filter_tsv.main(options), 0)
        self.assertSharedColumnsEqual(self.outputFile, primary)
        self.assertSharedColumnsEqual(self.outputFile, secondary)

    def test_fewer_secondary(self):
        primary = os.path.join(self.testDataDir, "sample1.100.cell_metadata.txt")
        secondary = os.path.join(self.testDataDir, "sample1.50.scPred.txt")
        options = self.options._replace(input=open(primary),
                                        join=[(secondary, "CELL_BARCODE", "CELL_BARCODE")])
        self.assertEqual(dropseq_aggregation.join_and_filter_tsv.main(options), 0)
        self.assertSharedColumnsEqual(self.outputFile, primary)
        self.assertSharedColumnsEqual(self.outputFile, secondary, wideRows=49)

    def test_fewer_primary(self):
        primary = os.path.join(self.testDataDir, "sample1.50.cell_metadata.txt")
        secondary = os.path.join(self.testDataDir, "sample1.100.scPred.txt")
        options = self.options._replace(input=open(primary),
                                        join=[(secondary, "CELL_BARCODE", "CELL_BARCODE")])
        self.assertEqual(dropseq_aggregation.join_and_filter_tsv.main(options), 0)
        self.assertSharedColumnsEqual(self.outputFile, primary)
        self.assertSharedColumnsEqual(self.outputFile, secondary, narrowRows=49)

    def test_additional_join(self):
        primary = os.path.join(self.testDataDir, "sample1.100.cell_metadata.txt")
        secondary1 = os.path.join(self.testDataDir, "sample1.100.scPred.txt")
        secondary2 = os.path.join(self.testDataDir, "sample1.donor_age.txt")
        options = self.options._replace(input=open(primary),
                                        join=[(secondary1, "CELL_BARCODE", "CELL_BARCODE"),
                                              (secondary2, "DONOR", "DONOR")])
        self.assertEqual(dropseq_aggregation.join_and_filter_tsv.main(options), 0)
        self.assertSharedColumnsEqual(self.outputFile, primary)
        self.assertSharedColumnsEqual(self.outputFile, secondary1)
        self.assertMultiJoin(self.outputFile, secondary2, "DONOR", "DONOR")

    def test_set(self):
        primary = os.path.join(self.testDataDir, "sample1.100.cell_metadata.txt")
        secondary = os.path.join(self.testDataDir, "sample1.100.scPred.txt")
        setTuples = [("PREFIX", "library1"), ("SUFFIX", "brary1")]
        options = self.options._replace(input=open(primary),
                                        join=[(secondary, "CELL_BARCODE", "CELL_BARCODE")],
                                        set=setTuples)
        self.assertEqual(dropseq_aggregation.join_and_filter_tsv.main(options), 0)
        outputDf = pd.read_csv(self.outputFile, sep='\t')
        for column, value in setTuples:
            self.assertTrue((outputDf[column] == value).all())
        self.assertSharedColumnsEqual(self.outputFile, primary, dropColumns = [column for column, _ in setTuples])
        self.assertSharedColumnsEqual(self.outputFile, secondary)

    def test_min(self):
        primary = os.path.join(self.testDataDir, "sample1.100.cell_metadata.txt")
        secondary = os.path.join(self.testDataDir, "sample1.100.scPred.txt")
        options = self.options._replace(input=open(primary),
                                        join=[(secondary, "CELL_BARCODE", "CELL_BARCODE")],
                                        min=[("max.prob", "0.8")])
        self.assertEqual(dropseq_aggregation.join_and_filter_tsv.main(options), 0)
        outputDf = pd.read_csv(self.outputFile, sep='\t')
        self.assertTrue((outputDf["max.prob"] >= 0.8).all())
        primaryDf = pd.read_csv(primary, sep='\t')
        self.assertEqual(len(primaryDf[primaryDf["max.prob"] >= 0.8]), len(outputDf))

    def test_max(self):
        primary = os.path.join(self.testDataDir, "sample1.100.cell_metadata.txt")
        secondary = os.path.join(self.testDataDir, "sample1.100.scPred.txt")
        options = self.options._replace(input=open(primary),
                                        join=[(secondary, "CELL_BARCODE", "CELL_BARCODE")],
                                        max=[("max.prob", "0.8")])
        self.assertEqual(dropseq_aggregation.join_and_filter_tsv.main(options), 0)
        outputDf = pd.read_csv(self.outputFile, sep='\t')
        self.assertTrue((outputDf["max.prob"] <= 0.8).all())
        primaryDf = pd.read_csv(primary, sep='\t')
        self.assertEqual(len(primaryDf[primaryDf["max.prob"] <= 0.8]), len(outputDf))

    def test_include_file(self):
        primary = os.path.join(self.testDataDir, "sample1.100.cell_metadata.txt")
        includeFile = os.path.join(self.testDataDir, "donor_subset.txt")
        options = self.options._replace(input=open(primary),
                                        include_file=[("DONOR", includeFile)])
        self.assertEqual(dropseq_aggregation.join_and_filter_tsv.main(options), 0)
        outputDf = pd.read_csv(self.outputFile, sep='\t')
        includeValues = pd.read_csv(includeFile, sep='\t', header=None).iloc[0]
        self.assertTrue((outputDf["DONOR"].isin(includeValues)).all())

    def test_exclude_file(self):
        primary = os.path.join(self.testDataDir, "sample1.100.cell_metadata.txt")
        excludeFile = os.path.join(self.testDataDir, "donor_subset.txt")
        options = self.options._replace(input=open(primary),
                                        exclude_file=[("DONOR", excludeFile)])
        self.assertEqual(dropseq_aggregation.join_and_filter_tsv.main(options), 0)
        outputDf = pd.read_csv(self.outputFile, sep='\t')
        excludeValues = pd.read_csv(excludeFile, sep='\t', header=None).iloc[0]
        self.assertFalse((outputDf["DONOR"].isin(excludeValues)).any())

    def test_include_exclude(self):
        primary = os.path.join(self.testDataDir, "sample1.100.cell_metadata.txt")
        donorsToInclude = ["donor1", "donor2", "donor3"]
        predClassesToExclude = ["gabaergic"]
        options = self.options._replace(input=open(primary),
                                        include=[["DONOR"] + donorsToInclude],
                                        exclude=[["predClass"] + predClassesToExclude])
        self.assertEqual(dropseq_aggregation.join_and_filter_tsv.main(options), 0)
        outputDf = pd.read_csv(self.outputFile, sep='\t')
        self.assertTrue((outputDf["DONOR"].isin(donorsToInclude)).all())
        self.assertFalse((outputDf["predClass"].isin(predClassesToExclude)).any())

    def test_negative_non_unique_join(self):
        primary = os.path.join(self.testDataDir, "sample1.100.cell_metadata.txt")
        secondary = os.path.join(self.testDataDir, "sample1.nonunique.scPred.txt")
        options = self.options._replace(input=open(primary),
                                        join=[(secondary, "CELL_BARCODE", "CELL_BARCODE")])
        self.assertEqual(dropseq_aggregation.join_and_filter_tsv.main(options), 1)

    def test_boolean(self):
        primary = os.path.join(self.testDataDir, "sample1.100.cell_metadata.txt")
        options = self.options._replace(input=open(primary),
                                        exclude=[["doublet", "true"]])
        self.assertEqual(dropseq_aggregation.join_and_filter_tsv.main(options), 0)
        outputDf = pd.read_csv(self.outputFile, sep='\t')
        self.assertFalse(outputDf["doublet"].any())


    def assertSharedColumnsEqual(self, wideFile, narrowFile, wideRows = None, narrowRows = None, dropColumns = None):
        wideDf = pd.read_csv(wideFile, sep='\t', index_col=False)
        if wideRows is not None:
            wideDf = wideDf.head(wideRows)
        narrowDf = pd.read_csv(narrowFile, sep='\t', index_col=False)
        if narrowRows is not None:
            narrowDf = narrowDf.head(narrowRows)
        # drop dropColumns
        if dropColumns is not None and len(dropColumns) > 0:
            dropColumns = [col for col in dropColumns if col in narrowDf.columns]
            narrowDf = narrowDf.drop(dropColumns, axis=1)
        # drop columns in wideDf that are not in narrowDf
        wideDf = wideDf.drop([col for col in wideDf.columns if col not in narrowDf.columns], axis=1)
        # make column order the same between the two dataframes
        wideDf = wideDf[narrowDf.columns]
        # Especially with subsetting, dtypes can differ, e.g. float64 vs. int64
        pandas.testing.assert_frame_equal(wideDf, narrowDf, check_dtype=False)

    def assertMultiJoin(self, wideFile, narrowFile, wideColumn, narrowColumn):
        wideDf = pd.read_csv(wideFile, sep='\t')
        narrowDf = pd.read_csv(narrowFile, sep='\t')
        # for each row in narrowDf, the value of rows in wideDf should match those in narrowDf, join column excluded
        otherColumns = narrowDf.columns[narrowDf.columns != narrowColumn]
        for index, row in narrowDf.iterrows():
            narrowValue = row[narrowColumn]
            matchingRows = wideDf[wideDf[wideColumn] == narrowValue]
            for index2, wideRow in matchingRows.iterrows():
                for column in otherColumns:
                    self.assertEqual(row[column], wideRow[column])
