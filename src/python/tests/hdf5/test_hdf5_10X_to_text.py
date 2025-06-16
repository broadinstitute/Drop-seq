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
import filecmp
import os.path
import shutil
import tempfile
import unittest
import collections
import dropseq.hdf5.hdf5_10X_to_text

OptionsTuple = collections.namedtuple("OptionsTuple", ["input", "output", "output_sizes",
                                                       "progress_interval", "analyzed_barcodes_only",
                                                       "limit", "header", "cbrb_log"],
                                      defaults=(None, 1000, False, 2000, None, None))
class TestHdf5_10X_to_text(unittest.TestCase):
    def setUp(self):
        self.testDataDir = "tests/data/hdf5/hdf5_10X_to_text"
        self.inputFile = os.path.join(self.testDataDir, "N701.h5")
        self.tmpDir = tempfile.mkdtemp(".tmp", "hdf5_10X_to_text.")
        self.sizesFile = os.path.join(self.tmpDir, "N701.digtal_expression_sizes.txt")
        self.options = OptionsTuple(self.inputFile,
                                    os.path.join(self.tmpDir, "N701.digital_expression.txt"),
                                    open(self.sizesFile, "w"))
    def tearDown(self):
        shutil.rmtree(self.tmpDir)
    def test_basic(self):
        self.assertEqual(dropseq.hdf5.hdf5_10X_to_text.run(self.options), 0)
        self.assertTrue(filecmp.cmp(self.options.output, os.path.join(self.testDataDir, "N701.2000.digital_expression.txt"), shallow=False))
        self.assertTrue(filecmp.cmp(self.sizesFile, os.path.join(self.testDataDir, "N701.2000.sizes.txt"), shallow=False))

    def test_analyzed_barcodes_only(self):
        options = self.options._replace(analyzed_barcodes_only=True)
        self.assertEqual(dropseq.hdf5.hdf5_10X_to_text.run(options), 0)
        self.assertTrue(filecmp.cmp(self.options.output, os.path.join(self.testDataDir, "N701.2000.analyzed_barcodes_only.digital_expression.txt"), shallow=False))
        self.assertTrue(filecmp.cmp(self.sizesFile, os.path.join(self.testDataDir, "N701.2000.analyzed_barcodes_only.sizes.txt"), shallow=False))

if __name__ == '__main__':
    unittest.main()
