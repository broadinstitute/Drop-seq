#!/usr/bin/env python3
# MIT License
# 
# Copyright 2025 Broad Institute
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
import csv
import os
import shutil
import tempfile
import unittest

import dropseq.metadata.read_gtf as read_gtf


class TestReadGTF(unittest.TestCase):
    def setUp(self):
        self.testDataDir = "tests/data/metadata"
        self.tmpDir = tempfile.mkdtemp(".tmp", "read_gtf.")
        self.outputFile = os.path.join(self.tmpDir, "output.gtf")

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def test_basic(self):
        df = read_gtf.read_gtf(os.path.join(self.testDataDir, "GRCh38-2020-A.100.gtf"))
        df.to_csv(self.outputFile, sep="\t", index=False, quotechar='', quoting=csv.QUOTE_NONE)
        with open(self.outputFile, 'r') as f1, open(os.path.join(self.testDataDir, "GRCh38-2020-A.100.parsed.gtf"), 'r') as f2:
            content1 = f1.read()
            content2 = f2.read()

        self.assertEqual(content1, content2, "Files are not equal")



