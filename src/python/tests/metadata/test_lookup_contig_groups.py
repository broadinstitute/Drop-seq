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
import collections
import shutil
import tempfile
import unittest
from pathlib import Path

import dropseq.metadata.lookup_contig_groups

OptionsTuple = collections.namedtuple('OptionsTuple', ['contig_groups', 'group', 'out'], defaults=[None])

class TestLookupContigGroups(unittest.TestCase):
    def setUp(self):
        self.tmpDir = Path(tempfile.mkdtemp(".tmp", "lookup_contig_groups."))
        self.outputFile = self.tmpDir / "output.txt"
        self.contigGroupsFile = "tests/data/metadata/contig_groups.yaml"
        self.options = OptionsTuple(
            contig_groups=self.contigGroupsFile,
            group="test_project",
            out=open(self.outputFile, 'w')
        )

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def test_mt(self):
        # Test group that returns a single contig
        options = self.options._replace(group=["MT"])
        self.assertEqual(dropseq.metadata.lookup_contig_groups.run(options), 0)
        with open(self.outputFile, 'r') as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 1)
            self.assertEqual(lines[0].strip(), "chrM")

    def test_non_autosomes(self):
        # Test group that returns multiple contigs
        options = self.options._replace(group=["non-autosome"])
        self.assertEqual(dropseq.metadata.lookup_contig_groups.run(options), 0)
        with open(self.outputFile, 'r') as f:
            lines = f.readlines()
        actual = set(line.strip() for line in lines)
        expected = set(['chrM', 'chrX', 'chrY'])
        self.assertEqual(actual, expected)
