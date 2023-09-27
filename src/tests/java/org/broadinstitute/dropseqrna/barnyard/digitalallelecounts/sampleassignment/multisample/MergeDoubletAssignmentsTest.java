/*
 * MIT License
 *
 * Copyright 2020 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample;

import java.io.File;
import java.util.Arrays;

import htsjdk.samtools.util.IOUtil;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.MergeDoubletAssignments;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;


public class MergeDoubletAssignmentsTest {

	private final static String rootDir="testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/sampleassignment/";
	
    private static File DOUBLETS_SPLIT_BAM0 = new File(rootDir+"/MergeDoubletAssignmentsTest.SplitBam.0.txt");
    private static File DOUBLETS_SPLIT_BAM1 = new File(rootDir+"/MergeDoubletAssignmentsTest.SplitBam.1.txt");
    private static File EXPECTED_DOUBLETS_FILE = new File (rootDir+"/MergeDoubletAssignmentsTest.txt");

    @Test
    public void testDoWork() {
        File mergedDoubletsFile = TestUtils.getTempReportFile("MergeDoubletAssignments", ".txt");

        MergeDoubletAssignments merger = new MergeDoubletAssignments();
        merger.INPUT = Arrays.asList(DOUBLETS_SPLIT_BAM0, DOUBLETS_SPLIT_BAM1);
        merger.OUTPUT = mergedDoubletsFile;
        Assert.assertEquals(merger.doWork(), 0);

        IOUtil.assertFilesEqual(mergedDoubletsFile, EXPECTED_DOUBLETS_FILE);
    }
}
