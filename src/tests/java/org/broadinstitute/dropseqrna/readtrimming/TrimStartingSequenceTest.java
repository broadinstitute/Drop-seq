/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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
package org.broadinstitute.dropseqrna.readtrimming;

import com.google.common.base.Strings;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class TrimStartingSequenceTest {
    private static final File INPUT = new File("testdata/org/broadinstitute/dropseq/readtrimming/N701.subset.tagged_filtered.sam");
    @Test
    public void testTrimStartingSequenceClp() {
        final TrimStartingSequence clp = new TrimStartingSequence();
        clp.INPUT = INPUT;
        clp.OUTPUT = TestUtils.getTempReportFile("TrimStartingSequenceTest.", ".sam");
        clp.OUTPUT_SUMMARY = TestUtils.getTempReportFile("TrimStartingSequenceTest.", ".summary");
        clp.SEQUENCE = "AAGCAGTGGTATCAACGCAGAGTGAATGGG";
        clp.MISMATCHES=0;
        clp.NUM_BASES=5;
        Assert.assertEquals(clp.doWork(), 0);

        // Confirm that the expected stuff is trimmed.
        final SamReader outputReader = SamReaderFactory.makeDefault().open(clp.OUTPUT);
        final SAMRecordIterator outputIterator = outputReader.iterator();
        final SamReader inputReader = SamReaderFactory.makeDefault().open(clp.INPUT);
        for (final SAMRecord inputRec: inputReader) {
            Assert.assertTrue(outputIterator.hasNext());
            final SAMRecord outputRec = outputIterator.next();
            final String inputBases = inputRec.getReadString();
            final String outputBases = outputRec.getReadString();
            Assert.assertTrue(inputBases.endsWith(outputBases));
            final String trimmedBases = inputBases.substring(0, inputBases.length() - outputBases.length());
            Assert.assertTrue(clp.SEQUENCE.endsWith(trimmedBases));
        }
        Assert.assertFalse(outputIterator.hasNext());
        CloserUtil.close(Arrays.asList(outputReader, inputReader));
    }

    /**
     * @param expectedResult 0: fully trimmed, i.e. bases unchanged but all Q2
     *                       -1: not trimmed
     *                       positive: trimmed this amount
     */
    @Test(dataProvider = "testTrimStartingSequenceDataProvider")
    public void testTrimStartingSequence(
            final String testName,
            final String read,
            final String adapter,
            final int maxMismatches,
            final int minMatchLength,
            final int expectedResult) {
        final TrimStartingSequence clp = new TrimStartingSequence();
        clp.INPUT = TestUtils.getTempReportFile("testTrimStartingSequence.input.", ".sam");
        clp.OUTPUT = TestUtils.getTempReportFile("testTrimStartingSequence.trimmed.", ".sam");
        clp.OUTPUT_SUMMARY = TestUtils.getTempReportFile("testTrimStartingSequence.", ".summary.txt");
        clp.SEQUENCE = adapter;
        clp.MISMATCHES = maxMismatches;
        clp.NUM_BASES = minMatchLength;
        SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(new SAMFileHeader(), false, clp.INPUT);
        final SAMRecord rawRec = new SAMRecord(writer.getFileHeader());
        rawRec.setReadName("trim.test.read");
        rawRec.setReadString(read);
        rawRec.setReadUnmappedFlag(true);
        rawRec.setBaseQualityString(Strings.repeat("H", read.length()));
        writer.addAlignment(rawRec);
        writer.close();
        Assert.assertEquals(clp.doWork(), 0);
        SamReader samReader = SamReaderFactory.makeDefault().open(clp.OUTPUT);
        final SAMRecord trimmedRec = samReader.iterator().next();
        CloserUtil.close(samReader);
        if (expectedResult == 0) {
            // Fully trimmed, so instead of trimming, just set all bases to Q3
            Assert.assertEquals(trimmedRec.getReadLength(), read.length());
            Assert.assertEquals(trimmedRec.getBaseQualityString(), Strings.repeat("$", read.length()));
            Assert.assertNull(trimmedRec.getAttribute(clp.TRIM_TAG));
        } else if (expectedResult == -1) {
            // Not trimmed
            Assert.assertEquals(trimmedRec.getReadLength(), read.length());
            Assert.assertEquals(trimmedRec.getBaseQualityString(), rawRec.getBaseQualityString());
            Assert.assertNull(trimmedRec.getAttribute(clp.TRIM_TAG));
        } else {
            // trimmed
            Assert.assertEquals(trimmedRec.getReadLength(), read.length() - expectedResult);
            Assert.assertEquals(trimmedRec.getIntegerAttribute(clp.TRIM_TAG).intValue(), expectedResult);
        }
    }

    @DataProvider(name = "testTrimStartingSequenceDataProvider")
    public Object[][] testTrimStartingSequenceDataProvider() {
        return new Object[][] {
                // test name, read, adapter, maxMismatches, minMatchLength, expectedResult
                {"basic exact match", "ACGTAAAACCCC", "ACGT", 0, 4, 4},
                {"mismatch no match", "AGGTAAAACCCC", "ACGT", 0, 4, -1},
                {"mismatch allowed", "AGGTAAAACCC", "ACGT", 1, 3, 4},
                {"too short no match", "ACGTAAAACCCC", "AAAACGT", 0, 5, -1},
                {"partial match", "ACGTAAAACCCC", "AAAACGT", 0, 4, 4},
                {"full trim exact match", "AAAAAGGGACGT", "AAAAAGGGACGT", 0, 4, 0},
                {"full trim mismatch allowed", "AAAAAGGGACCT", "AAAAAGGGACGT", 1, 4, 0},
                {"full trim too many mismatches", "AAAAAGGGACCT", "AAAAAGGGACGT", 0, 4, -1},
                {"exact match with N", "ACGTAAAACCCC", "ACNT", 0, 4, 4},
        };
    }
}
