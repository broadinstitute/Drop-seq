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
    public static final String LENGTH_TAG = "Zl";

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
     * @param expectedResult -1: fully trimmed, i.e. bases unchanged but all Q2
     *                       0: not trimmed
     *                       positive: trimmed this amount
     */
    @Test(dataProvider = "testFixedMismatchesDataProvider")
    public void testFixedMismatches(
            final String testName,
            final String read,
            final String adapter,
            final int maxMismatches,
            final int minMatchLength,
            final int expectedResult) {
        final SAMFileHeader samFileHeader = new SAMFileHeader();
        final SAMRecord rawRec = makeUntrimmedRecord(read, samFileHeader);
        final SAMRecord trimmedRec = trimRead(rawRec, samFileHeader, adapter, maxMismatches, minMatchLength,
                null, false);
        if (expectedResult == -1) {
            assertFullyTrimmed(read, trimmedRec, false);
        } else if (expectedResult == 0) {
            assertNotTrimmed(read, rawRec, trimmedRec);
        } else {
            assertTrimmed(read, trimmedRec, expectedResult, expectedResult, false);
        }
    }

    private void assertTrimmed(String read, SAMRecord trimmedRec, int trimEndPosition, int trimLength, boolean legacy) {
        Assert.assertEquals(trimmedRec.getReadLength(), read.length() - trimEndPosition);
        Assert.assertEquals(trimmedRec.getIntegerAttribute(TrimStartingSequence.DEFAULT_TRIM_TAG).intValue(), trimEndPosition);
        if (!legacy) {
            Assert.assertEquals(trimmedRec.getIntegerAttribute(LENGTH_TAG).intValue(), trimLength);
        }
    }

    private void assertNotTrimmed(String read, SAMRecord rawRec, SAMRecord trimmedRec) {
        Assert.assertEquals(trimmedRec.getReadLength(), read.length());
        Assert.assertEquals(trimmedRec.getBaseQualityString(), rawRec.getBaseQualityString());
        Assert.assertNull(trimmedRec.getAttribute(TrimStartingSequence.DEFAULT_TRIM_TAG));
        Assert.assertNull(trimmedRec.getAttribute(LENGTH_TAG));
    }

    private void assertFullyTrimmed(String read, SAMRecord trimmedRec, boolean legacy) {
        // Fully trimmed, so instead of trimming, just set all bases to Q3
        Assert.assertEquals(trimmedRec.getReadLength(), read.length());
        Assert.assertEquals(trimmedRec.getBaseQualityString(), Strings.repeat("$", read.length()));
        Assert.assertNull(trimmedRec.getAttribute(TrimStartingSequence.DEFAULT_TRIM_TAG));
        if (!legacy) {
            Assert.assertNotNull(trimmedRec.getAttribute(LENGTH_TAG));
        }
    }

    private SAMRecord trimRead(final SAMRecord rawRec, final SAMFileHeader samFileHeader,
                               final String adapter, final Integer maxMismatches, final int minMatchLength,
                               final Double mismatchRate, final boolean legacy) {
        final TrimStartingSequence clp = new TrimStartingSequence();
        clp.INPUT = TestUtils.getTempReportFile("testTrimStartingSequence.input.", ".sam");
        clp.OUTPUT = TestUtils.getTempReportFile("testTrimStartingSequence.trimmed.", ".sam");
        clp.OUTPUT_SUMMARY = TestUtils.getTempReportFile("testTrimStartingSequence.", ".summary.txt");
        clp.SEQUENCE = adapter;
        clp.MISMATCHES = maxMismatches;
        clp.NUM_BASES = minMatchLength;
        clp.LENGTH_TAG = LENGTH_TAG;
        clp.MISMATCH_RATE = mismatchRate;
        clp.LEGACY = legacy;
        SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(samFileHeader, false, clp.INPUT);
        writer.addAlignment(rawRec);
        writer.close();
        Assert.assertEquals(clp.doWork(), 0);
        SamReader samReader = SamReaderFactory.makeDefault().open(clp.OUTPUT);
        final SAMRecord trimmedRec = samReader.iterator().next();
        CloserUtil.close(samReader);
        return trimmedRec;
    }

    private SAMRecord makeUntrimmedRecord(String read, SAMFileHeader header) {
        final SAMRecord rawRec = new SAMRecord(header);
        rawRec.setReadName("trim.test.read");
        rawRec.setReadString(read);
        rawRec.setReadUnmappedFlag(true);
        rawRec.setBaseQualityString(Strings.repeat("H", read.length()));
        return rawRec;
    }

    @DataProvider(name = "testFixedMismatchesDataProvider")
    public Object[][] testFixedMismatchesDataProvider() {
        return new Object[][] {
                // test name, read, adapter, maxMismatches, minMatchLength, expectedResult
                {"basic exact match", "ACGTAAAACCCC", "ACGT", 0, 4, 4},
                {"mismatch no match", "AGGTAAAACCCC", "ACGT", 0, 4, 0},
                {"mismatch allowed", "AGGTAAAACCC", "ACGT", 1, 3, 4},
                {"too short no match", "ACGTAAAACCCC", "AAAACGT", 0, 5, 0},
                {"partial match", "ACGTAAAACCCC", "AAAACGT", 0, 4, 4},
                {"full trim exact match", "AAAAAGGGACGT", "AAAAAGGGACGT", 0, 4, -1},
                {"full trim mismatch allowed", "AAAAAGGGACCT", "AAAAAGGGACGT", 1, 4, -1},
                {"full trim too many mismatches", "AAAAAGGGACCT", "AAAAAGGGACGT", 0, 4, 0},
                {"exact match with N", "ACGTAAAACCCC", "ACNT", 0, 4, 4},
                {"TSO longer than read", "AAGCAGTGGTATCAACGCAGAGTACATGGGGC",
                        "AAGCAGTGGTATCAACGCAGAGTACATGGGAAGCAGTGGTATCAACGCAGAGTACATGGGAAGCAGTGGTATCAACGCAGAGTACATGGG",
                        3, 4, 30}
        };
    }

    @Test(dataProvider = "testMismatchRateDataProvider")
    public void testMismatchRate(
            final String testName,
            final String read,
            final String adapter,
            final double mismatchRate,
            final int minMatchLength,
            final int trimEndPosition,
            final int trimLength) {
        final SAMFileHeader samFileHeader = new SAMFileHeader();
        final SAMRecord rawRec = makeUntrimmedRecord(read, samFileHeader);
        final SAMRecord trimmedRec = trimRead(rawRec, samFileHeader, adapter, null, minMatchLength,
                mismatchRate, false);
        if (trimEndPosition == -1) {
            assertFullyTrimmed(read, trimmedRec, false);
        } else if (trimEndPosition == 0) {
            assertNotTrimmed(read, rawRec, trimmedRec);
        } else {
            assertTrimmed(read, trimmedRec, trimEndPosition, trimLength, false);
        }
    }
    @DataProvider(name = "testMismatchRateDataProvider")
    public Object[][] testMismatchRateDataProvider() {
        return new Object[][]{
                // test name, read sequence, adapter sequence, mismatch rate, min match length,
                // expected trim end position, expected trim length
                {"multiple TSOs",
                        "AAGCAGTGGTATCAACGCAGAGTACATGGGAAGCAGTGGTATCAACGCAGAGTACATGGGAAGCAGTGGTATCAACGCAGAGTACATGGG",
                        "AAGCAGTGGTATCAACGCAGAGTACATGGG", 0, 4, -1, -1},
                {"multiple TSOs mismatches allowed",
                        "AAGCAGTGGTATCAACGCAGAGTACATGGGAAGCAGTGGTATCAACGCAGAGTACATGGGAAGCAGTGGTATCAACGCAGAGTACACGGG",
                        "AAGCAGTGGTATCAACGCAGAGTACATGGG", 0.1, 4, -1, -1},
                {"multiple TSOs too many mismatches",
                        "AAGCAGTGGTATCAACGCAGAGTACACGGGAAGCAGTGGTATCAACGCAGAGTACACGGGAAGCAGTGGTATCAACGCAGAGTACACGGG",
                        "AAGCAGTGGTATCAACGCAGAGTACATGGG", 0.01, 4, 0, 0},
                {"1.5 TSOs",
                        "CGCAGAGTACATGGGAAGCAGTGGTATCAACGCAGAGTACATGGGACGTACGTACGTACGTACGT",
                        "AAGCAGTGGTATCAACGCAGAGTACATGGG", 0, 4, 45, 30},
                {"1.5 TSOs mismatches allowed",
                        "CGCAGAGTACATGGGAAGCAGTGGTATCAACGCAGAGTACACGGGACGTACGTACGTACGTACGT",
                        "AAGCAGTGGTATCAACGCAGAGTACATGGG", 0.1, 4, 45, 30},
                {"1.5 TSOs too many mismatches",
                        "CGCAGAGTACACGGGAAGCAGTGGTATCAACGCAGAGTACACGGGACGTACGTACGTACGTACGT",
                        "AAGCAGTGGTATCAACGCAGAGTACATGGG", 0.01, 4, 0, 0},
                {"partial TSOs mismatches allowed",
                        "GAGTACACGGGACGTACGTACGTACGTACGT",
                        "AAGCAGTGGTATCAACGCAGAGTACATGGG", 0.1, 4, 11, 11},
                {"partial TSOs mismatches too many mismatches",
                        "GAGTACCCGGGACGTACGTACGTACGTACGT",
                        "AAGCAGTGGTATCAACGCAGAGTACATGGG", 0.1, 4, 0, 0},
                {"basic exact match", "ACGTAAAACCCC", "ACGT", 0.0, 4, 4, 4},
                {"mismatch no match", "AGGTAAAACCCC", "ACGT", 0.1, 4, 0, 0},
                {"mismatch allowed", "AGGTAAAACCC", "ACGT", 0.25, 4, 4, 4},
                {"too short no match", "ACGTAAAACCCC", "AAAACGT", 0, 5, 0, 0},
                {"partial match", "ACGTAAAACCCC", "AAAACGT", 0, 4, 4, 4},
                {"full trim exact match", "AAAAAGGGACGT", "AAAAAGGGACGT", 0, 4, -1, -1},
                {"full trim mismatch allowed", "AAAAAGGGACCT", "AAAAAGGGACGT", 0.1, 4, -1, -1},
                {"full trim too many mismatches", "AAAAAGGGACCT", "AAAAAGGGACGT", 0.01, 4, 0, 0},
                {"exact match with N", "ACGTAAAACCCC", "ACNT", 0, 4, 4, 4},
                {"TSO longer than read", "GCAGTGGTATCAACGCAGAGTACATGGGGC",
                        "AAGCAGTGGTATCAACGCAGAGTACATGGG",
                        0, 4, 28, 28}

        };
    }
    /**
     * @param expectedResult -1: fully trimmed, i.e. bases unchanged but all Q2
     *                       0: not trimmed
     *                       positive: trimmed this amount
     */
    @Test(dataProvider = "testLegacyDataProvider")
    public void testLegacy(
            final String testName,
            final String read,
            final String adapter,
            final int maxMismatches,
            final int minMatchLength,
            final int expectedResult) {
        final SAMFileHeader samFileHeader = new SAMFileHeader();
        final SAMRecord rawRec = makeUntrimmedRecord(read, samFileHeader);
        final SAMRecord trimmedRec = trimRead(rawRec, samFileHeader, adapter, maxMismatches, minMatchLength,
                null, true);
        if (expectedResult == -1) {
            assertFullyTrimmed(read, trimmedRec, true);
        } else if (expectedResult == 0) {
            assertNotTrimmed(read, rawRec, trimmedRec);
        } else {
            assertTrimmed(read, trimmedRec, expectedResult, expectedResult, true);
        }
    }

    @DataProvider(name = "testLegacyDataProvider")
    public Object[][] testLegacyDataProvider() {
        return new Object[][] {
                // test name, read, adapter, maxMismatches, minMatchLength, expectedResult
                {"basic exact match", "ACGTAAAACCCC", "ACGT", 0, 4, 4},
                {"mismatch no match", "AGGTAAAACCCC", "ACGT", 0, 4, 0},
                {"mismatch allowed", "AGGTAAAACCC", "ACGT", 1, 3, 4},
                {"too short no match", "ACGTAAAACCCC", "AAAACGT", 0, 5, 0},
                {"partial match", "ACGTAAAACCCC", "AAAACGT", 0, 4, 4},
                {"full trim exact match", "AAAAAGGGACGT", "AAAAAGGGACGT", 0, 4, -1},
                {"full trim mismatch allowed", "AAAAAGGGACCT", "AAAAAGGGACGT", 1, 4, -1},
                {"full trim too many mismatches", "AAAAAGGGACCT", "AAAAAGGGACGT", 0, 4, 0},
                {"exact match with N", "ACGTAAAACCCC", "ACNT", 0, 4, 4},
                // This is incorrect behavior but it is what the old trimmer does
                {"TSO longer than read", "AAGCAGTGGTATCAACGCAGAGTACATGGGGC",
                        "AAGCAGTGGTATCAACGCAGAGTACATGGGAAGCAGTGGTATCAACGCAGAGTACATGGGAAGCAGTGGTATCAACGCAGAGTACATGGG",
                        3, 4, -1}
        };
    }
}
