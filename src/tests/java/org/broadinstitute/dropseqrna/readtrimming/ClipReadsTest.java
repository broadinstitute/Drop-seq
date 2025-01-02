/*
 * MIT License
 *
 * Copyright 2025 Broad Institute
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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;

public class ClipReadsTest {
    private static final File TESTDATA_DIR = new File("testdata/org/broadinstitute/dropseq/readtrimming");
    private static final File PAIRED_INPUT = new File(TESTDATA_DIR, "paired_end.28_technical.sam");
    private static final File PAIRED_INPUT_SHORT_READ = new File(TESTDATA_DIR, "paired_end.28_technical.short_read.sam");

    @Test
    public void testClipReads() {
        final File tempDir = TestUtils.createTempDirectory("ClipReadsTest.");
        final ClipReads clp = new ClipReads();
        clp.INPUT = PAIRED_INPUT;
        clp.OUTPUT = TestUtils.getTempReportFile("ClipReadsTest.", ".sam");
        clp.OUTPUT.deleteOnExit();
        clp.TMP_DIR = Collections.singletonList(tempDir);
        clp.BASE_RANGE = "1-16:17-28";
        clp.WHICH_READ = Collections.singletonList(AbstractTrimmerClp.FIRST_OF_PAIR);
        Assert.assertEquals(clp.doWork(), 0);

        final SamReader inputReader = SamReaderFactory.makeDefault().open(clp.INPUT);
        final SamReader actualReader = SamReaderFactory.makeDefault().open(clp.OUTPUT);
        final SAMRecordIterator inputIterator = inputReader.iterator();
        final SAMRecordIterator actualIterator = actualReader.iterator();
        while (inputIterator.hasNext() && actualIterator.hasNext()) {
            final SAMRecord inputRecord = inputIterator.next();
            final SAMRecord actualRecord = actualIterator.next();
            Assert.assertEquals(actualRecord.getReadName(), inputRecord.getReadName());
            Assert.assertEquals(actualRecord.getFirstOfPairFlag(), inputRecord.getFirstOfPairFlag(), actualRecord.getReadName());
            if (actualRecord.getFirstOfPairFlag()) {
                final String inputRead = inputRecord.getReadString().substring(28);
                final String actualRead = actualRecord.getReadString();
                Assert.assertEquals(actualRead, inputRead, actualRecord.getReadName());
            } else {
                Assert.assertEquals(actualRecord, inputRecord);
            }
        }
        Assert.assertFalse(inputIterator.hasNext());
        Assert.assertFalse(actualIterator.hasNext());
        CloserUtil.close(inputReader);
        CloserUtil.close(actualReader);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testShortRead() {
        final File tempDir = TestUtils.createTempDirectory("ClipReadsTest.");
        final ClipReads clp = new ClipReads();
        clp.INPUT = PAIRED_INPUT_SHORT_READ;
        clp.OUTPUT = TestUtils.getTempReportFile("ClipReadsTest.", ".sam");
        clp.OUTPUT.deleteOnExit();
        clp.TMP_DIR = Collections.singletonList(tempDir);
        clp.BASE_RANGE = "1-16:17-28";
        clp.WHICH_READ = Collections.singletonList(AbstractTrimmerClp.FIRST_OF_PAIR);
        clp.doWork();

    }
}
