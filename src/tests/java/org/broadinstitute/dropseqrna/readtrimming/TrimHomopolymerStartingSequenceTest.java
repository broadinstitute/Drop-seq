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
import java.util.Arrays;
import java.util.Collections;

public class TrimHomopolymerStartingSequenceTest {
    private static final File TESTDATA_DIR = new File("testdata/org/broadinstitute/dropseq/readtrimming");
    private static final File INPUT = new File(TESTDATA_DIR, "prePolyTTrim.paired.sam");

    @Test
    public void testBasic() {
        final File tempDir = TestUtils.createTempDirectory("TrimHomopolymerStartingSequenceTest.");
        final TrimHomopolymerStartingSequence clp = new TrimHomopolymerStartingSequence();
        clp.INPUT = INPUT;
        clp.OUTPUT = TestUtils.getTempReportFile("TrimHomopolymerStartingSequenceTest.", ".sam");
        clp.WHICH_READ = Collections.singletonList(AbstractTrimmerClp.FIRST_OF_PAIR);
        clp.TMP_DIR = Collections.singletonList(tempDir);
        Assert.assertEquals(clp.doWork(), 0);

        final SamReader inputReader = SamReaderFactory.makeDefault().open(clp.INPUT);
        final SamReader actualReader = SamReaderFactory.makeDefault().open(clp.OUTPUT);
        final SAMRecordIterator inputIterator = inputReader.iterator();
        final SAMRecordIterator actualIterator = actualReader.iterator();
        while (inputIterator.hasNext() && actualIterator.hasNext()) {
            final SAMRecord inputRecord = inputIterator.next();
            final SAMRecord actualRecord = actualIterator.next();
            final String readName = actualRecord.getReadName();
            Assert.assertEquals(readName, inputRecord.getReadName());
            Assert.assertEquals(actualRecord.getFirstOfPairFlag(), inputRecord.getFirstOfPairFlag(), readName);
            if (actualRecord.getFirstOfPairFlag()) {
                final String[] readNameFields = readName.split(":");
                if (readNameFields[0].equals("trimmed")) {
                    final int trimLength = Integer.parseInt(readNameFields[1]);
                    final String inputRead = inputRecord.getReadString().substring(trimLength);
                    Assert.assertEquals(actualRecord.getReadString(), inputRead, readName);
                    final String inputQual = inputRecord.getBaseQualityString().substring(trimLength);
                    Assert.assertEquals(actualRecord.getBaseQualityString(), inputQual, readName);
                } else if (readNameFields[0].equals("notrim")) {
                    Assert.assertEquals(actualRecord, inputRecord);
                } else if (readNameFields[0].equals("fulltrim")) {
                    Assert.assertEquals(actualRecord.getReadString(), inputRecord.getReadString(), readName);
                    byte[] expectedQuals = new byte[inputRecord.getReadLength()];
                    Arrays.fill(expectedQuals, TrimHomopolymerStartingSequence.FULL_TRIM_QUALITY_SCORE);
                    Assert.assertEquals(actualRecord.getBaseQualities(), expectedQuals, readName);
                } else {
                    Assert.fail("Unexpected read name: " + readName);
                }
            } else {
                Assert.assertEquals(actualRecord, inputRecord);
            }
        }
        Assert.assertFalse(inputIterator.hasNext());
        Assert.assertFalse(actualIterator.hasNext());
        CloserUtil.close(inputReader);
        CloserUtil.close(actualReader);

    }
}
