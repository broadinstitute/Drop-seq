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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class TrimStartingSequenceTest {
    private static final File INPUT = new File("testdata/org/broadinstitute/dropseq/readtrimming/N701.subset.tagged_filtered.sam");
    @Test
    public void testTrimStartingSequenceClp() throws IOException {
        final TrimStartingSequence clp = new TrimStartingSequence();
        clp.INPUT = INPUT;
        clp.OUTPUT = File.createTempFile("TrimStartingSequenceTest.", ".sam");
        clp.OUTPUT.deleteOnExit();
        clp.OUTPUT_SUMMARY = File.createTempFile("TrimStartingSequenceTest.", ".summary");
        clp.OUTPUT_SUMMARY.deleteOnExit();
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
}
