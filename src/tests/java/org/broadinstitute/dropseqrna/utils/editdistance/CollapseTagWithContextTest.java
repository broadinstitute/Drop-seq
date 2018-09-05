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
package org.broadinstitute.dropseqrna.utils.editdistance;

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class CollapseTagWithContextTest {

    private static final File TEST_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.bam");

    private Random random = new Random(0);

    /**
     * TODO: Very lame test -- just confirms that program doesn't crash and return 0 exit status.  I don't think
     * any collapsing happens on this input.
     */
    @Test
    public void testTrivial() throws IOException {
        final CollapseTagWithContext clp = new CollapseTagWithContext();
        clp.INPUT = TEST_FILE;
        clp.COLLAPSE_TAG="XM";
        clp.OUT_TAG="XN";
        clp.CONTEXT_TAGS = Arrays.asList("XC", "XG");
        clp.OUTPUT = File.createTempFile("CollapseTagWithContextTest.", ".bam");
        clp.OUTPUT.deleteOnExit();
        clp.ADAPTIVE_ED_METRICS_FILE = File.createTempFile("CollapseTagWithContextTest.", ".adaptive_ed_metrics");
        clp.ADAPTIVE_ED_METRICS_FILE.deleteOnExit();
        Assert.assertEquals(clp.doWork(), 0);
    }

    private final SAMRecordSetBuilder createUnmappedFragments(int numRecords) {
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder();
        for (int i = 0; i < numRecords; ++i) {
            builder.addUnmappedFragment("read" + (i+1));
        }
        return builder;
    }

    private final String makeRandomBaseString(final int length) {
        final byte[] bases = new byte[length];
        for(int i = 0; i < length; ++i) {
            bases[i] = getRandomBase();
        }
        return StringUtil.bytesToString(bases);
    }

    private byte getRandomBase() {
        return SequenceUtil.VALID_BASES_UPPER[random.nextInt(SequenceUtil.VALID_BASES_UPPER.length)];
    }

    private final byte alterBase(final byte base) {
        byte ret = getRandomBase();
        while (ret == base) {
            ret = getRandomBase();
        }
        return ret;
    }

    private final String alterBaseString(final String baseString, final int numChanges) {
        final byte[] bases = StringUtil.stringToBytes(baseString);
        if (numChanges > baseString.length()) {
            throw new IllegalArgumentException("Too many changes requested");
        }
        final Set<Integer> mutatedPositions = new HashSet<>();
        int changesSoFar = 0;
        while (changesSoFar < numChanges) {
            int positionToChange = random.nextInt(bases.length);
            while (mutatedPositions.contains(positionToChange)) {
                positionToChange = random.nextInt(bases.length);
            }
            mutatedPositions.add(positionToChange);
            bases[positionToChange] = alterBase(bases[positionToChange]);
            ++changesSoFar;
        }
        return StringUtil.bytesToString(bases);
    }



    @Test
    public void testNoCountTags() throws IOException, CloneNotSupportedException {
        final CollapseTagWithContext clp = new CollapseTagWithContext();
        clp.INPUT = File.createTempFile("CollapseTagWithContextTest.input.", ".sam");;
        clp.INPUT.deleteOnExit();
        clp.COLLAPSE_TAG="XM";
        clp.OUT_TAG="XN";
        clp.CONTEXT_TAGS = Arrays.asList("XC", "XG");
        clp.OUTPUT = File.createTempFile("CollapseTagWithContextTest.output.", ".sam");
        clp.OUTPUT.deleteOnExit();
        clp.READ_MQ = 0; // Do not filter on mapping quality (these are unmapped reads)
        final SAMRecordSetBuilder builder = createUnmappedFragments(4);
        final ArrayList<SAMRecord> records = new ArrayList<>(builder.getRecords());
        final SAMRecord r1 = records.get(0);
        final SAMRecord r2 = records.get(1);
        final SAMRecord r3 = records.get(2);
        final SAMRecord r4 = records.get(3);
        final String tagValue = makeRandomBaseString(8);
        final String[] contextStrings = new String[clp.CONTEXT_TAGS.size()];
        for (int i = 0; i < contextStrings.length; ++i) {
            contextStrings[i] = makeRandomBaseString(6);
            // These 3 reads have same context
            r1.setAttribute(clp.CONTEXT_TAGS.get(i), contextStrings[i]);
            r2.setAttribute(clp.CONTEXT_TAGS.get(i), contextStrings[i]);
            r3.setAttribute(clp.CONTEXT_TAGS.get(i), contextStrings[i]);
        }
        // This read has different context
        r4.setAttribute(clp.CONTEXT_TAGS.get(0), contextStrings[0]);
        r4.setAttribute(clp.CONTEXT_TAGS.get(1), alterBaseString(contextStrings[1], 1));

        // read 2 has ED close enough, and same context
        // read 3 has ED too far, and same context
        // read 4 has ED close enough, but different context
        r1.setAttribute(clp.COLLAPSE_TAG, tagValue);
        r2.setAttribute(clp.COLLAPSE_TAG, alterBaseString(tagValue, clp.EDIT_DISTANCE));
        r3.setAttribute(clp.COLLAPSE_TAG, alterBaseString(tagValue, clp.EDIT_DISTANCE + 1));
        r4.setAttribute(clp.COLLAPSE_TAG, alterBaseString(tagValue, clp.EDIT_DISTANCE));

        final SAMFileHeader header = builder.getHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(header, true, clp.INPUT, null);
        for (final SAMRecord rec : records) {
            writer.addAlignment(rec);
        }
        // Add one more read with same value for collapse tag as read 1, so that value will be the one that wins.
        final SAMRecord r5 = (SAMRecord)r1.clone();
        r5.setReadName("read5");
        writer.addAlignment(r5);
        writer.close();

        Assert.assertEquals(clp.doWork(), 0);

        final Map<String, SAMRecord> collapsedRecords = new HashMap<>();
        final SamReader samReader = SamReaderFactory.makeDefault().open(clp.OUTPUT);
        for (final SAMRecord rec : samReader) {
            collapsedRecords.put(rec.getReadName(), rec);
        }
        final SAMRecord r1Collapsed = collapsedRecords.get(r1.getReadName());
        Assert.assertEquals(r1Collapsed.getAttribute(clp.COLLAPSE_TAG), r1.getAttribute(clp.COLLAPSE_TAG));
        Assert.assertEquals(r1Collapsed.getAttribute(clp.OUT_TAG), r1Collapsed.getAttribute(clp.COLLAPSE_TAG));
        final SAMRecord r2Collapsed = collapsedRecords.get(r2.getReadName());
        Assert.assertEquals(r2Collapsed.getAttribute(clp.OUT_TAG), r1Collapsed.getAttribute(clp.COLLAPSE_TAG));
        final SAMRecord r3Collapsed = collapsedRecords.get(r3.getReadName());
        Assert.assertEquals(r3Collapsed.getAttribute(clp.OUT_TAG), r3.getAttribute(clp.COLLAPSE_TAG));
        final SAMRecord r4Collapsed = collapsedRecords.get(r4.getReadName());
        Assert.assertEquals(r4Collapsed.getAttribute(clp.OUT_TAG), r4.getAttribute(clp.COLLAPSE_TAG));
    }
}
