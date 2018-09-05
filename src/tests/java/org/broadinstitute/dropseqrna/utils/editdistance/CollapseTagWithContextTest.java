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
import htsjdk.samtools.util.CloserUtil;
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

    @Test
    public void testWithCountTags() throws IOException {
        final CollapseTagWithContext clp = new CollapseTagWithContext();
        clp.INPUT = File.createTempFile("CollapseTagWithContextTest.input.", ".sam");;
        clp.INPUT.deleteOnExit();
        clp.COLLAPSE_TAG="XM";
        clp.OUT_TAG="XN";
        clp.CONTEXT_TAGS = Arrays.asList("XC");
        clp.COUNT_TAGS = Arrays.asList("XG");
        clp.OUTPUT = File.createTempFile("CollapseTagWithContextTest.output.", ".sam");
        clp.OUTPUT.deleteOnExit();
        clp.READ_MQ = 0; // Do not filter on mapping quality (these are unmapped reads)
        final SAMRecordSetBuilder builder = createUnmappedFragments(12);
        final ArrayList<SAMRecord> records = new ArrayList<>(builder.getRecords());
        final String collapseTagValue = makeRandomBaseString(8);
        final String contextTagValue = makeRandomBaseString(6);
        final String sharedCountTagValue = makeRandomBaseString(10);

        // 8 records with same collapseTagValue, and 2 unrelated countTagValues
        for (int i = 0; i < 4; ++i) {
            final SAMRecord rec = records.get(i);
            rec.setAttribute(clp.COLLAPSE_TAG, collapseTagValue);
            rec.setAttribute(clp.CONTEXT_TAGS.get(0), contextTagValue);
            rec.setAttribute(clp.COUNT_TAGS.get(0), sharedCountTagValue);
        }
        final String sharedCountTagValue2 = makeRandomBaseString(10);
        for (int i = 4; i < 8; ++i) {
            final SAMRecord rec = records.get(i);
            rec.setAttribute(clp.COLLAPSE_TAG, collapseTagValue);
            rec.setAttribute(clp.CONTEXT_TAGS.get(0), contextTagValue);
            rec.setAttribute(clp.COUNT_TAGS.get(0), sharedCountTagValue2);
        }

        // 4 records with same ed1CollapseTagValue, 2 identical countTagValues, and 2 that are ED 1
        final String ed1CollapseTagValue = alterBaseString(collapseTagValue, 1);
        String countTagValue = makeRandomBaseString(10);
        for (int i = 8; i < 10; ++i) {
            final SAMRecord rec = records.get(i);
            rec.setAttribute(clp.COLLAPSE_TAG, ed1CollapseTagValue);
            rec.setAttribute(clp.CONTEXT_TAGS.get(0), contextTagValue);
            rec.setAttribute(clp.COUNT_TAGS.get(0), countTagValue);
        }
        for (int i = 10; i < records.size(); ++i) {
            final SAMRecord rec = records.get(i);
            rec.setAttribute(clp.COLLAPSE_TAG, ed1CollapseTagValue);
            rec.setAttribute(clp.CONTEXT_TAGS.get(0), contextTagValue);
            rec.setAttribute(clp.COUNT_TAGS.get(0), alterBaseString(countTagValue, 1));
        }
        final SAMFileHeader header = builder.getHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(header, true, clp.INPUT, null);
        for (final SAMRecord rec : records) {
            writer.addAlignment(rec);
        }
        writer.close();
        Assert.assertEquals(clp.doWork(), 0);
        // ed1CollapseTagValue should be selected, because it has 4 distinct values for the COUNT_TAG
        SamReader samReader = SamReaderFactory.makeDefault().open(clp.OUTPUT);
        for (final SAMRecord rec : samReader) {
            Assert.assertEquals(rec.getAttribute(clp.OUT_TAG), ed1CollapseTagValue, rec.getSAMString());
        }
        CloserUtil.close(samReader);

        // Test the same input but without COUNT_TAGS, and confirm different result.
        final CollapseTagWithContext clp2 = new CollapseTagWithContext();
        clp2.INPUT = clp.INPUT;
        clp2.COLLAPSE_TAG = clp.COLLAPSE_TAG;
        clp2.OUT_TAG = clp.OUT_TAG;
        clp2.CONTEXT_TAGS = clp.CONTEXT_TAGS;
        clp2.OUTPUT = clp.OUTPUT;
        clp2.READ_MQ = clp.READ_MQ;
        Assert.assertEquals(clp2.doWork(), 0);
        // collapseTagValue should be selected, because it has 3 reads
        samReader = SamReaderFactory.makeDefault().open(clp.OUTPUT);
        for (final SAMRecord rec : samReader) {
            Assert.assertEquals(rec.getAttribute(clp.OUT_TAG), collapseTagValue, rec.getSAMString());
        }
        CloserUtil.close(samReader);

        // Test the same input, but collapse COUNT_TAG values with ED=1
        final CollapseTagWithContext clp3 = new CollapseTagWithContext();
        clp3.INPUT = clp.INPUT;
        clp3.COLLAPSE_TAG = clp.COLLAPSE_TAG;
        clp3.OUT_TAG = clp.OUT_TAG;
        clp3.CONTEXT_TAGS = clp.CONTEXT_TAGS;
        clp3.COUNT_TAGS = clp.COUNT_TAGS;
        clp3.COUNT_TAGS_EDIT_DISTANCE = 1;
        clp3.OUTPUT = clp.OUTPUT;
        clp3.READ_MQ = clp.READ_MQ;
        clp3.MIN_COUNT = 1;
        Assert.assertEquals(clp3.doWork(), 0);
        // collapseTagValue should be selected, because it has 3 reads
        samReader = SamReaderFactory.makeDefault().open(clp.OUTPUT);
        for (final SAMRecord rec : samReader) {
            Assert.assertEquals(rec.getAttribute(clp.OUT_TAG), collapseTagValue, rec.getSAMString());
        }
        CloserUtil.close(samReader);
    }
}
