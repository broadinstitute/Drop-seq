/*
 * MIT License
 *
 * Copyright 2019 Broad Institute
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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class SequenceDictionaryIntersectionTest {

    private static final File TESTDATA_DIR = new File("testdata/org/broadinstitute/dropseq/utils/SequenceDictionaryIntersectionTest");
    private static final String CHR_BASENAME = "chr.";
    private static final String NO_CHR_BASENAME = "no_chr.";
    private static final int NUM_SEQUENCES = 85;
    private static final int NUM_NON_CHR_SEQUENCES = 60;

    @Test(dataProvider = "ssdiDataProvider")
    public void testWithDescriptions(
            final Object o1, final Object o2,
            final String description1, final String description2,
            final int expectedIntersectionSize) {
        final SequenceDictionaryIntersection sdi =
                new SequenceDictionaryIntersection(o1, description1, o2, description2);
        Assert.assertEquals(sdi.getIntersection().size(), expectedIntersectionSize);
        System.out.println("Terse: " + sdi.message(false));
        System.out.println("Verbose: " + sdi.message(true));
    }

    @Test(dataProvider = "ssdiDataProvider")
    public void testWithoutDescriptions(
            final Object o1, final Object o2,
            final String description1, final String description2,
            final int expectedIntersectionSize) {
        final SequenceDictionaryIntersection sdi =
                new SequenceDictionaryIntersection(o1, o2);
        Assert.assertEquals(sdi.getIntersection().size(), expectedIntersectionSize);
        System.out.println("Terse: " + sdi.message(false));
        System.out.println("Verbose: " + sdi.message(true));
    }

    @DataProvider(name="ssdiDataProvider")
    public Object[][] ssdiDataProvider() {
        final ArrayList<Object[]> ret = new ArrayList<>();
        final File leftDictFile = new File(TESTDATA_DIR, CHR_BASENAME + "sam");
        final SamReader leftDict = SamReaderFactory.makeDefault().open(leftDictFile);

        for (final boolean chr: new boolean[]{true, false}) {
            final String basename;
            final int expectedIntersection;
            if (chr) {
                basename = CHR_BASENAME;
                expectedIntersection = NUM_SEQUENCES;
            } else {
                basename = NO_CHR_BASENAME;
                expectedIntersection = NUM_NON_CHR_SEQUENCES;
            }
            final File samFile = new File(TESTDATA_DIR, basename + "sam");
            final SamReader samReader = SamReaderFactory.makeDefault().open(samFile);
            ret.add(new Object[]{leftDict, samReader, leftDictFile.getName(), samFile.getName(), expectedIntersection});
            ret.add(new Object[]{leftDict, samReader.getFileHeader(), leftDictFile.getName(), samFile.getName(), expectedIntersection});
            ret.add(new Object[]{leftDict, samReader.getFileHeader().getSequenceDictionary(), leftDictFile.getName(), samFile.getName(), expectedIntersection});
            final File vcfFile = new File(TESTDATA_DIR, basename + "vcf");
            final VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);
            ret.add(new Object[]{leftDict, vcfReader, leftDictFile.getName(), vcfFile.getName(), expectedIntersection});
            final File intervalListFile = new File(TESTDATA_DIR, basename + "interval_list");
            final IntervalList intervalList = IntervalList.fromFile(intervalListFile);
            ret.add(new Object[]{leftDict, intervalList, leftDictFile.getName(), intervalListFile.getName(), expectedIntersection});
        }

        return ret.toArray(new Object[ret.size()][]);
    }

    @Test
    public void testNoMatch() {
        final File leftDictFile = new File(TESTDATA_DIR, CHR_BASENAME + "sam");
        final SamReader leftDict = SamReaderFactory.makeDefault().open(leftDictFile);
        try {
            final SAMSequenceDictionary nonMatchingDict = new SAMSequenceDictionary(
                    leftDict.getFileHeader().getSequenceDictionary().getSequences().stream().
                            map((s) -> new SAMSequenceRecord("xyz_" + s.getSequenceName(), s.getSequenceLength())).
                            collect(Collectors.toList()));
            SequenceDictionaryIntersection sdi = new SequenceDictionaryIntersection(
                    leftDict, leftDictFile.getName(), nonMatchingDict, "with prefixes");
            Assert.assertEquals(sdi.getIntersection().size(), 0);
            System.out.println("Terse with descriptions: " + sdi.message(false));
            System.out.println("Verbose with descriptions: " + sdi.message(true));

            sdi = new SequenceDictionaryIntersection(leftDict, nonMatchingDict);
            Assert.assertEquals(sdi.getIntersection().size(), 0);
            System.out.println("Terse without descriptions: " + sdi.message(false) + "\n");
            System.out.println("Verbose without descriptions: " + sdi.message(true) + "\n");
        } finally {
            CloserUtil.close(leftDict);
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBadObject() {
        final File leftDictFile = new File(TESTDATA_DIR, CHR_BASENAME + "sam");
        final SamReader leftDict = SamReaderFactory.makeDefault().open(leftDictFile);
        try {
            SequenceDictionaryIntersection sdi = new SequenceDictionaryIntersection(
                    leftDict, leftDictFile.getName(), "Just a string", "not a sequence dictionary");

        } finally {
            CloserUtil.close(leftDict);;
        }
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullObject() {
        final File leftDictFile = new File(TESTDATA_DIR, CHR_BASENAME + "sam");
        final SamReader leftDict = SamReaderFactory.makeDefault().open(leftDictFile);
        try {
            SequenceDictionaryIntersection sdi = new SequenceDictionaryIntersection(
                    leftDict, leftDictFile.getName(),null, "null");

        } finally {
            CloserUtil.close(leftDict);;
        }
    }
}
