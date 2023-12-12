/*
 * MIT License
 *
 * Copyright 2023 Broad Institute
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
package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.stream.StreamSupport;

public class SamFileMergeUtilTest {
    private static File COORDINATE_SORTED_SAM = new File("testdata/org/broadinstitute/dropseq/barnyard/MarkChimericReads.input2.sam");
    private static File QUERYNAME_SORTED_SAM = new File("testdata/org/broadinstitute/dropseq/readtrimming/N701.subset.tagged_filtered.sam");

    private long countSamRecords(final Iterator<SAMRecord> iterator) {
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(iterator, Spliterator.ORDERED), false).count();
    }

    @Test(dataProvider = "testMergeInputsDataProvider")
    public void testMergeInputs(final SAMFileHeader.SortOrder sortOrder) {
        final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
        final SamReader coordinateSortedSam = samReaderFactory.open(COORDINATE_SORTED_SAM);
        final SAMSequenceDictionary sequenceDictionary = coordinateSortedSam.getFileHeader().getSequenceDictionary();
        final long numCoordinateSortedRecords = countSamRecords(coordinateSortedSam.iterator());
        final SamReader querynameSortedSam = samReaderFactory.open(QUERYNAME_SORTED_SAM);
        final long numQuerynameSortedRecords = countSamRecords(querynameSortedSam.iterator());
        CloserUtil.close(Arrays.asList(coordinateSortedSam, querynameSortedSam));
        final SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(
                Arrays.asList(COORDINATE_SORTED_SAM, QUERYNAME_SORTED_SAM), sortOrder, samReaderFactory, true);
        Assert.assertEquals(headerAndIterator.header.getSequenceDictionary(), sequenceDictionary);
        SAMRecord prev = null;
        long numRecords = 0;
        final SAMRecordComparator comparator = sortOrder.getComparatorInstance();
        while (headerAndIterator.iterator.hasNext()) {
            ++numRecords;
            final SAMRecord thisRec = headerAndIterator.iterator.next();
            if (prev != null) {
                Assert.assertTrue(comparator.compare(prev, thisRec) <= 0);
            }
            prev = thisRec;
        }
        CloserUtil.close(headerAndIterator.iterator);
        Assert.assertEquals(numRecords, numCoordinateSortedRecords + numQuerynameSortedRecords);
    }

    @DataProvider(name = "testMergeInputsDataProvider")
    private Object[][] testBasicDataProvider() {
        return new Object[][]{
                {SAMFileHeader.SortOrder.queryname},
                {SAMFileHeader.SortOrder.coordinate}
        };
    }

}
