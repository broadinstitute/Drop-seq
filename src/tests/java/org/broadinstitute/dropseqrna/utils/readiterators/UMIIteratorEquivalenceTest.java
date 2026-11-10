/*
 * MIT License
 *
 * Copyright 2026 Broad Institute
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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorStrategy;
import org.broadinstitute.dropseqrna.barnyard.DigitalExpressionTestUtil;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.annotation.LocusFunction;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.LinkedHashMap;

/**
 * When reads are not retained UMIIterator sorts (cell, gene, UMI, count) records instead of SAMRecords.
 * It must give exactly the same UMICollections, in the same order, as the SAMRecord path that
 * retainReads(true) uses.  The reads themselves are not compared.
 */
public class UMIIteratorEquivalenceTest {

    private static final File IN_FILE = DigitalExpressionTestUtil.IN_FILE;
    private static final List<LocusFunction> CODING_UTR = Arrays.asList(LocusFunction.CODING, LocusFunction.UTR);
    private static final List<LocusFunction> WITH_INTRONIC =
            Arrays.asList(LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC);

    @DataProvider(name = "matrix")
    public Object[][] matrix() {
        final List<Object[]> rows = new ArrayList<>();
        for (final int readMQ : new int[]{0, 3, 10})
            for (final boolean cellFirst : new boolean[]{true, false})
                for (final boolean assignAll : new boolean[]{true, false})
                    for (final boolean cellList : new boolean[]{true, false})
                        for (final boolean intervals : new boolean[]{true, false})
                            for (final int numInputs : new int[]{1, 2})
                                for (final boolean record : new boolean[]{true, false})
                                    for (final boolean intronic : new boolean[]{true, false})
                                        rows.add(new Object[]{readMQ, cellFirst, assignAll, cellList, intervals, numInputs, record, intronic});
        return rows.toArray(new Object[0][]);
    }

    private static class Result {
        final List<String> collections = new ArrayList<>();
        Set<String> cellsSeen;
    }

    private static Result run(final boolean retainReads, final int readMQ, final boolean cellFirst, final boolean assignAll,
                              final boolean cellList, final boolean intervals, final int numInputs,
                              final boolean record, final boolean intronic) {
        final UMIIterator.UMIIteratorBuilder builder = new UMIIterator.UMIIteratorBuilder(
                SamFileMergeUtil.mergeInputs(Collections.nCopies(numInputs, IN_FILE), false),
                "gn", "gs", "gf", StrandStrategy.SENSE, intronic ? WITH_INTRONIC : CODING_UTR,
                FunctionalDataProcessorStrategy.STARSOLO, "XC", "XM", readMQ)
                .assignReadsToAllGenes(assignAll).cellFirstSort(cellFirst).recordCellsInInput(record).retainReads(retainReads);
        if (cellList)
            builder.setCellBarcodes(Arrays.asList(Arrays.copyOf(DigitalExpressionTestUtil.barcodes, 3)));
        if (intervals) {
            final SAMFileHeader header = SamReaderFactory.makeDefault().getFileHeader(IN_FILE);
            final IntervalList il = new IntervalList(header);
            il.add(new Interval("HUMAN_3", 42642106, 42650000));
            builder.setIntervals(il);
        }
        final Result result = new Result();
        final UMIIterator it = builder.build();
        UMICollection c;
        while (it.hasNext() && (c = it.next()) != null) {
            // LinkedHashMap: the iteration order of the UMIs must match too, because it shows up in some output files
            final LinkedHashMap<String, Integer> umis = new LinkedHashMap<>();
            for (final String umi : c.getMolecularBarcodes())
                umis.put(umi, c.getMolecularBarcodeCounts().getCountForKey(umi));
            result.collections.add(c.getCellBarcode() + "|" + c.getGeneName() + "|" + umis);
        }
        result.cellsSeen = it.getCellBarcodesSeen();
        it.close();
        return result;
    }

    @Test(dataProvider = "matrix")
    public void testSameAsSamRecordSort(final int readMQ, final boolean cellFirst, final boolean assignAll,
                                        final boolean cellList, final boolean intervals, final int numInputs,
                                        final boolean record, final boolean intronic) {
        final Result expected = run(true, readMQ, cellFirst, assignAll, cellList, intervals, numInputs, record, intronic);
        final Result actual = run(false, readMQ, cellFirst, assignAll, cellList, intervals, numInputs, record, intronic);
        Assert.assertFalse(expected.collections.isEmpty(), "test data should produce UMI collections");
        Assert.assertEquals(actual.collections, expected.collections);
        Assert.assertEquals(actual.cellsSeen, expected.cellsSeen);
    }

    /** Checks the matrix actually exercises the multi-mapper path: low MQ keeps reads that MQ 10 rejects. */
    @Test
    public void testLowMapQualityKeepsMoreReads() {
        final Result low = run(false, 0, true, false, false, false, 1, false, false);
        final Result high = run(false, 10, true, false, false, false, 1, false, false);
        Assert.assertNotEquals(low.collections, high.collections);
    }
}
