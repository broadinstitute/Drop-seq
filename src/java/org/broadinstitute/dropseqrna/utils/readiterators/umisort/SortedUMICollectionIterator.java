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
package org.broadinstitute.dropseqrna.utils.readiterators.umisort;

import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.SortingIteratorFactory;
import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.utils.TransformingIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionProcessor;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;
import org.broadinstitute.dropseqrna.utils.sorting.CombiningSortingCollection;

import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * The fast path of {@link UMIIterator} for when reads do not need to be retained in the {@link UMICollection}s.
 * <p>
 * Reads that have passed the SAMRecord filters are reduced to {@link UMIRead}s, and what is sorted and spilled to disk is
 * (cell, gene, UMI, count) rather than whole SAMRecords.  Duplicates are combined as the sort proceeds, see
 * {@link CombiningSortingCollection}.  When low map quality reads are kept (readMQ &lt;= 3), reads are first
 * sorted by name so that reads that map to more than one place can be resolved to at most one gene.
 * <p>
 * The UMICollections, and the order of the UMIs inside them, are the same as when SAMRecords are sorted.
 */
public class SortedUMICollectionIterator implements CloseableIterator<UMICollection> {

    private static final Log log = Log.getInstance(SortedUMICollectionIterator.class);

    /** Number of reads held in memory before a chunk of UMI keys is combined and spilled. */
    private static final int UMI_KEY_MAX_RECORDS_IN_RAM = 2_000_000;

    private final GroupingIterator<UMIKeyCount> groups;
    private final StringInterner geneCache = new StringInterner();
    private final StringInterner umiCache = new StringInterner();
    private final String cellBarcodeTag;
    private final String geneTag;
    private final String molecularBarcodeTag;
    private final ProgressLogger progress;

    /**
     * @param reads the reads after the SAMRecord filters, each tagged with a single gene
     * @param geneFunctionProcessor resolves reads that map to more than one place, when readMQ &lt;= 3
     * @param readMQ if &lt;= 3, non primary reads are in the input and are resolved by read name
     * @param cellFirstSort if true, groups are ordered by cell then gene, otherwise by gene then cell
     * @param progress logs the progress of the reads as they are reduced
     */
    public SortedUMICollectionIterator(final Iterator<SAMRecord> reads, final GeneFunctionProcessor geneFunctionProcessor,
                                       final String cellBarcodeTag, final String geneTag, final String geneStrandTag,
                                       final String geneFunctionTag, final String molecularBarcodeTag,
                                       final int readMQ, final boolean cellFirstSort, final ProgressLogger progress) {
        this.cellBarcodeTag = cellBarcodeTag;
        this.geneTag = geneTag;
        this.molecularBarcodeTag = molecularBarcodeTag;
        this.progress = progress;

        // The read name, strand and function are only needed to resolve reads that map to more than one place.
        final boolean multiMappersKept = readMQ <= 3;
        Iterator<UMIRead> umiReads = new SAMToUMIReadIterator(reads, multiMappersKept, geneStrandTag, geneFunctionTag);
        if (multiMappersKept)
            umiReads = resolveMultiMappers(umiReads, geneFunctionProcessor);

        final Iterator<UMIKeyCount> keys = new TransformingIterator<UMIRead, UMIKeyCount>(umiReads) {
            private long position = 0;

            @Override
            public UMIKeyCount next() {
                final UMIRead r = this.underlyingIterator.next();
                return new UMIKeyCount(r.cell, r.gene, r.umi, 1, position++);
            }
        };
        final CloseableIterator<UMIKeyCount> sorted = SortingIteratorFactory.createCombining(keys,
                UMIKeyCount.sortOrder(cellFirstSort), new UMIKeyCountCodec(), UMIKeyCount::combine,
                UMI_KEY_MAX_RECORDS_IN_RAM, null);
        log.info("Sorting finished.");
        this.groups = new GroupingIterator<>(sorted, UMIKeyCount.groupOrder(cellFirstSort));
    }

    /**
     * Slim version of the low map quality handling of UMIIterator: sort by read name, resolve the reads of each name to
     * at most one read, and return those for further sorting.
     */
    private static Iterator<UMIRead> resolveMultiMappers(final Iterator<UMIRead> reads, final GeneFunctionProcessor gfp) {
        final Comparator<UMIRead> byName = Comparator.comparing(r -> r.readName);
        final CloseableIterator<UMIRead> nameSorted = SortingIteratorFactory.create(UMIRead.class, reads, byName,
                new UMIReadCodec(), SAMFileWriterImpl.getDefaultMaxRecordsInRam(), null);
        final GroupingIterator<UMIRead> byNameGroups = new GroupingIterator<>(nameSorted, byName);
        return new Iterator<UMIRead>() {
            private UMIRead next = null;

            @Override
            public boolean hasNext() {
                while (next == null && byNameGroups.hasNext())
                    next = gfp.processReads(byNameGroups.next());
                return next != null;
            }

            @Override
            public UMIRead next() {
                if (!hasNext())
                    throw new NoSuchElementException();
                final UMIRead ret = next;
                next = null;
                return ret;
            }
        };
    }

    /** Reduces SAMRecords to {@link UMIRead}s, and logs progress as the SAMRecord sort would have. */
    private class SAMToUMIReadIterator extends TransformingIterator<SAMRecord, UMIRead> {
        private final boolean keepReadDetail;
        private final String geneStrandTag;
        private final String geneFunctionTag;

        SAMToUMIReadIterator(final Iterator<SAMRecord> underlyingIterator, final boolean keepReadDetail,
                             final String geneStrandTag, final String geneFunctionTag) {
            super(underlyingIterator);
            this.keepReadDetail = keepReadDetail;
            this.geneStrandTag = geneStrandTag;
            this.geneFunctionTag = geneFunctionTag;
        }

        @Override
        public UMIRead next() {
            final SAMRecord r = this.underlyingIterator.next();
            progress.record(r);
            final String cell = Utils.getCellBC(r, cellBarcodeTag);
            final String gene = geneCache.intern(r.getStringAttribute(geneTag));
            final String umi = r.getStringAttribute(molecularBarcodeTag);
            if (keepReadDetail)
                return new UMIRead(r.getReadName(), r.getReadNegativeStrandFlag(), cell, gene,
                        r.getStringAttribute(geneStrandTag), r.getStringAttribute(geneFunctionTag), umi);
            return new UMIRead(null, r.getReadNegativeStrandFlag(), cell, gene, null, null, umi);
        }
    }

    @Override
    public boolean hasNext() {
        return groups.hasNext();
    }

    /** @return all the UMIs of the next (cell, gene), or (gene, cell) when not cell first. */
    @Override
    public UMICollection next() {
        final List<UMIKeyCount> keys = groups.next();
        final UMIKeyCount first = keys.getFirst();
        final UMICollection umi = new UMICollection(first.cell, first.gene);
        // Sorting put the UMIs in alphabetical order.  Put them back in the order they first appeared in the input,
        // which is the order they were added when SAMRecords were sorted.  The UMICollection's counter is a
        // HashMap so its iteration order depends on this, and that order shows up in some output files.
        keys.sort(UMIKeyCount.firstSeenOrder());
        for (final UMIKeyCount k : keys)
            umi.incrementMolecularBarcodeCount(umiCache.intern(k.umi), k.getCount());
        return umi;
    }

    @Override
    public void close() {
        CloserUtil.close(groups);
    }
}
