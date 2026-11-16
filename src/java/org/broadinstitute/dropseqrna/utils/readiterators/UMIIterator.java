/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.*;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorStrategy;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.*;
import org.broadinstitute.dropseqrna.utils.readiterators.umisort.ReadIntervalSummary;
import org.broadinstitute.dropseqrna.utils.readiterators.umisort.SortedUMICollectionIterator;
import picard.annotation.LocusFunction;

import java.util.*;

/**
 * Generates one {@link UMICollection} per (cell, gene), or (gene, cell), from a BAM file.  Reads are filtered, reduced to
 * (cell, gene, UMI) and sorted by {@link SortedUMICollectionIterator}; whole SAMRecords are never sorted.
 */
public class UMIIterator implements CloseableIterator<UMICollection>  {

    private static final Log log = Log.getInstance(UMIIterator.class);
    private static final ProgressLogger prog = new ProgressLogger(log);

	private final SortedUMICollectionIterator umiKeyIter;
	private final String cellBarcodeTag;
    private final Set<String> cellBarcodesSeen;
    // The records being read.  Closed with this iterator: a multi-threaded reader holds decode threads until it is closed.
    private final CloseableIterator<SAMRecord> source;

	/**
	 * Construct an object that generates UMI objects from a BAM file
	 * @param headerAndIterator The BAM records to extract UMIs from
	 * @param geneTag The geneExon tag on BAM records
	 * @param cellBarcodeTag The cell barcode tag on BAM records
	 * @param molecularBarcodeTag The molecular barcode tag on BAM records
	 * @param geneStrandTag The strand tag on BAM records
	 * @param readMQ The minimum map quality of the reads
	 * @param assignReadsToAllGenes Should records tagged with multiple genes be double counted, once for each gene?
	 * @param strandStrategy should the gene and read strand match for the read to be accepted
	 * @param cellBarcodes The list of cell barcode tag values that match the <cellBarcodeTag> tag on the BAM records.
     *                     Only reads with these values will be used.  If set to null, all cell barcodes are used.
	 * @param cellFirstSort if true, then cell barcodes are sorted first, followed by gene/exon tags.
     *                      If false, then gene/exon tags are sorted first, followed by cells.  false is the default and used in the other constructor.
	 * @param recordCellsInInput While sorting the input, keep track of what cells appear in the input.  This record
	 *                           is not complete until iteration is started.
	 * @param collectReadIntervals If true, each UMI of the UMICollections carries a {@link ReadIntervalSummary} of
	 *                             where its reads aligned (see {@link UMICollection#getReadInterval}).
	 */
	private UMIIterator(final SamHeaderAndIterator headerAndIterator,
					   final String geneTag,
                       final String geneStrandTag,
                       final String geneFunctionTag,
                       final StrandStrategy strandStrategy,
                       final Collection <LocusFunction> acceptedLociFunctions,
					   FunctionalDataProcessorStrategy functionStrategy,
                       final String cellBarcodeTag,
                       final String molecularBarcodeTag,
                       final int readMQ,
                       final boolean assignReadsToAllGenes,
                       final Collection<String> cellBarcodes,
                       final boolean cellFirstSort,
					   final boolean recordCellsInInput,
					   final boolean collectReadIntervals,
					   final IntervalList intervals) {

        this.source=headerAndIterator.iterator;
        this.cellBarcodeTag=cellBarcodeTag;

        // Filter records before sorting, to reduce I/O
		Iterator<SAMRecord> filteringIterator =
                new MissingTagFilteringIterator(headerAndIterator.iterator, cellBarcodeTag, geneTag, molecularBarcodeTag);

		// Filter out reads that STARsolo has marked as chimeric
		filteringIterator = new STARSoloChimericReadFilteringIterator(filteringIterator, molecularBarcodeTag);

		// Filter reads on map quality.  Optionally keep non-primary reads at low map quality.
		boolean rejectNonPrimaryReads = readMQ>3;
		if (!rejectNonPrimaryReads)
			log.info("Detected map quality threshold <=3, retaining non primary reads");

		filteringIterator = new MapQualityFilteredIterator(filteringIterator, readMQ, rejectNonPrimaryReads);

		if (recordCellsInInput) {
			this.cellBarcodesSeen = new HashSet<>();
			filteringIterator = new CellBarcodeRecorder(filteringIterator);
		} else {
			this.cellBarcodesSeen = null;
		}

		// Filter reads on if the read contains a cell barcode, if cell barcodes have been specified.
		if (cellBarcodes != null) {
			filteringIterator =
					new TagValueFilteringIterator<>(filteringIterator, cellBarcodeTag, cellBarcodes);
		}
		if (intervals != null) {
			filteringIterator = new IntervalFilteringIterator(filteringIterator, intervals, true);
		}

		// Filter/assign reads based on functional annotations
		GeneFunctionIteratorWrapper wrapper = new GeneFunctionIteratorWrapper(filteringIterator, geneTag,
				geneStrandTag, geneFunctionTag, assignReadsToAllGenes, strandStrategy, acceptedLociFunctions, functionStrategy);

		// Multi-mapping reads (readMQ <= 3) are resolved on the compact reads, not by sorting SAMRecords by name.
		this.umiKeyIter = new SortedUMICollectionIterator(wrapper, wrapper.getGeneFunctionProcessor(),
				cellBarcodeTag, geneTag, geneStrandTag, geneFunctionTag, molecularBarcodeTag, readMQ, cellFirstSort,
				collectReadIntervals, prog);
	}

    /**
	 * Gets the next UMI Collection - all the molecular barcodes for a particular gene/cell.
	 * @return Null if there are no reads left in the iterator.  Otherwise, returns a UMICollection.
	 */
	@Override
	public UMICollection next () {
		return this.umiKeyIter.hasNext() ? this.umiKeyIter.next() : null;
	}

	@Override
	public void remove() {
		this.umiKeyIter.remove();
	}

	@Override
	public void close() {
        CloserUtil.close(this.umiKeyIter);
        CloserUtil.close(this.source);
	}

	@Override
	public boolean hasNext() {
		return this.umiKeyIter.hasNext();
	}

	public Set<String> getCellBarcodesSeen() {
		return cellBarcodesSeen;
	}

	private class CellBarcodeRecorder
	implements Iterator<SAMRecord> {
		private final Iterator<SAMRecord> underlyingIterator;
		private final short cellBarcodeBinaryTag;

		public CellBarcodeRecorder(Iterator<SAMRecord> underlyingIterator) {
			this.underlyingIterator = underlyingIterator;
			cellBarcodeBinaryTag = SAMTag.makeBinaryTag(cellBarcodeTag);
		}

		@Override
		public boolean hasNext() {
			return underlyingIterator.hasNext();
		}

		@Override
		public SAMRecord next() {
			final SAMRecord ret = underlyingIterator.next();
			cellBarcodesSeen.add((String)ret.getAttribute(cellBarcodeBinaryTag));
			return ret;
		}
	}

	public static class UMIIteratorBuilder {
		// required parameters
		private final SamHeaderAndIterator headerAndIterator;
		private final String geneTag;
		private final String geneStrandTag;
		private final String geneFunctionTag;
		private final StrandStrategy strandStrategy;
		private final Collection <LocusFunction> acceptedLociFunctions;
		private final FunctionalDataProcessorStrategy functionStrategy;
		private final String cellBarcodeTag;
		private final String molecularBarcodeTag;
		private final int readMQ;

		// parameters with default values
		private boolean assignReadsToAllGenes=false;
		private boolean cellFirstSort=false;
		private boolean recordCellsInInput=false;
		private boolean collectReadIntervals=false;

		// nullable parameters
		private Collection<String> cellBarcodes;
		private IntervalList intervals;

		public UMIIteratorBuilder(SamHeaderAndIterator headerAndIterator, String geneTag, String geneStrandTag,
								  String geneFunctionTag, StrandStrategy strandStrategy,
								  Collection<LocusFunction> acceptedLociFunctions,
								  FunctionalDataProcessorStrategy functionStrategy, String cellBarcodeTag,
								  String molecularBarcodeTag, int readMQ) {
			this.headerAndIterator = headerAndIterator;
			this.geneTag = geneTag;
			this.geneStrandTag = geneStrandTag;
			this.geneFunctionTag = geneFunctionTag;
			this.strandStrategy = strandStrategy;
			this.acceptedLociFunctions = acceptedLociFunctions;
			this.functionStrategy = functionStrategy;
			this.cellBarcodeTag = cellBarcodeTag;
			this.molecularBarcodeTag = molecularBarcodeTag;
			this.readMQ = readMQ;
		}

		public boolean isAssignReadsToAllGenes() {
			return assignReadsToAllGenes;
		}

		public UMIIteratorBuilder assignReadsToAllGenes(boolean assignReadsToAllGenes) {
			this.assignReadsToAllGenes = assignReadsToAllGenes;
			return this;
		}

		public boolean isCellFirstSort() {
			return cellFirstSort;
		}

		public UMIIteratorBuilder cellFirstSort(boolean cellFirstSort) {
			this.cellFirstSort = cellFirstSort;
			return this;
		}

		public boolean isRecordCellsInInput() {
			return recordCellsInInput;
		}

		public UMIIteratorBuilder recordCellsInInput(boolean recordCellsInInput) {
			this.recordCellsInInput = recordCellsInInput;
			return this;
		}

		public boolean isCollectReadIntervals() {
			return collectReadIntervals;
		}

		/**
		 * @param collectReadIntervals If true, each UMI of the UMICollections carries a {@link ReadIntervalSummary} of
		 *                             where its reads aligned.  Off by default, as it makes the sort larger.
		 */
		public UMIIteratorBuilder collectReadIntervals(boolean collectReadIntervals) {
			this.collectReadIntervals = collectReadIntervals;
			return this;
		}

		public Collection<String> getCellBarcodes() {
			return cellBarcodes;
		}

		public UMIIteratorBuilder setCellBarcodes(Collection<String> cellBarcodes) {
			this.cellBarcodes = cellBarcodes;
			return this;
		}

		public IntervalList getIntervals() {
			return intervals;
		}

		public UMIIteratorBuilder setIntervals(IntervalList intervals) {
			this.intervals = intervals;
			return this;
		}

		public UMIIterator build() {
			return new UMIIterator(headerAndIterator, geneTag, geneStrandTag, geneFunctionTag, strandStrategy,
					acceptedLociFunctions, functionStrategy, cellBarcodeTag, molecularBarcodeTag, readMQ,
					assignReadsToAllGenes, cellBarcodes, cellFirstSort, recordCellsInInput, collectReadIntervals, intervals);
		}
	}
}
