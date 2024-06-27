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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.*;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorStrategy;
import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.*;
import picard.annotation.LocusFunction;

import java.util.*;
import java.util.stream.StreamSupport;

public class UMIIterator implements CloseableIterator<UMICollection>  {

    private static final Log log = Log.getInstance(UMIIterator.class);
    private static final ProgressLogger prog = new ProgressLogger(log);

	private final GroupingIterator<SAMRecord> atoi;
	private final String geneTag;
	private final String cellBarcodeTag;
	private final String molecularBarcodeTag;
    private final StringInterner stringCache = new StringInterner();
    private final Set<String> cellBarcodesSeen;
    private final boolean retainReads;
	/**
	 * Construct an object that generates UMI objects from a BAM file
     * @param headerAndIterator The BAM records to extract UMIs from
	 * @param geneTag The gene tag on BAM records
	 * @param cellBarcodeTag The cell barcode tag on BAM records
	 * @param molecularBarcodeTag The molecular barcode tag on BAM records
	 * @param geneStrandTag The strand tag on BAM records
	 * @param readMQ The minimum map quality of the reads
	 * @param assignReadsToAllGenes Should records tagged with multiple genes be double counted, once for each gene?
	 * @param strandStrategy should the gene and read strand match for the read to be accepted
	 * @param cellBarcodes The list of cell barcode tag values that match the <cellBarcodeTag> tag on the BAM records.
     *                     Only reads with these values will be used.  If set to null, all cell barcodes are used.
	 */
	public UMIIterator(final SamHeaderAndIterator headerAndIterator,
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
                       final Collection<String> cellBarcodes) {
		this(headerAndIterator, geneTag, geneStrandTag, geneFunctionTag, strandStrategy, acceptedLociFunctions, functionStrategy,
				cellBarcodeTag, molecularBarcodeTag, readMQ, assignReadsToAllGenes, cellBarcodes, false, false);
	}

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
	 */
	public UMIIterator(final SamHeaderAndIterator headerAndIterator,
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
					   final boolean recordCellsInInput) {
		
		this(headerAndIterator, geneTag, geneStrandTag, geneFunctionTag, strandStrategy, acceptedLociFunctions, functionStrategy,
				cellBarcodeTag, molecularBarcodeTag, readMQ, assignReadsToAllGenes, cellBarcodes, cellFirstSort, recordCellsInInput, false,
				null);
		
	}
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
	 * @param retainReads If true the SAMRecords are added to the UMICollection.  This uses more memory than the standard UMICollection so exercise caution.  
	 * This is false in other method signatures by default.  If false, reads can be simplified for faster serialization.
	 * 
	 */
	public UMIIterator(final SamHeaderAndIterator headerAndIterator,
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
					   final boolean retainReads,
					   final IntervalList intervals) {
		
		this.retainReads=retainReads;
        this.geneTag=geneTag;
        this.cellBarcodeTag=cellBarcodeTag;
        this.molecularBarcodeTag=molecularBarcodeTag;

        final StringTagComparator cellBarcodeTagComparator = new StringTagComparator(cellBarcodeTag);
        final StringTagComparator geneExonTagComparator = new StringTagComparator(geneTag);
        final MultiComparator<SAMRecord> multiComparator;
		if (cellFirstSort)
			multiComparator = new MultiComparator<>(cellBarcodeTagComparator, geneExonTagComparator);
		else
			multiComparator = new MultiComparator<>(geneExonTagComparator, cellBarcodeTagComparator);
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
					new TagValueFilteringIterator<>(filteringIterator, this.cellBarcodeTag, cellBarcodes);
		}
		if (intervals != null) {
			filteringIterator = new IntervalFilteringIterator(filteringIterator, intervals, true);
		}

		// Filter/assign reads based on functional annotations
		GeneFunctionIteratorWrapper wrapper = new GeneFunctionIteratorWrapper(filteringIterator, geneTag,
				geneStrandTag, geneFunctionTag, assignReadsToAllGenes, strandStrategy, acceptedLociFunctions, functionStrategy);
		Iterator<SAMRecord> samRecordIter = wrapper;

		// Strip down the reads to a more minimal set of TAGS, set reads to be empty on request
		if (!retainReads) {
			List<String> requiredTags = Arrays.asList(geneTag, geneStrandTag, geneFunctionTag, cellBarcodeTag, molecularBarcodeTag);
			samRecordIter = new SimplifySAMRecordIterator(samRecordIter, requiredTags);
		}

		// if map quality < unique handle low map quality reads.
		if (readMQ <=3) {
			samRecordIter= processLowMapQuality (samRecordIter, headerAndIterator.header, wrapper.getGeneFunctionProcessor());
		}

		// now sort by bam tags.
		CloseableIterator<SAMRecord> sortedAlignmentIterator = SamRecordSortingIteratorFactory.create(
				headerAndIterator.header, samRecordIter, multiComparator, prog);

		// Not really -- merge sort is ongoing.
        log.info("Sorting finished.");

		// get reads from the sink, sort by multiComparator, group
		this.atoi = new GroupingIterator<>(sortedAlignmentIterator, multiComparator);

	}

	/**
	 * If map quality includes non-unique reads:
	 * 1. sort the reads in query name order
	 * 2. filter multimapping reads that are map to more than one gene
	 * 3. return an interator for further sorting.
	 * @return An iterator of SAMRecord objects
	 */
	private Iterator<SAMRecord> processLowMapQuality (Iterator<SAMRecord> readIter, SAMFileHeader header, GeneFunctionProcessor gfp) {
		SimpleQueryNameComparator comp = new SimpleQueryNameComparator();

		CloseableIterator<SAMRecord> queryNameSortedIterator = SamRecordSortingIteratorFactory.create(
				header, readIter, comp, prog);

		// now that the data is sorting in read name order you can group and process.
		GroupingIterator<SAMRecord> qnGroup = new GroupingIterator<>(queryNameSortedIterator, comp);
		Iterator<List<SAMRecord>> iter = new MultiMapFilteringIterator(qnGroup, gfp);
		Iterator<SAMRecord> flatIterator = StreamSupport.stream(Spliterators.spliteratorUnknownSize(iter, Spliterator.ORDERED), false).flatMap(Collection::stream).iterator();
		return flatIterator;
	}



    /**
	 * Gets the next UMI Collection - all the molecular barcodes and reads on those for a particular gene/cell.
	 * @return Null if there are no reads left in the iterator.  Otherwise, returns a UMICollection.
	 */
	@Override
	public UMICollection next () {
		if (!this.atoi.hasNext())
			return null;

		Collection<SAMRecord> records = this.atoi.next();
		PeekableIterator<SAMRecord> recordCollectionIter = new PeekableIterator<>(records.iterator());

		// the next record is the first of the "batch"
		SAMRecord r = recordCollectionIter.peek();
		String currentGene = r.getStringAttribute(geneTag);
		String currentCell = Utils.getCellBC(r, cellBarcodeTag);
		UMICollection umi = new UMICollection(currentCell, currentGene);

		while (recordCollectionIter.hasNext()) {
			// if there's a next read, set up the current gene/cell variables.
			r=recordCollectionIter.next();
            // Cache the UMI.  Note that within a single (cell, gene), the ObjectCounter effectively caches the UMI string,
            // but they are not cached *across* (cell, gene), so this should help.
			String molecularBarcode = stringCache.intern(r.getStringAttribute(this.molecularBarcodeTag));
			umi.incrementMolecularBarcodeCount(molecularBarcode);
			// optionally store the reads for a cell/gene 
			if (retainReads)
				umi.addRead(molecularBarcode, r);
			
		}
		// I ran out of reads
		recordCollectionIter.close();
		return (umi);
	}

	@Override
	public void remove() {
		this.atoi.remove();
	}

	@Override
	public void close() {
        CloserUtil.close(this.atoi);
	}

	@Override
	public boolean hasNext() {
		return this.atoi.hasNext();
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

	/**
	 * Given an iterator that serves reads grouped by query name, process reads that map multiple times.
	 * This tests if the interpreted GTF tags map to one gene, or more than one.  If only one gene, then
	 * that read is selected for downstream analysis.  If reads are ambiguous, all reads are rejected.
	 * <p>
	 * This should enable capture of expression for reads that map to multiple locations, but only once to a gene.
	 * This replicates how expression is counted in STARSolo.
	 */
	// change to TransformingIterator
	private static class MultiMapFilteringIterator extends CountChangingIteratorWrapper<List<SAMRecord>> {

		private final GeneFunctionProcessor gfp;
		protected MultiMapFilteringIterator(Iterator<List<SAMRecord>> underlyingIterator, GeneFunctionProcessor gfp) {
			super(underlyingIterator);
			this.gfp = gfp;
		}

		@Override
		protected void processRecord(List<SAMRecord> rec) {
			SAMRecord result = gfp.processReads(rec);
			if (result!=null)
				queueRecordForOutput(Collections.singletonList(result));
		}

	}

	private static class SimplifySAMRecordIterator extends TransformingIterator<SAMRecord, SAMRecord> {

		private final Set<String> requiredTags;
		public SimplifySAMRecordIterator(Iterator<SAMRecord> underlyingIterator, Collection<String> requiredTags) {
			super(underlyingIterator);
			this.requiredTags = new HashSet<>(requiredTags);
		}

		@Override
		public SAMRecord next() {
			SAMRecord r= this.underlyingIterator.next();
			for (SAMRecord.SAMTagAndValue v: r.getAttributes()) {
				if(!requiredTags.contains(v.tag))
					r.setAttribute(v.tag, null);
			}
			r.setBaseQualities(SAMRecord.NULL_QUALS);
			r.setReadBases(SAMRecord.NULL_SEQUENCE);
			r.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
			return r;
		}
	}

	private static class SimpleQueryNameComparator implements Comparator<SAMRecord> {
		SAMRecordQueryNameComparator comp = new SAMRecordQueryNameComparator();
		@Override
		public int compare(SAMRecord o1, SAMRecord o2) {
			return comp.fileOrderCompare(o1, o2);
		}
	}
}
