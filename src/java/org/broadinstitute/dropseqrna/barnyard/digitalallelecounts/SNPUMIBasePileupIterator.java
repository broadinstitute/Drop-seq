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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.io.File;
import java.util.Collection;
import java.util.List;

import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;
import org.broadinstitute.dropseqrna.utils.MultiComparator;
import org.broadinstitute.dropseqrna.utils.ProgressLoggingIterator;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionIteratorWrapper;
import org.broadinstitute.dropseqrna.utils.readiterators.MapQualityFilteredIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.MissingTagFilteringIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.annotation.LocusFunction;


public class SNPUMIBasePileupIterator implements CloseableIterator<SNPUMIBasePileup>{
    private static final Log LOG = Log.getInstance(SNPUMIBasePileupIterator.class);
	private final GroupingIterator<SAMRecord> atoi;
	private final String geneTag;
	private final String cellBarcodeTag;
	private final String molecularBarcodeTag;
	private final String snpTag;
	private final IntervalList snpIntervals;
	private final SortOrder sortOrder;
	private final String functionTag;

	/**
	 * Construct an object that generates SNPUMIBasePileup objects from a BAM file.
	 * Each of these objects constructs a pileup of the bases and qualities
	 * @param bamFile The BAM file to extract UMIs from
	 * @param geneTag The geneExon tag on BAM records
	 * @param cellBarcodeTag The cell barcode tag on BAM records
	 * @param molecularBarcodeTag The molecular barcode tag on BAM records
	 * @param strandTag The strand tag on BAM records
	 * @param readMQ The minimum map quality of the reads
	 * @param assignReadsToAllGenes Should records tagged with multiple genes be double counted, once for each gene?
	 * @param useStrandInfo should the gene and read strand match for the read to be accepted
	 * @param cellBarcodes The list of cell barcode tag values that match the <cellBarcodeTag> tag on the BAM records.  Only reads with these values will be used.  If set to null, all cell barcodes are used.
	 */
	public SNPUMIBasePileupIterator (final File bamFile, final IntervalList snpIntervals, final String geneTag, final String geneStrandTag, final String geneFunctionTag,
									 final Collection <LocusFunction> acceptedLociFunctions, final StrandStrategy strandStrategy, final String cellBarcodeTag,
                                     final String molecularBarcodeTag, final String snpTag, final String functionTag, final int readMQ,
                                     final boolean assignReadsToAllGenes, final List<String> cellBarcodes) {
		this(bamFile, snpIntervals, geneTag, geneStrandTag, geneFunctionTag, acceptedLociFunctions, strandStrategy, cellBarcodeTag, molecularBarcodeTag, snpTag, functionTag, readMQ, assignReadsToAllGenes, cellBarcodes, SortOrder.SNP_GENE);
	}


	/**
	 * Construct an object that generates SNPUMIBasePileup objects from a BAM file.
	 * Each of these objects constructs a pileup of the bases and qualities
	 * @param bamFile The BAM file to extract UMIs from
	 * @param geneTag The geneExon tag on BAM records
	 * @param cellBarcodeTag The cell barcode tag on BAM records
	 * @param molecularBarcodeTag The molecular barcode tag on BAM records
	 * @param strandTag The strand tag on BAM records
	 * @param readMQ The minimum map quality of the reads
	 * @param assignReadsToAllGenes Should records tagged with multiple genes be double counted, once for each gene?
	 * @param useStrandInfo should the gene and read strand match for the read to be accepted
	 * @param cellBarcodes The list of cell barcode tag values that match the <cellBarcodeTag> tag on the BAM records.  Only reads with these values will be used.  If set to null, all cell barcodes are used.
	 */
	public SNPUMIBasePileupIterator (final File bamFile, final IntervalList snpIntervals, final String geneTag, final String geneStrandTag, final String geneFunctionTag,
									 final Collection <LocusFunction> acceptedLociFunctions, final StrandStrategy strandStrategy, final String cellBarcodeTag,
                                     final String molecularBarcodeTag, final String snpTag, final String functionTag, final int readMQ,
                                     final boolean assignReadsToAllGenes, final List<String> cellBarcodes, final SortOrder order) {
		
		
		// this is a big ugly.
        this(new SamHeaderAndIterator(bamFile), snpIntervals, geneTag, geneStrandTag, geneFunctionTag, acceptedLociFunctions, strandStrategy, cellBarcodeTag, molecularBarcodeTag, snpTag, functionTag, readMQ, assignReadsToAllGenes, cellBarcodes, order);        
	}
	
	/**
	 * Construct an object that generates SNPUMIBasePileup objects from a BAM file.
	 * Each of these objects constructs a pileup of the bases and qualities
	 * @param headerAndIter The SAM header and iterator to read from.
	 * @param geneTag The geneExon tag on BAM records
	 * @param cellBarcodeTag The cell barcode tag on BAM records
	 * @param molecularBarcodeTag The molecular barcode tag on BAM records
	 * @param strandTag The strand tag on BAM records
	 * @param readMQ The minimum map quality of the reads
	 * @param assignReadsToAllGenes Should records tagged with multiple genes be double counted, once for each gene?
	 * @param useStrandInfo should the gene and read strand match for the read to be accepted
	 * @param cellBarcodes The list of cell barcode tag values that match the <cellBarcodeTag> tag on the BAM records.  Only reads with these values will be used.  If set to null, all cell barcodes are used.
	 */
	public SNPUMIBasePileupIterator (final SamHeaderAndIterator headerAndIter, final IntervalList snpIntervals, final String geneTag, final String geneStrandTag, final String geneFunctionTag,
			 final Collection <LocusFunction> acceptedLociFunctions, final StrandStrategy strandStrategy, final String cellBarcodeTag,
            final String molecularBarcodeTag, final String snpTag, final String functionTag, final int readMQ,
            final boolean assignReadsToAllGenes, final List<String> cellBarcodes, final SortOrder order) {
		
		this.geneTag=geneTag;
		this.cellBarcodeTag=cellBarcodeTag;
		this.molecularBarcodeTag=molecularBarcodeTag;
		this.snpTag = snpTag;
		this.snpIntervals=snpIntervals;
		this.sortOrder=order;
		this.functionTag = functionTag;
		
		final ProgressLogger logger = new ProgressLogger(LOG);

        // add a progress logger.
        ProgressLoggingIterator loggingIterator = new ProgressLoggingIterator (headerAndIter.iterator, logger);

        // Filter records before sorting, to reduce I/O
        final MissingTagFilteringIterator filteringIterator =
                new MissingTagFilteringIterator(loggingIterator, geneTag, cellBarcodeTag, molecularBarcodeTag);

        // Filter reads on map quality
     	MapQualityFilteredIterator filteringIterator2 = new MapQualityFilteredIterator(filteringIterator, readMQ, true);

     	// Filter/assign reads based on functional annotations
     	GeneFunctionIteratorWrapper gfteratorWrapper = new GeneFunctionIteratorWrapper(filteringIterator2, geneTag, geneStrandTag, geneFunctionTag, assignReadsToAllGenes, strandStrategy, acceptedLociFunctions);

        SNPUMICellReadIteratorWrapper snpumiCellReadIterator = new SNPUMICellReadIteratorWrapper(gfteratorWrapper, snpIntervals, cellBarcodeTag, cellBarcodes, geneTag, snpTag, readMQ);

        // create comparators in the order the data should be sorted
        SAMSequenceDictionary sd = snpIntervals.getHeader().getSequenceDictionary();
        @SuppressWarnings("unchecked")

        MultiComparator<SAMRecord> multiComparator;

        if (order==SortOrder.SNP_GENE)
			multiComparator = new MultiComparator<>(
        			new IntervalTagComparator(snpTag, sd),
                    new StringTagComparator(geneTag),
                    new StringTagComparator(cellBarcodeTag),
                    new StringTagComparator(molecularBarcodeTag));
		else if (order==SortOrder.SNP_CELL)
			multiComparator = new MultiComparator<>(
        			new IntervalTagComparator(snpTag, sd),
        			new StringTagComparator(cellBarcodeTag),
                    new StringTagComparator(geneTag),
                    new StringTagComparator(molecularBarcodeTag));
		else if (order==SortOrder.CELL_SNP)
			multiComparator = new MultiComparator<>(
        			new StringTagComparator(cellBarcodeTag),
        			new IntervalTagComparator(snpTag, sd),
                    new StringTagComparator(geneTag),
                    new StringTagComparator(molecularBarcodeTag));
		else {
			multiComparator=null;
			throw new IllegalArgumentException("Sort order " + order +" unsupported");
		}
		final CloseableIterator<SAMRecord> sortingIterator =
		                SamRecordSortingIteratorFactory.create(headerAndIter.header, snpumiCellReadIterator, multiComparator, logger);

		this.atoi = new GroupingIterator<>(sortingIterator, multiComparator);
	}

	public SortOrder getSortOrder() {
		return this.sortOrder;
	}

	@Override
	public SNPUMIBasePileup next() {
		if (!this.atoi.hasNext())
			return null;

		Collection<SAMRecord> records = this.atoi.next();
		// the first read and all other reads share the same tags, so can extract these properties out of the first read.
		SAMRecord firstRead = records.iterator().next();

		SNPUMIBasePileup result = getInitialPileup(firstRead);

		// get the base and quality for each read.
		for (SAMRecord r: records) {
			result.addRead(r);
			result=addLocusFunction(result, r);
		}
		return result;
	}

	private SNPUMIBasePileup getInitialPileup (final SAMRecord rec) {
		String snpID=rec.getStringAttribute(this.snpTag);
		Interval snpInterval = IntervalTagComparator.fromString(snpID);
		String gene=rec.getStringAttribute(this.geneTag);
		String cell =rec.getStringAttribute(this.cellBarcodeTag);
		String molecularBarcode = rec.getStringAttribute(this.molecularBarcodeTag);
		SNPUMIBasePileup result =new SNPUMIBasePileup(snpInterval, gene, cell, molecularBarcode);
		result=addLocusFunction(result, rec);
		return result;
	}

	/**
	 * Since this read is being added to the pileup, record the function of the read.
	 * @param result
	 * @param rec
	 * @return
	 */
	private SNPUMIBasePileup addLocusFunction (final SNPUMIBasePileup result, final SAMRecord rec) {
		if (this.functionTag==null) return result; // no op.
		String lfValue = rec.getStringAttribute(this.functionTag);
		if (lfValue!=null) {
			LocusFunction lf = LocusFunction.valueOf(lfValue);
			result.addLocusFunction(lf);
		}
		return result;
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


	@Override
	public String toString () {
		return "";
	}


	public String getGeneTag() {
		return geneTag;
	}


	public String getCellBarcodeTag() {
		return cellBarcodeTag;
	}


	public String getMolecularBarcodeTag() {
		return molecularBarcodeTag;
	}

	public String getSnpTag() {
		return snpTag;
	}

	public IntervalList getSNPIntervals () {
		return this.snpIntervals;
	}



}
