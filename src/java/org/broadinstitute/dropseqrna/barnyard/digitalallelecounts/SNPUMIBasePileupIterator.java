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

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorStrategy;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;
import org.broadinstitute.dropseqrna.utils.MultiComparator;
import org.broadinstitute.dropseqrna.utils.ProgressLoggingIterator;
import org.broadinstitute.dropseqrna.utils.SamWriterSink;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
import org.broadinstitute.dropseqrna.utils.readiterators.*;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
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
	private SamWriterSink sink;
	
	
	/**
	 * Construct an object that generates SNPUMIBasePileup objects from a BAM file.
	 * Each of these objects constructs a pileup of the bases and qualities
	 * @param headerAndIter The SAM header and iterator to read from.
	 * @param geneTag The geneExon tag on BAM records
	 * @param cellBarcodeTag The cell barcode tag on BAM records
	 * @param molecularBarcodeTag The molecular barcode tag on BAM records
	 * @param geneStrandTag The strand tag on BAM records
	 * @param readMQ The minimum map quality of the reads
	 * @param assignReadsToAllGenes Should records tagged with multiple genes be double counted, once for each gene?
	 * @param cellBarcodes The list of cell barcode tag values that match the <cellBarcodeTag> tag on the BAM records.  Only reads with these values will be used.  If set to null, all cell barcodes are used.
	 * @param meanGenotypeQuality A map from SNP Intervals to the mean genotype quality.  This is an optional parameter.  If supplied, in cases where a read has 
	 * multiple SNPs, the SNP with the highest average genotype quality will be selected to avoid read "double counting".
	 * @param order How the reads groups will be ordered by this iterator
	 * @param failFastThreshold If the iterator sees at least this many UMIs without encountering a transcribed SNP the iterator will fail.  Set to -1 to disable
	 */
	public SNPUMIBasePileupIterator (final SamHeaderAndIterator headerAndIter, final IntervalList snpIntervals, final String geneTag, final String geneStrandTag, final String geneFunctionTag,
			 final Collection <LocusFunction> acceptedLociFunctions, final StrandStrategy strandStrategy, FunctionalDataProcessorStrategy functionStrategy,
									 final String cellBarcodeTag, final String molecularBarcodeTag, final String snpTag, final String functionTag, final int readMQ,
            final boolean assignReadsToAllGenes, final List<String> cellBarcodes, final Map<Interval, Double> meanGenotypeQuality, final SortOrder order,
            int failFastThreshold) {
				
		this.geneTag=geneTag; 
		this.cellBarcodeTag=cellBarcodeTag;
		this.molecularBarcodeTag=molecularBarcodeTag;
		this.snpTag = snpTag;
		this.snpIntervals=snpIntervals;
		this.sortOrder=order;
		this.functionTag = functionTag;		
		
		final ProgressLogger logger = new ProgressLogger(LOG);

        // add a progress logger.
        final ProgressLoggingIterator loggingIterator = new ProgressLoggingIterator (headerAndIter.iterator, logger);

        // filter out to contigs in the interval list.
        Set<String> contigsToFilter = snpIntervals.getIntervals().stream().map(x -> x.getContig()).collect(Collectors.toSet());
        ChromosomeFilteringIterator cfi = new ChromosomeFilteringIterator(loggingIterator, contigsToFilter, false);
        
        // Filter reads on map quality
     	MapQualityFilteredIterator filteringIterator = new MapQualityFilteredIterator(cfi, readMQ, true);

        // Filter records before sorting, to reduce I/O
        final MissingTagFilteringIterator filteringIterator2 =
                new MissingTagFilteringIterator(filteringIterator, geneTag, cellBarcodeTag, molecularBarcodeTag);

		final STARSoloChimericReadFilteringIterator chimericReadFilteringIterator = new STARSoloChimericReadFilteringIterator(filteringIterator2, molecularBarcodeTag);
     	// Filter/assign reads based on functional annotations
     	GeneFunctionIteratorWrapper gfteratorWrapper = new GeneFunctionIteratorWrapper(chimericReadFilteringIterator, geneTag, geneStrandTag, geneFunctionTag, assignReadsToAllGenes, strandStrategy, acceptedLociFunctions, functionStrategy);
     	
     	/**
     	 * In this section, data is sorted into the standard cell / gene / umi order and grouped so that all reads for a single UMI can be evaluated at once.
     	 * After evaluation by SNPUMICellReadIteratorWrapper2, the groups of reads (List<SAMRECORD>) are flattened back to an iterator of SAMRecords,
     	 * which is then pushed into the next sorting collection.
     	 * 
     	 * An alternative to this process would be to gather all reads that are within a 10kb window, and adjudicate which reads should be tagged with SNPs such 
     	 * that each UMI only sees a single SNP.  This would happen with a sliding window where as each SNP fell "out" of the window it would be written.  Each SNP
     	 * would memorize the cell/gene/umi of the reads, and look for conflicts with new incoming data.  
     	 * That would prevent having to write records to disk to sort them, but would be far more complex to code.  Edge cases where there were long introns 
     	 * would likely prevent all SNPs from being filtered, but might be a faster approach. 
     	 */
     	
     	
     	MultiComparator<SAMRecord> cellGeneMolBcComparator=new MultiComparator<>(
     			new StringTagComparator(cellBarcodeTag),
                new StringTagComparator(geneTag),
                new StringTagComparator(molecularBarcodeTag));
     	
     	final CloseableIterator<SAMRecord> sortingIterator1 =
                SamRecordSortingIteratorFactory.create(headerAndIter.header, gfteratorWrapper, cellGeneMolBcComparator, logger);
     	
     	GroupingIterator<SAMRecord> cellGeneMolBCGroupingIter= new GroupingIterator<>(sortingIterator1, cellGeneMolBcComparator);
     	// retaining the reference to this to modify the fail fast threshold.
     	SNPUMICellReadIteratorWrapper2 snpumiCellReadIterator = new SNPUMICellReadIteratorWrapper2(cellGeneMolBCGroupingIter, snpIntervals, cellBarcodeTag, cellBarcodes, geneTag, snpTag, readMQ, meanGenotypeQuality);
     	snpumiCellReadIterator.setFailFastThreshold(failFastThreshold);
     	
     	Iterator<SAMRecord> flattenedIterator = StreamSupport.stream(snpumiCellReadIterator.spliterator(), false).flatMap(x->x.stream()).iterator();     	          	
     	     	
     	/**
     	 * Return to the usual sorting by SNP/etc to generate pileups.
     	 */
     	     	
        // If you comment out the block above and change the sorting iterator below to use this iter, you can swap back to the "old" code.
     	// SNPUMICellReadIteratorWrapper snpumiCellReadIterator = new SNPUMICellReadIteratorWrapper(gfteratorWrapper, snpIntervals, cellBarcodeTag, cellBarcodes, geneTag, snpTag, readMQ, meanGenotypeQuality);
     	
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
		                SamRecordSortingIteratorFactory.create(headerAndIter.header, flattenedIterator, multiComparator, logger);

		this.atoi = new GroupingIterator<>(sortingIterator, multiComparator);
	}
	
	public SNPUMIBasePileupIterator (final SamHeaderAndIterator headerAndIter, final IntervalList snpIntervals, final String geneTag, final String geneStrandTag, final String geneFunctionTag,
			 final Collection <LocusFunction> acceptedLociFunctions, final StrandStrategy strandStrategy, FunctionalDataProcessorStrategy functionStrategy,
									 final String cellBarcodeTag,
           final String molecularBarcodeTag, final String snpTag, final String functionTag, final int readMQ,
           final boolean assignReadsToAllGenes, final List<String> cellBarcodes, final Map<Interval, Double> meanGenotypeQuality, final SortOrder order) {
		
		this(headerAndIter, snpIntervals, geneTag, geneStrandTag, geneFunctionTag, acceptedLociFunctions, strandStrategy, functionStrategy, cellBarcodeTag, molecularBarcodeTag, snpTag, functionTag, readMQ,
				assignReadsToAllGenes, cellBarcodes, meanGenotypeQuality, order, -1);
	}
	
	
	
	public void addReadSink (SamWriterSink sink) {
		this.sink=sink;
	}
	
	public void closeSink () {
		if (this.sink!=null) {
			this.sink.writer.close();
		}
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
			// if there's a BAM sink, add the informative read.
			if (this.sink!=null) 
				sink.add(r); 
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
		this.closeSink();
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
