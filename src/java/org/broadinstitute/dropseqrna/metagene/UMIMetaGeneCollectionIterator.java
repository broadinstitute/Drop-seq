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
package org.broadinstitute.dropseqrna.metagene;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorStrategy;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.MultiComparator;
import org.broadinstitute.dropseqrna.utils.ObjectSink;
import org.broadinstitute.dropseqrna.utils.SamWriterSink;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionProcessor;
import org.broadinstitute.dropseqrna.utils.readiterators.MapQualityPredicate;
import org.broadinstitute.dropseqrna.utils.readiterators.ReadEditDistancePredicate;
import org.broadinstitute.dropseqrna.utils.readiterators.RequiredTagPredicate;
import org.broadinstitute.dropseqrna.utils.readiterators.RequiredTagStringValuePredicate;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.annotation.LocusFunction;

/**
 * Gathers a group of reads for a cell barcode and UMI barcode across all genes.
 * @author nemesh
 *
 */
public class UMIMetaGeneCollectionIterator implements CloseableIterator<UMIMetaGeneCollection>{

	private static final Log log = Log.getInstance(UMIMetaGeneCollectionIterator.class);

	private final String cellBarcodeTag;
	private final String molecularBarcodeTag;
	private final String geneTag;
	private final GroupingIterator<SAMRecord> atoi;
	private final int mapQualityUnique;
	
	// Filtering iterators used by both modes of iterator.
	
	private final GeneFunctionProcessor gfteratorWrapper;    
    private final RequiredTagPredicate missingTagPredicate;
	private final MapQualityPredicate mapQualityPredicate;
	private final RequiredTagStringValuePredicate cellBarcodePredicate;
	private final ReadEditDistancePredicate editDistancePredicate;
	
    private final SamWriterSink sink;
    
    /**
	 * Construct an object that generates UMI objects from a BAM file for informative reads.
	 * Informative reads are those reads which have cell barcodes from the input set, have appropriate tags [cell, molecular, gene name/strand/function tags].
	 * Additionally, informative reads pass the gene function strategy passed in, having an acceptable strand,gene function, and unambiguous assignment.
	 * 
	 * If a sink is passed in, uninformative reads are added to the sink.
	 * 
	 * @param headerAndIterator The BAM records to extract UMIs from
	 * @param geneTag The gene name tag on BAM records
	 * @param geneStrandTag The gene strand tag
	 * @param geneFunctionTag The gene function tag
	 * @param strandStrategy The strand strategy to match the gene strand to the read strand.
	 * @param acceptedLociFunctions The list of gene function types that are accepted.
	 * @param cellBarcodeTag The cell barcode tag on BAM records
	 * @param molecularBarcodeTag The molecular barcode tag on BAM records
	 * @param minMapQuality The minimum read map quality of the reads.  This controls how many secondary reads are considered.  For STAR, a setting of 3 includes pairs of reads [1 secondary], 
	 * and a setting of 2 includes read triplets [2 secondary].
	 * @param mapQualityUnique The read map quality that defines a read as uniquely mapped.
	 * @param maxEditDistance THe maximum edit distance from the aligned read to the reference genome as computed by the aligner.
	 * @param cellBarcodes The list of cell barcode tag values that match the cellBarcodeTag tag on the BAM records.
     *                     Only reads with these values will be used.  If set to null, all cell barcodes are used.
     * @param sink If a sink is set with a non-null value, uninformative reads are added to the sink.  Otherwise they are discarded.                     
     * @param assignReadsToAllGenes Should records tagged with multiple genes be double counted, once for each gene?  Useful to discover ambiguous gene models, but can confuse metagene discovery where
     * reads that map uniquely to multiple genes can confuse UMI assignment.                     
	 */
		
	public UMIMetaGeneCollectionIterator(final SamHeaderAndIterator headerAndIterator, 
			final String geneTag,
            final String geneStrandTag,
            final String geneFunctionTag,
            final StrandStrategy strandStrategy,
            final Collection <LocusFunction> acceptedLociFunctions,
			final FunctionalDataProcessorStrategy functionStrategy,
            final String cellBarcodeTag,
            final String molecularBarcodeTag, final int minMapQuality, final int mapQualityUnique, final Integer maxEditDistance,
            final List<String> cellBarcodes, final SamWriterSink sink, final boolean assignReadsToAllGenes) {

		this.cellBarcodeTag=cellBarcodeTag;
		this.molecularBarcodeTag=molecularBarcodeTag;
		this.mapQualityUnique=mapQualityUnique;
		this.geneTag=geneTag;
		this.sink=sink;
		
		final StringTagComparator cellBarcodeTagComparator = new StringTagComparator(cellBarcodeTag);
		final StringTagComparator molBarcodeTagComparator = new StringTagComparator(molecularBarcodeTag);
		final MultiComparator<SAMRecord> multiComparator = new MultiComparator<>(cellBarcodeTagComparator, molBarcodeTagComparator);

		ProgressLogger pl = new ProgressLogger(this.log);

		CloseableIterator<SAMRecord> sortedAlignmentIterator;

		this.missingTagPredicate = new RequiredTagPredicate (cellBarcodeTag, molecularBarcodeTag, geneTag);
		this.mapQualityPredicate = new MapQualityPredicate(minMapQuality, false);
		this.cellBarcodePredicate = new RequiredTagStringValuePredicate(cellBarcodeTag, cellBarcodes, false);
		this.editDistancePredicate = new ReadEditDistancePredicate(maxEditDistance);
		InformativeReadFilter irf = new InformativeReadFilter(headerAndIterator.iterator, sink, mapQualityPredicate, missingTagPredicate, cellBarcodePredicate, editDistancePredicate);
		
		gfteratorWrapper = new GeneFunctionProcessor(geneTag, geneStrandTag, geneFunctionTag, assignReadsToAllGenes, strandStrategy, acceptedLociFunctions, functionStrategy);
		// need to build an informative reads FilteredIterator that gets the sink passed to it.		
		sortedAlignmentIterator = SamRecordSortingIteratorFactory.create(headerAndIterator.header, irf.iterator(), multiComparator, pl);

		this.atoi = new GroupingIterator<>(sortedAlignmentIterator, multiComparator);
	}

	/**
	 * Get a collection of queryname sorted reads encapsulated in an object designed to help make decisions
	 * about which genes are meta genes.
	 * 
	 */
	@Override
	public UMIMetaGeneCollection next () {		
		// keep pulling read sets until you have informative reads
		while (this.atoi.hasNext()) {
			Set<SAMRecord> allReads = new HashSet<> (this.atoi.next());
			Set<SAMRecord> informativeReads = getRecordsPassingFilters(allReads);

			// allReads - informative reads = uninformative reads.
			allReads.removeAll(informativeReads);
			
			// sink uninformative reads to the writer if it isn't null
			if (sink!=null) allReads.stream().forEach(x -> sink.add(x));
			// informative reads found!
			if (informativeReads.size()>0) {
				UMIMetaGeneCollection c = new UMIMetaGeneCollection(this.geneTag, this.molecularBarcodeTag, this.cellBarcodeTag, informativeReads, this.mapQualityUnique);
				c.populateGeneSets();
				return (c);						
			}				
		}
		
		// no reads left, return null.
		return null;						 		 
	} 

	private Set<SAMRecord> getRecordsPassingFilters (final Collection <SAMRecord> recs) {		
		Set<SAMRecord> result = recs.stream()				
				.map(x -> this.gfteratorWrapper.processRead(x))
				.flatMap(list -> list.stream())
				.collect(Collectors.toSet());
						
		return result;
	}
	
	/**
	 * 
	 * @author nemesh
	 *
	 */
	private class InformativeReadFilter extends FilteredIterator<SAMRecord> {
    	private final MapQualityPredicate mapQualityPredicate;
    	private final RequiredTagPredicate requiredTagPredicate;
    	private final RequiredTagStringValuePredicate cellBarcodePredicate;
    	private final ReadEditDistancePredicate editDistancePredicate;    	
    	
    	/**
    	 * Filter reads that don't pass predicates.  If the filteredReadSink is not null, push filtered into the sink.
    	 * 
    	 * @param underlyingIterator The SAMRecord iterator to filter
    	 * @param filteredReadSink If not null, filtered reads are sunk here.
    	 * @param mapQualityPredicate A filter on map quality
    	 * @param requiredTagPredicate A filter on required BAM tags
    	 * @param cellBarcodePredicate A filter on cell barcodes
    	 * @param editDistancePredicate A filter on edit distance.
    	 */
		protected InformativeReadFilter(Iterator<SAMRecord> underlyingIterator, ObjectSink<SAMRecord> filteredReadSink, 
				MapQualityPredicate mapQualityPredicate, RequiredTagPredicate requiredTagPredicate, 
				RequiredTagStringValuePredicate cellBarcodePredicate, ReadEditDistancePredicate editDistancePredicate) {
			super(underlyingIterator, filteredReadSink);
			this.mapQualityPredicate=mapQualityPredicate;
			this.requiredTagPredicate=requiredTagPredicate;
			this.cellBarcodePredicate=cellBarcodePredicate;
			this.editDistancePredicate=editDistancePredicate;						
		}

		@Override
		public boolean filterOut(SAMRecord rec) {			
			// filter out read if either test fails.
			return (! mapQualityPredicate.test(rec) || !requiredTagPredicate.test(rec) || 
					! cellBarcodePredicate.test(rec) || ! editDistancePredicate.test (rec));
		}
										    
    }


	@Override
	public boolean hasNext() {
		return this.atoi.hasNext();
	}

	@Override
	public void close() {
		CloserUtil.close(this.atoi);
	}

}
