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
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorStrategy;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.SamWriterSink;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloserUtil;
import picard.annotation.LocusFunction;

/**
 * Gathers a group of reads for a cell barcode and UMI barcode across all genes.
 * All informative reads are sorted as SAMRecords by cell and UMI, which allows uninformative reads to be returned to a sink.
 * See {@link UMIMetaGeneJointIterator} for a faster iterator when uninformative reads are not needed.
 * @author nemesh
 *
 */
public class UMIMetaGeneCollectionIterator extends AbstractUMIMetaGeneCollectionIterator {

	private final GroupingIterator<SAMRecord> atoi;
	private final SamWriterSink sink;

	/**
	 * Construct an object that generates UMI objects from a BAM file for informative reads.
	 * See {@link AbstractUMIMetaGeneCollectionIterator} for the parameters and the definition of informative reads.
	 * If a sink is passed in, uninformative reads are added to the sink.
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
		super(headerAndIterator, geneTag, geneStrandTag, geneFunctionTag, strandStrategy, acceptedLociFunctions, functionStrategy, cellBarcodeTag,
				molecularBarcodeTag, minMapQuality, mapQualityUnique, maxEditDistance, cellBarcodes, sink, assignReadsToAllGenes);
		this.sink=sink;
		this.atoi = new GroupingIterator<>(SamRecordSortingIteratorFactory.create(headerAndIterator.header, this.informativeReads, this.cellAndUmiComparator, this.progress),
				this.cellAndUmiComparator);
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
			// Records are tracked by identity.  Hashing a SAMRecord hashes all of its attributes, which was the largest
			// calling-thread cost of this loop.
			List<SAMRecord> allReads = this.atoi.next();
			List<SAMRecord> informativeReads = getRecordsPassingFilters(allReads);

			// allReads - informative reads = uninformative reads, which go to the sink if there is one.
			if (sink!=null) {
				Set<SAMRecord> informativeSet = Collections.newSetFromMap(new IdentityHashMap<>());
				informativeSet.addAll(informativeReads);
				for (SAMRecord r: allReads)
					if (!informativeSet.contains(r)) sink.add(r);
			}
			// informative reads found!
			if (!informativeReads.isEmpty()) {
				UMIMetaGeneCollection c = new UMIMetaGeneCollection(this.geneTag, this.molecularBarcodeTag, this.cellBarcodeTag, informativeReads, this.mapQualityUnique);
				c.populateGeneSets();
				return (c);
			}
		}

		// no reads left, return null.
		return null;
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
