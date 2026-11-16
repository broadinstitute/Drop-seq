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
package org.broadinstitute.dropseqrna.metagene;

import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorStrategy;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.utils.readiterators.MapQualityPredicate;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.broadinstitute.dropseqrna.utils.readiterators.umisort.UMIKeyCount;
import org.broadinstitute.dropseqrna.utils.readiterators.umisort.UMIKeyCountCodec;
import org.broadinstitute.dropseqrna.utils.sorting.CombiningSortingCollection;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.SortingCollection;
import picard.annotation.LocusFunction;

/**
 * Gathers a group of reads for a cell barcode and UMI barcode across all genes, without sorting uniquely mapped reads as SAMRecords.
 * <p>
 * A read that is uniquely mapped (map quality of at least mapQualityUnique) and resolves to one gene is a read group of one,
 * and only contributes its gene to the uniquely mapped genes of its cell and UMI.  Such reads are reduced to (cell, UMI, gene)
 * and sorted compactly, while all other informative reads are sorted as SAMRecords.  The two sorted streams are then joined
 * on (cell, UMI), in the manner of {@link org.broadinstitute.dropseqrna.censusseq.VCFPileupJointIterator}.
 * <p>
 * Uniquely mapped reads are not returned in the collections, so this cannot return uninformative reads to a sink.
 * It requires that a uniquely mapped read never shares its read name with another alignment, which holds for STAR.
 */
public class UMIMetaGeneJointIterator extends AbstractUMIMetaGeneCollectionIterator {

	/** Number of unique-read keys held in memory before a chunk is combined and spilled. */
	private static final int UNIQUE_KEY_MAX_RECORDS_IN_RAM = 2_000_000;

	/** The reads that are not uniquely mapped, grouped by cell and UMI. */
	private final PeekableIterator<List<SAMRecord>> multiMappedUmis;
	/** The genes of uniquely mapped reads, one group per cell and UMI. */
	private final PeekableIterator<List<UMIKeyCount>> uniqueUmis;
	private final boolean emitUniqueOnlyUmis;
	private UMIMetaGeneCollection pending=null;

	/**
	 * Construct an object that generates UMI objects from a BAM file for informative reads.
	 * See {@link AbstractUMIMetaGeneCollectionIterator} for the parameters and the definition of informative reads.
	 * All of the reads are consumed by the constructor.
	 *
	 * @param emitUniqueOnlyUmis If true, a UMI with only uniquely mapped reads is returned as a collection without reads
	 * that carries its uniquely mapped genes (needed to count unique genes).  If false, only UMIs with other reads are returned.
	 */
	public UMIMetaGeneJointIterator(final SamHeaderAndIterator headerAndIterator,
			final String geneTag,
			final String geneStrandTag,
			final String geneFunctionTag,
			final StrandStrategy strandStrategy,
			final Collection <LocusFunction> acceptedLociFunctions,
			final FunctionalDataProcessorStrategy functionStrategy,
			final String cellBarcodeTag,
			final String molecularBarcodeTag, final int minMapQuality, final int mapQualityUnique, final Integer maxEditDistance,
			final List<String> cellBarcodes, final boolean assignReadsToAllGenes, final boolean emitUniqueOnlyUmis) {
		super(headerAndIterator, geneTag, geneStrandTag, geneFunctionTag, strandStrategy, acceptedLociFunctions, functionStrategy, cellBarcodeTag,
				molecularBarcodeTag, minMapQuality, mapQualityUnique, maxEditDistance, cellBarcodes, null, assignReadsToAllGenes);
		this.emitUniqueOnlyUmis = emitUniqueOnlyUmis;

		final SortingCollection<SAMRecord> multiMappedSort = SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(headerAndIterator.header),
				this.cellAndUmiComparator, SAMFileWriterImpl.getDefaultMaxRecordsInRam());
		final CombiningSortingCollection<UMIKeyCount> uniqueSort = new CombiningSortingCollection<>(new UMIKeyCountCodec(),
				UMIKeyCount.cellUmiSortOrder(), UMIKeyCount::combine, UNIQUE_KEY_MAX_RECORDS_IN_RAM);
		splitReads(multiMappedSort, uniqueSort);
		this.multiMappedUmis = new PeekableIterator<>(new GroupingIterator<>(multiMappedSort.iterator(), this.cellAndUmiComparator));
		this.uniqueUmis = new PeekableIterator<>(new GroupingIterator<>(uniqueSort.iterator(), UMIKeyCount.cellUmiGroupOrder()));
	}

	/**
	 * Sends each informative read to one of the sorts.
	 * @param multiMappedSort receives the reads that are not uniquely mapped, as records
	 * @param uniqueSort receives (cell, UMI, gene) for the uniquely mapped reads
	 */
	private void splitReads (final SortingCollection<SAMRecord> multiMappedSort, final CombiningSortingCollection<UMIKeyCount> uniqueSort) {
		final MapQualityPredicate uniquePredicate = new MapQualityPredicate(this.mapQualityUnique, false);
		final StringInterner geneCache = new StringInterner();
		while (this.informativeReads.hasNext()) {
			final SAMRecord r = this.informativeReads.next();
			this.progress.record(r);
			if (uniquePredicate.test(r)) {
				// not informative if the read does not resolve to exactly one gene, as in processRead.
				final FunctionalData fd = this.gfteratorWrapper.getSingleFunction(r);
				if (fd!=null)
					uniqueSort.add(new UMIKeyCount(r.getStringAttribute(this.cellBarcodeTag), geneCache.intern(fd.getGene()),
							r.getStringAttribute(this.molecularBarcodeTag), 1));
			} else
				multiMappedSort.add(r);
		}
		CloserUtil.close(this.informativeReads);
	}

	@Override
	public boolean hasNext() {
		if (this.pending==null)
			this.pending = nextJoined();
		return this.pending!=null;
	}

	@Override
	public UMIMetaGeneCollection next () {
		if (!hasNext()) throw new NoSuchElementException();
		final UMIMetaGeneCollection result = this.pending;
		this.pending = null;
		return result;
	}

	/**
	 * Join the two sorted streams on (cell, UMI).  Advances to the next (cell, UMI) in either stream that gives a result.
	 * @return the next collection, or null if the streams are exhausted.
	 */
	private UMIMetaGeneCollection nextJoined () {
		while (this.multiMappedUmis.hasNext() || this.uniqueUmis.hasNext()) {
			final int cmp;
			if (!this.uniqueUmis.hasNext()) cmp = -1;
			else if (!this.multiMappedUmis.hasNext()) cmp = 1;
			else cmp = compareCellAndUmi(this.multiMappedUmis.peek().getFirst(), this.uniqueUmis.peek().getFirst());

			List<SAMRecord> reads = Collections.emptyList();
			List<UMIKeyCount> uniqueKeys = Collections.emptyList();
			if (cmp<=0) reads = this.multiMappedUmis.next();
			if (cmp>=0) uniqueKeys = this.uniqueUmis.next();

			final Set<String> uniqueGenes = new HashSet<>();
			for (UMIKeyCount k: uniqueKeys) uniqueGenes.add(k.gene);

			final List<SAMRecord> informativeReads = getRecordsPassingFilters(reads);
			final UMIMetaGeneCollection c;
			if (!informativeReads.isEmpty())
				c = UMIMetaGeneCollection.withUniqueGenes(this.geneTag, this.molecularBarcodeTag, this.cellBarcodeTag, informativeReads, this.mapQualityUnique, uniqueGenes);
			else if (this.emitUniqueOnlyUmis && !uniqueKeys.isEmpty()) {
				// A UMI with only uniquely mapped reads has nothing to tag, but counts towards the unique genes of discovery.
				final UMIKeyCount first = uniqueKeys.getFirst();
				c = UMIMetaGeneCollection.uniqueGenesOnly(this.geneTag, this.molecularBarcodeTag, this.cellBarcodeTag, this.mapQualityUnique, first.cell, first.umi, uniqueGenes);
			} else
				continue;
			c.populateGeneSets();
			return c;
		}
		return null;
	}

	private int compareCellAndUmi (final SAMRecord read, final UMIKeyCount key) {
		final int cmp = read.getStringAttribute(this.cellBarcodeTag).compareTo(key.cell);
		if (cmp!=0) return cmp;
		return read.getStringAttribute(this.molecularBarcodeTag).compareTo(key.umi);
	}

	@Override
	public void close() {
		CloserUtil.close(this.multiMappedUmis);
		CloserUtil.close(this.uniqueUmis);
	}

}
