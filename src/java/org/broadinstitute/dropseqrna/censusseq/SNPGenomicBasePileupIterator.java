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
package org.broadinstitute.dropseqrna.censusseq;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.utils.CountChangingIteratorWrapper;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.ProgressLoggingIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.BamTagCountingIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.ChromosomeFilteringIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.MapQualityFilteredIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.PCRDuplicateFilteringIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.ProgressLogger;

/**
 * A genomic SNP pileup iterator serves up pileups of bases and qualities of reads on intervals, in the genomic order of the
 * intervals in the IntervalList provided.
 * @author nemesh
 *
 */
public class SNPGenomicBasePileupIterator implements CloseableIterator<SNPGenomicBasePileUp> {

	private static final Log log = Log.getInstance(SNPGenomicBasePileupIterator.class);
	private final GroupingIterator<SAMRecord> iter;
	private final String snpTag;
	private final String knownDonorTag;
	private final BamTagCountingIterator countingIter;
	private final Byte minBaseQuality;
	
	// private final OverlapDetector<Interval> snpOverlapDetector;
	public SNPGenomicBasePileupIterator(SamHeaderAndIterator headerAndIter, final IntervalList snpIntervals, final String snpTag,
            final int readMQ, final Collection<String> contigsToFilter, final String knownDonorTag, final Map<Interval, Double> meanGenotypeQuality) {
		this(headerAndIter, snpIntervals, snpTag, readMQ, contigsToFilter, knownDonorTag, meanGenotypeQuality, null);
	}
	
	// if meanGenotypeQuality is provided, reads that touch multiple SNPs will be filtered to the highest quality SNP per read.
	public SNPGenomicBasePileupIterator(final SamHeaderAndIterator headerAndIter, final IntervalList snpIntervals, final String snpTag,
            final int readMQ, final Collection<String> contigsToFilter, final String knownDonorTag, 
            final Map<Interval, Double> meanGenotypeQuality, final Integer minBaseQuality) {
		
		final ProgressLogger logger = new ProgressLogger(log);
		ProgressLoggingIterator loggingIterator = new ProgressLoggingIterator (headerAndIter.iterator, logger);
		OverlapDetector<Interval> snpOverlapDetector = new OverlapDetector<>(0, 0);
		snpOverlapDetector.addAll(snpIntervals.getIntervals(), snpIntervals.getIntervals());
		this.snpTag=snpTag;
		this.knownDonorTag = knownDonorTag;
		
		if (minBaseQuality==null)
			this.minBaseQuality=null; 
		else
			this.minBaseQuality=(byte) minBaseQuality.intValue();
		
		Iterator<SAMRecord> filteringIterator = new MapQualityFilteredIterator(loggingIterator, readMQ, true); 
		filteringIterator = new PCRDuplicateFilteringIterator(filteringIterator);
		// add chromosome filter.
		Collection<String> contigsToInclude = getContigs(snpIntervals);

		filteringIterator = new ChromosomeFilteringIterator(filteringIterator, contigsToInclude, false);
		// count the tags on the records flowing through.
		countingIter = new BamTagCountingIterator(filteringIterator, knownDonorTag);

		SNPTaggingIterator snpTaggingIterator = new SNPTaggingIterator(countingIter, snpOverlapDetector, snpTag, meanGenotypeQuality);

		SAMSequenceDictionary sd = snpIntervals.getHeader().getSequenceDictionary();

		IntervalTagComparator snpTagComparator = new IntervalTagComparator(this.snpTag, sd);

		final CloseableIterator<SAMRecord> sortingIterator =
                SamRecordSortingIteratorFactory.create(headerAndIter.header, snpTaggingIterator, snpTagComparator, logger);
		iter = new GroupingIterator<>(sortingIterator, snpTagComparator);

	}

	/**
	 * Get the list of contigs covered by the interval list of SNPs.
	 * @return
	 */
	public Collection<String> getContigs (final IntervalList snpIntervals) {
		Set<String> contigs = new HashSet<>();
		for (Interval i: snpIntervals.getIntervals())
			contigs.add(i.getContig());
		return contigs;
	}

	/**
	 * If a non-null knownDonorTag was passed into the constructor, this will return a non-null object counter
	 * that has counts of the number of reads per donor observed.
	 * @return
	 */
	public ObjectCounter<String> getKnownDonorCounts () {
		return this.countingIter.getCounts();
	}

	public String getKnownDonorTag () {
		return this.knownDonorTag;
	}

	@Override
	public boolean hasNext() {
		return this.iter.hasNext();
	}

	@Override
	/**
	 * Returns the next pileup.
	 * Pileups can contain 0 reads if all reads are filtered out.
	 * @return
	 */
	public SNPGenomicBasePileUp next() {
		if (!this.iter.hasNext())
			return null;

		Collection<SAMRecord> records = this.iter.next();
		// the first read and all other reads share the same tags, so can extract these properties out of the first read.
		SAMRecord firstRead = records.iterator().next();

		String snpID=firstRead.getStringAttribute(this.snpTag);
		Interval snpInterval = IntervalTagComparator.fromString(snpID);

		SNPGenomicBasePileUp p = new SNPGenomicBasePileUp(snpInterval);
		// add subsequent reads
		for (SAMRecord r: records)
			p.addRead(r);
		
		p.finalizePileup();
		if (this.minBaseQuality!=null)
			p=p.getFilteredPileupByBaseQuality(this.minBaseQuality);
		
		return p;
	}

	@Override
	public void close() {
		 CloserUtil.close(this.iter);
	}

	@Override
	public void remove () {
		this.iter.remove();
	}

	/**
	 * As reads flow through this iterator, emit the reads that overlap a SNP, tagged by the snpTag.
	 * @author nemesh
	 *
	 */
	private class SNPTaggingIterator extends CountChangingIteratorWrapper<SAMRecord> {

		private OverlapDetector<Interval> snpOverlapDetector;
		private final String snpTag;
		private final Map<Interval, Double> meanGenotypeQuality;
		
		protected SNPTaggingIterator(final Iterator<SAMRecord> underlyingIterator, final OverlapDetector<Interval> snpOverlapDetector, 
				final String snpTag, final Map<Interval, Double> meanGenotypeQuality) {
			super(underlyingIterator);
			this.snpOverlapDetector=snpOverlapDetector;
			this.snpTag = snpTag;
			this.meanGenotypeQuality=meanGenotypeQuality;
		}

		@Override
		// TODO: Refactor?  This is now a copy of SNPUMICellReadIteratorWrapper, but without the additional UMI features.
		protected void processRecord(final SAMRecord r) { 
			List<AlignmentBlock> blocks = r.getAlignmentBlocks();

			Collection<Interval> snpIntervals = new HashSet<>();
			
			for (AlignmentBlock b: blocks) {
				int start = b.getReferenceStart();
				int end = start + b.getLength() -1;

				Interval i = new Interval(r.getReferenceName(), start, end);
				Collection<Interval> overlaps = this.snpOverlapDetector.getOverlaps(i);
				for (Interval o: overlaps)
					snpIntervals.addAll(overlaps);
					// snps.add(IntervalTagComparator.toString(o));
			}
			
			// exit early if no SNPs found.
			if (snpIntervals.size()==0) return;
									
			// exit without cloning if not needed.
			if (snpIntervals.size()==1) {
				r.setAttribute(this.snpTag, IntervalTagComparator.toString(snpIntervals.iterator().next()));
				queueRecordForOutput(r);
				return;
			}
			
			if (this.meanGenotypeQuality!=null) {
				Interval bestSNP=getBestSNP(snpIntervals);				
				r.setAttribute(this.snpTag, IntervalTagComparator.toString(bestSNP));
				queueRecordForOutput(r);
				return;
			}
			
			// 1 read per SNP.  This is the old behavior where a read touching multiple SNPs treated each SNP observation as being independent.
			for (Interval snp:snpIntervals) {
				SAMRecord rr = Utils.getClone(r);
				rr.setAttribute(this.snpTag, IntervalTagComparator.toString(snp));
				queueRecordForOutput(rr);
			}
			
		}
		
		private Interval getBestSNP (Collection<Interval> snpIntervals) {
			// set the worse score to be worse than the missing value of -1.
			double maxGQ=-2d;
			Interval best=null;
			
			for (Interval i: snpIntervals) {
				double gq = this.meanGenotypeQuality.get(i);
				if (gq>maxGQ) {
					maxGQ=gq;
					best=i;
				}
			}
			return (best);
		}

	}

}
