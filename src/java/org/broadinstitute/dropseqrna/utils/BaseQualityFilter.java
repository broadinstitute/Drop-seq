/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils;

import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.dropseqrna.TranscriptomeException;

import htsjdk.samtools.SAMRecord;

/**
 * Finds if specified potions of a read are below a minimum base quality.
 * @author nemesh
 *
 */
public class BaseQualityFilter {

	private final List<BaseRange> baseRanges;
	private final int baseQualityThrehsold;
	private FailedBaseMetric metric = null;

	public BaseQualityFilter (final List<BaseRange> baseRanges, final int baseQualityThrehsold) {
		this.baseQualityThrehsold=baseQualityThrehsold;
		this.baseRanges=baseRanges;
		this.metric = new FailedBaseMetric(BaseRange.getTotalRangeSize(baseRanges));
	}

	public int scoreBaseQuality(final SAMRecord barcodedRead) {
		int numBasesBelowQuality=0;
		byte [] qual= barcodedRead.getBaseQualities();
		char [] seq = barcodedRead.getReadString().toUpperCase().toCharArray();

		for (BaseRange b: baseRanges)
			for (int i=b.getStart()-1; i<b.getEnd(); i++) {
				if (i> (seq.length-1))
					throw new TranscriptomeException("Base [" + Integer.toString(i+1) + "] was requested, but the read isn't long enough ["+ barcodedRead.getReadString()+"]");

				byte q = qual[i];
				char s = seq[i];

				if (q < this.baseQualityThrehsold || s=='N')
					numBasesBelowQuality++;
			}

		this.metric.addFailedBase(numBasesBelowQuality);
		return (numBasesBelowQuality);
	}

	public FailedBaseMetric getMetric () {
		return this.metric;
	}

	public List<BaseRange> getBaseRanges () {
		return this.baseRanges;
	}

	public int getBaseQualityThreshold () {
		return this.baseQualityThrehsold;
	}

	public static class FailedBaseMetric {
		List<Integer> data = null;

		public FailedBaseMetric (final Integer length){
			data=new ArrayList<>(length+1);
			for (int i=0; i<=length; i++)
				data.add(Integer.valueOf(0));
		}

		public void addFailedBase(final int numBasesFailed) {
			addMultipleFailedBase(numBasesFailed, 1);
		}

		/**
		 * Increment numReads for the given number of failed bases by an arbitrary number
		 * @param numBasesFailed
		 * @param numReads
		 */
		public void addMultipleFailedBase(final int numBasesFailed, final int numReads) {
			while (getLength() <= numBasesFailed) {
				data.add(0);
			}
			data.set(numBasesFailed, data.get(numBasesFailed) + numReads);
		}

		public int getNumFailedBases(final int position) {
			return (data.get(position));
		}

		public int getLength() {
			return (data.size());
		}

		public void merge(final FailedBaseMetric other) {
			while (getLength() < other.getLength()) {
				data.add(0);
			}
			for (int i = 0; i < other.getLength(); ++i) {
				addMultipleFailedBase(i, other.getNumFailedBases(i));
			}
		}

	}


}
