package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.metrics.MetricBase;

/**
 * A structure to track how many reads have been accepted or rejected
 * @author nemesh
 *
 */
public class FilteredReadsMetric extends MetricBase {
	public long READS_REJECTED;
	public long READS_ACCEPTED;
	
	public FilteredReadsMetric () {
		this.READS_ACCEPTED=0;
		this.READS_REJECTED=0;
	}

	public void accumulate(final FilteredReadsMetric other) {
		this.READS_ACCEPTED += other.READS_ACCEPTED;
		this.READS_REJECTED += other.READS_REJECTED;
	}
}
