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
}
