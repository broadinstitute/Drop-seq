package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.util.Histogram;

/**
 * A structure to track how many reads have been accepted or rejected
 * @author nemesh
 *
 */
public class FilteredReadsMetric extends MetricBase {
	public long READS_REJECTED;
	public long READS_ACCEPTED;
	
	/** The distribution of TAG_VALUES accepted or rejected  */
	private Histogram <String> histogramTagsRejected = null;
	private Histogram <String> histogramTagsAccepted = null;
	
	
	public FilteredReadsMetric () {
		this.READS_ACCEPTED=0;
		this.READS_REJECTED=0;
		this.histogramTagsRejected = new Histogram<>("TAG_VALUE", "READS_REJECTED");
		this.histogramTagsAccepted = new Histogram<>("TAG_VALUE", "READS_ACCEPTED");
	}

	public FilteredReadsMetric copy () {
		return new FilteredReadsMetric();
	}
	
	
	public void incrementTagAccepted (final String tagValue) {
		if (tagValue!=null)
			histogramTagsAccepted.increment(tagValue);
	}
	
	public void incrementTagRejected (final String tagValue) {
		if (tagValue!=null)
			histogramTagsRejected.increment(tagValue);
	}
	
	public Histogram<String> getHistogramTagsRejected() {
		return histogramTagsRejected;
	}

	public Histogram<String> getHistogramTagsAccepted() {
		return histogramTagsAccepted;
	}

	public void merge(final FilteredReadsMetric other) {
		this.READS_ACCEPTED += other.READS_ACCEPTED;
		this.READS_REJECTED += other.READS_REJECTED;
		this.histogramTagsAccepted.addHistogram(other.histogramTagsAccepted);
		this.histogramTagsRejected.addHistogram(other.histogramTagsRejected);
	}
	
	
}
