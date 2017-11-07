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
package org.broadinstitute.dropseqrna.metrics;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.util.Histogram;

public class ReadQualityMetrics extends MetricBase {
	
	/** the way the data is aggregated, if at all */
	public String aggregate="";
	
	/** The total number of reads */
	public long totalReads;
	
	/** The count of mapped reads */
	public long mappedReads;
	
	/** The count of high quality mapped reads - HQ is a score of 10 or more. */
	public long hqMappedReads;
	
	/** The number of high quality mapped reads that are not PCR duplicates */
	public long hqMappedReadsNoPCRDupes;
	
	/** The distribution of high quality mapped reads that are not PCR duplicates */
	private Histogram <Integer>histogram = null;
	
	private int mapQuality;
	
	/**
	 * @param mapQuality The map quality of a read to be considered high quality mapping.
	 * @param aggregate If the data should be aggregated at a tag level, this is the name of that aggregate level.
	 */
	public ReadQualityMetrics (int mapQuality, String aggregate, boolean gatherQualityHistogram) {
		this.mapQuality = mapQuality;
		this.aggregate=aggregate;
		
		if (gatherQualityHistogram) {
			histogram = new Histogram<Integer>("read quality", "num reads");
		}
	}

	/** No-arg ctor needed for instantiating with MetricsFile.read */
	public ReadQualityMetrics() {
	}

	public Histogram<Integer> getHistogram() {
		return histogram;
	}
	
	public void addRead (SAMRecord r) {
		// skip secondary of supplemental reads.
		if (r.isSecondaryOrSupplementary()) {
			return;
		}
		
		boolean isDupe = r.getDuplicateReadFlag();
		int mapQuality = r.getMappingQuality();
		boolean unmapped = r.getReadUnmappedFlag();
		if (histogram!=null) {
			histogram.increment(mapQuality);
		}
		
		totalReads++;
		if (!unmapped) {
			mappedReads++;
			if (mapQuality >= this.mapQuality) {
				hqMappedReads++;
				if (!isDupe) {
					hqMappedReadsNoPCRDupes++;					
				}
			}
		}
	}
	
	
}
