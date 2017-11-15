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
