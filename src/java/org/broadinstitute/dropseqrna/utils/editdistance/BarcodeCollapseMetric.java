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
package org.broadinstitute.dropseqrna.utils.editdistance;

import htsjdk.samtools.metrics.MetricBase;

public class BarcodeCollapseMetric extends MetricBase {

	public int editDistance;
	public int count;
	
	public BarcodeCollapseMetric() {
		
	}
	
	public BarcodeCollapseMetric(int editDistance, int count) {
		this.editDistance = editDistance;
		this.count = count;
	}
	
	public int getEditDistance() {
		return editDistance;
	}
	public void setEditDistance(int editDistance) {
		this.editDistance = editDistance;
	}
	public int getCount() {
		return count;
	}
	public void setCount(int count) {
		this.count = count;
	}
	
	public String toString () {
		StringBuilder b = new StringBuilder();
		b.append("Edit Distance: " + this.editDistance + " ");
		b.append("Num Barcodes: " + this.count);
		return (b.toString());
	}
	
	
	
}
