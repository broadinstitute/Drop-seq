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
