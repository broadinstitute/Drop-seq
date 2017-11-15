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
