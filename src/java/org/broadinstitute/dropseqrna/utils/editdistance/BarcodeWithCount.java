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

import java.util.Comparator;

import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;


public class BarcodeWithCount {

	private String barcode;
	private int count;
	

	public BarcodeWithCount(String barcode, int count) {
		this.barcode = barcode;
		this.count = count;
	}

	public String getBarcode() {
		return barcode;
	}

	public int getCount() {
		return count;
	}
	
	public boolean equals(Object obj) {
		if (obj == null) { return false; }
		if (obj == this) { return true; }
		if (obj.getClass() != getClass()) {
			return false;
		}
		
		BarcodeWithCount rhs = (BarcodeWithCount) obj;
		return new EqualsBuilder()
			.appendSuper(super.equals(obj))
		    .append(barcode, rhs.barcode)
		    .isEquals();
		
	}
	
	public int hashCode() {
		return new HashCodeBuilder(17, 37).
	    append(barcode).toHashCode();
	}
	
	/**
     * Comparator to compare two barcodes with counts.
     * Sort them in descending order.
     */
    public static class CountComparator implements Comparator<BarcodeWithCount> {

        @Override
        public int compare(final BarcodeWithCount f1, final BarcodeWithCount f2) {
            final Integer v1 = new Integer(f1.getCount());
            final Integer v2 = new Integer(f2.getCount());
            return v1.compareTo(v2);
        }
    }
	

	
}
