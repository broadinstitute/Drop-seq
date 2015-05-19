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
