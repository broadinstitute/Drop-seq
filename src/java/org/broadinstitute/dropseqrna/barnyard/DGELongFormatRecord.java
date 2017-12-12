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
package org.broadinstitute.dropseqrna.barnyard;

import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

import com.google.common.collect.ComparisonChain;

/**
 * A simple structure to hold the cell/gene/UMI count.
 * @author nemesh
 *
 */
public class DGELongFormatRecord {

	private final String cell;
	private final String gene;
	private final int count;

	public DGELongFormatRecord (final String cell, final String gene, final int count) {
		this.cell=cell;
		this.gene=gene;
		this.count=count;
	}

	public String getCell() {
		return cell;
	}

	public String getGene() {
		return gene;
	}

	public int getCount() {
		return count;
	}

	/**
	 * Orders DGELongFormatRecords first by cell barcodes (in the order specified in the list), then by UMI count of genes, then finally by gene name.
	 * @author nemesh
	 *
	 */
	public static class CellBarcodeOrderComparator implements Comparator<DGELongFormatRecord>{

		private final Map<String, Integer> cellBarcodeOrder;

		public CellBarcodeOrderComparator (final Collection<String> cellBarcodes) {
			cellBarcodeOrder = new HashMap<>(cellBarcodes.size());
			int counter=0;
			for (String string : cellBarcodes) {
				cellBarcodeOrder.put(string, counter);
				counter++;
			}
		}

		@Override
		public int compare(final DGELongFormatRecord o1, final DGELongFormatRecord o2) {
			Integer pos1=cellBarcodeOrder.get(o1.getCell());
			Integer pos2=cellBarcodeOrder.get(o2.getCell());
			if (pos1==null || pos2==null)
				throw new IllegalStateException ("Trying to compare cell barcode records where the ordering is not known.  Please have all cell barcodes in input list.");
			int result = ComparisonChain.start().compare(pos1, pos2).compare(o2.getCount(), o1.getCount()).compare(o1.getGene(), o2.getGene()).result();
			return (result);
		}

	}



}
