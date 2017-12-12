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
