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
package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.SAMRecord;

/**
 * Filter out all reads that don't have on of the listed cell barcode values.
 * @author nemesh
 *
 */
public class CellBarcodeFilteringIterator extends FilteredIterator<SAMRecord> {

	private final Set<String> cellBarcodes;
	private final String cellBarcodeTag;

	public CellBarcodeFilteringIterator (final Iterator<SAMRecord> underlyingIterator, final String cellBarcodeTag, final Collection<String> cellBarcodes) {
		super(underlyingIterator);
		this.cellBarcodeTag=cellBarcodeTag;
		if (cellBarcodes==null) this.cellBarcodes=null;
		else this.cellBarcodes=new HashSet<>(cellBarcodes);
	}

	@Override
	public boolean filterOut(final SAMRecord rec) {
		if (cellBarcodes==null) return false;
		String cellBarcode = rec.getStringAttribute(this.cellBarcodeTag);
		if (cellBarcode==null) return true;
		return (!this.cellBarcodes.contains(cellBarcode));
	}


}
