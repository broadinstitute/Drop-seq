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
	protected boolean filterOut(final SAMRecord rec) {
		if (cellBarcodes==null) return false;
		String cellBarcode = rec.getStringAttribute(this.cellBarcodeTag);
		if (cellBarcode==null) return true;
		return (!this.cellBarcodes.contains(cellBarcode));
	}


}
