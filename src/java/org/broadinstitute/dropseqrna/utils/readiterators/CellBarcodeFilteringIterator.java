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
package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.PredicateFilteredIterator;

import java.util.*;

/**
 * Filter out all reads that don't have on of the listed cell barcode values.
 * @author nemesh
 *
 */
public class CellBarcodeFilteringIterator extends FilteredIterator<SAMRecord> {

	private final RequiredTagStringValuePredicate predicate;
	private final Log log = Log.getInstance(CellBarcodeFilteringIterator.class);

	/**
	 * Filters out reads that do not contain cell barcodes in the collection.
	 * If the collection has entries, then only reads with cell barcodes matching those entries will be retained.
	 * If the collection is null or empty, then this class will only filter out cell barcodes that match the missing value.
	 * The constructor that does not including the missing value sets it to "-" by default, which is the CellRanger / Starsolo default missing value.
	 * @param underlyingIterator An iterator of SAMRecords
	 * @param cellBarcodeTag A tag on the SAMRecord that encodes the cell barcode
	 * @param cellBarcodes A collection of cell barcodes, or a null/empty list to not filter cell barcodes
	 * @param missingValue A value of a cell barcode to indicate when that cell barcode is missing or failed cell barcode repair.  Used by CellRanger/Starsolo, defaults to "-".
	 */
	public CellBarcodeFilteringIterator (final Iterator<SAMRecord> underlyingIterator, final String cellBarcodeTag, final Collection<String> cellBarcodes, String missingValue) {
		super(underlyingIterator);
		// instantiate the predicate
		if (cellBarcodes==null || cellBarcodes.isEmpty()) {
			this.predicate = new RequiredTagStringValuePredicate(cellBarcodeTag, List.of(missingValue), true);
			String msg = String.format("No cell barcodes supplied to filter, explicitly filtering missing value [" + missingValue +"]");
			log.info(msg);
		}
		else
			this.predicate = new RequiredTagStringValuePredicate(cellBarcodeTag, cellBarcodes, false);
	}

	public CellBarcodeFilteringIterator (final Iterator<SAMRecord> underlyingIterator, final String cellBarcodeTag, final Collection<String> cellBarcodes) {
		this(underlyingIterator, cellBarcodeTag, cellBarcodes, "-");
	}

	/**
	 * @param rec The SAMRecord to test.
	 * @return true if record should be skipped
	 */
	@Override
	public boolean filterOut(SAMRecord rec) {
		// predicates return true when the test passes, but we want to filter out reads.
		return !predicate.test(rec);
	}

	@Override
	public void logFilterResults() {
		String msg = String.format("Cell barcode records pass [%d] records fail [%d] ",this.getRecordsPassed(), this.getRecordsFailed());
		log.info(msg);

	}

}
