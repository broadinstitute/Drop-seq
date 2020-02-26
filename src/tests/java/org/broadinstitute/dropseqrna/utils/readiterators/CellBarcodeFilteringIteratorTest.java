package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class CellBarcodeFilteringIteratorTest {

	@Test
	public void filterOut() {

		Iterator<SAMRecord> underlyingIterator = Collections.emptyIterator();
		Collection<String> cellBarcodesToKeep = Arrays.asList("CELL1", "CELL2");
		String cellBCTag ="XC";

		CellBarcodeFilteringIterator f = new CellBarcodeFilteringIterator(underlyingIterator, cellBCTag, cellBarcodesToKeep);
		SAMRecordSetBuilder b = new SAMRecordSetBuilder();
		// second argument is the contig index which is 0 based.  So contig index=0 -> chr1.  index=2 -> chr3, etc.
		b.addFrag("1", 0, 1, false);
		SAMRecord r1 = b.getRecords().iterator().next();

		r1.setAttribute(cellBCTag, "CELL1");
		Assert.assertFalse(f.filterOut(r1));

		r1.setAttribute(cellBCTag, "CELL2");
		Assert.assertFalse(f.filterOut(r1));

		r1.setAttribute(cellBCTag, "FOO");
		Assert.assertTrue(f.filterOut(r1));

	}
	
	@Test
	public void filterOut2() {
		Iterator<SAMRecord> underlyingIterator = Collections.emptyIterator();
		Collection<String> cellBarcodesToKeep = null;
		String cellBCTag ="XC";

		CellBarcodeFilteringIterator f = new CellBarcodeFilteringIterator(underlyingIterator, cellBCTag, cellBarcodesToKeep);
		SAMRecordSetBuilder b = new SAMRecordSetBuilder();
		// second argument is the contig index which is 0 based.  So contig index=0 -> chr1.  index=2 -> chr3, etc.
		b.addFrag("1", 0, 1, false);
		SAMRecord r1 = b.getRecords().iterator().next();

		r1.setAttribute(cellBCTag, "CELL1");
		Assert.assertFalse(f.filterOut(r1));

		r1.setAttribute(cellBCTag, "CELL2");
		Assert.assertFalse(f.filterOut(r1));

	}
	
	@Test
	public void filterOut3() {
		Iterator<SAMRecord> underlyingIterator = Collections.emptyIterator();
		Collection<String> cellBarcodesToKeep = Collections.EMPTY_LIST;
		String cellBCTag ="XC";

		CellBarcodeFilteringIterator f = new CellBarcodeFilteringIterator(underlyingIterator, cellBCTag, cellBarcodesToKeep);
		SAMRecordSetBuilder b = new SAMRecordSetBuilder();
		// second argument is the contig index which is 0 based.  So contig index=0 -> chr1.  index=2 -> chr3, etc.
		b.addFrag("1", 0, 1, false);
		SAMRecord r1 = b.getRecords().iterator().next();

		r1.setAttribute(cellBCTag, "CELL1");
		Assert.assertFalse(f.filterOut(r1));

		r1.setAttribute(cellBCTag, "CELL2");
		Assert.assertFalse(f.filterOut(r1));

	}
	
	
}
