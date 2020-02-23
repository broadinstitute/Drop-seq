package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Collections;
import java.util.Iterator;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class EditDistanceFilteringIteratorTest {
	@Test
	public void f() {
		Iterator<SAMRecord> underlyingIterator = Collections.emptyIterator();
		EditDistanceFilteringIterator f = new EditDistanceFilteringIterator(underlyingIterator, 2);

		SAMRecordSetBuilder b = new SAMRecordSetBuilder();
		// second argument is the contig index which is 0 based.  So contig index=0 -> chr1.  index=2 -> chr3, etc.
		b.addFrag("1", 0, 1, false);
		SAMRecord r = b.getRecords().iterator().next();

		r.setAttribute(ReadEditDistancePredicate.EDIT_DISTANCE_TAG, 1);
		Assert.assertFalse(f.filterOut(r));
		r.setAttribute(ReadEditDistancePredicate.EDIT_DISTANCE_TAG, 2);
		Assert.assertFalse(f.filterOut(r));
		r.setAttribute(ReadEditDistancePredicate.EDIT_DISTANCE_TAG, 3);
		Assert.assertTrue(f.filterOut(r));
	}
}
