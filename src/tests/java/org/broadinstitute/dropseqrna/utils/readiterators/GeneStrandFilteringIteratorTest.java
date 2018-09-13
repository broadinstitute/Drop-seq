package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Collections;
import java.util.Iterator;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class GeneStrandFilteringIteratorTest {

	@Test
	public void filterOut() {
		String strandTag = "GS";
		Iterator<SAMRecord> underlyingIterator = Collections.emptyIterator();
		GeneStrandFilteringIterator f = new GeneStrandFilteringIterator(underlyingIterator, strandTag);

		SAMRecordSetBuilder b = new SAMRecordSetBuilder();
		// positive strand read.
		b.addFrag("1", 1, 1, false);
		b.addFrag("1", 1, 1, true);
		Iterator<SAMRecord> recs = b.getRecords().iterator();

		SAMRecord r1 = recs.next();
		SAMRecord r2 = recs.next();

		// test all 4 possibilities.
		r1.setAttribute(strandTag, "+");
		boolean t1 = f.filterOut(r1);
		Assert.assertFalse(t1);

		r1.setAttribute(strandTag, "-");
		t1 = f.filterOut(r1);
		Assert.assertTrue(t1);


		r2.setAttribute(strandTag, "+");
		t1 = f.filterOut(r2);
		Assert.assertTrue(t1);

		r2.setAttribute(strandTag, "-");
		t1 = f.filterOut(r2);
		Assert.assertFalse(t1);

	}
}
