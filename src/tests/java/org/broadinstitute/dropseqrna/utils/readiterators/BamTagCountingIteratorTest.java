package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Collections;
import java.util.Iterator;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class BamTagCountingIteratorTest {
	@Test
	public void testFilter() {
		Iterator<SAMRecord> underlyingIterator = Collections.emptyIterator();
		String tagName = "XC";
		// filter always returns false.
		BamTagCountingIterator f = new BamTagCountingIterator(underlyingIterator, tagName);


		SAMRecordSetBuilder b = new SAMRecordSetBuilder();
		// second argument is the contig index which is 0 based.  So contig index=0 -> chr1.  index=2 -> chr3, etc.
		b.addFrag("1", 0, 1, false);
		SAMRecord r1 = b.getRecords().iterator().next();
		r1.setAttribute(tagName, "A");
		// this one counts.
		Assert.assertFalse(f.filterOut(r1));

		// wrong tag set, correct tag is not set.
		r1.setAttribute("XZ", "A");
		r1.setAttribute(tagName, null);
		Assert.assertFalse(f.filterOut(r1));

		r1.setAttribute(tagName, "A");
		// this one counts.
		Assert.assertFalse(f.filterOut(r1));

		r1.setAttribute(tagName, "B");
		// this one counts.
		Assert.assertFalse(f.filterOut(r1));

		r1.setAttribute(tagName, "C");
		// this one counts.
		Assert.assertFalse(f.filterOut(r1));

		ObjectCounter<String> expected = f.getCounts();
		expected.incrementByCount("A", 2);
		expected.incrementByCount("B", 1);
		expected.incrementByCount("C", 1);
		ObjectCounter<String> actual = f.getCounts();

		Assert.assertEquals(actual, expected);



	}
}
