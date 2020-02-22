package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;

import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import junit.framework.Assert;

public class ChromosomeFilteringIteratorTest {
	@Test
	public void testFilterExclude() {
		Iterator<SAMRecord> underlyingIterator = Collections.emptyIterator();
		Collection<String> contigsToFilter = Arrays.asList("chr1", "chr2");

		ChromosomeFilteringIterator f = new ChromosomeFilteringIterator(underlyingIterator, contigsToFilter);

		SAMRecordSetBuilder b = new SAMRecordSetBuilder();
		// second argument is the contig index which is 0 based.  So contig index=0 -> chr1.  index=2 -> chr3, etc.
		b.addFrag("chr1", 0, 1, false);
		b.addFrag("chr2", 1, 1, false);
		b.addFrag("chr3", 2, 1, false);
		b.addFrag("chr4", 3, 1, false);
		Iterator<SAMRecord> recs = b.getRecords().iterator();

		Assert.assertTrue(f.filterOut(recs.next()));
		Assert.assertTrue(f.filterOut(recs.next()));
		// 3 and 4 not on the list.
		Assert.assertFalse(f.filterOut(recs.next()));
		Assert.assertFalse(f.filterOut(recs.next()));


	}

	@Test
	public void testFilterInclude() {
		Iterator<SAMRecord> underlyingIterator = Collections.emptyIterator();
		Collection<String> contigsToFilter = Arrays.asList("chr1", "chr2");

		ChromosomeFilteringIterator f = new ChromosomeFilteringIterator(underlyingIterator, contigsToFilter, false);

		SAMRecordSetBuilder b = new SAMRecordSetBuilder();
		// second argument is the contig index which is 0 based.  So contig index=0 -> chr1.  index=2 -> chr3, etc.
		b.addFrag("chr1", 0, 1, false);
		b.addFrag("chr2", 1, 1, false);
		b.addFrag("chr3", 2, 1, false);
		b.addFrag("chr4", 3, 1, false);
		Iterator<SAMRecord> recs = b.getRecords().iterator();

		Assert.assertFalse(f.filterOut(recs.next()));
		Assert.assertFalse(f.filterOut(recs.next()));
		// 3 and 4 not on the list.  Reject.
		Assert.assertTrue(f.filterOut(recs.next()));
		Assert.assertTrue(f.filterOut(recs.next()));


	}


	@Test
	public void testFilterEmpty() {
		Iterator<SAMRecord> underlyingIterator = Collections.emptyIterator();
		Collection<String> contigsToFilter = null;

		ChromosomeFilteringIterator f = new ChromosomeFilteringIterator(underlyingIterator, contigsToFilter, true);

		SAMRecordSetBuilder b = new SAMRecordSetBuilder();
		// second argument is the contig index which is 0 based.  So contig index=0 -> chr1.  index=2 -> chr3, etc.
		b.addFrag("1", 0, 1, false);
		b.addFrag("2", 1, 1, false);
		b.addFrag("3", 2, 1, false);
		b.addFrag("4", 3, 1, false);
		Iterator<SAMRecord> recs = b.getRecords().iterator();

		Assert.assertFalse(f.filterOut(recs.next()));
		Assert.assertFalse(f.filterOut(recs.next()));
		// 3 and 4 not on the list.
		Assert.assertFalse(f.filterOut(recs.next()));
		Assert.assertFalse(f.filterOut(recs.next()));


	}
}
