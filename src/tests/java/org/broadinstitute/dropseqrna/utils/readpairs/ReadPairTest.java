package org.broadinstitute.dropseqrna.utils.readpairs;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.junit.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class ReadPairTest {
	@Test
	public void f() {
		List<SAMRecord> recs = getPairedRead("test", 1, 10, 20);

		Assert.assertNotNull(recs);
		ReadPair p = new ReadPair (recs.get(0), recs.get(1));
		SAMRecord r1 = p.getFirstRead();
		Assert.assertTrue(r1.getFirstOfPairFlag());
		SAMRecord r2 = p.getSecondRead();
		Assert.assertFalse(r2.getFirstOfPairFlag());

		SAMRecord left = p.getLeftRead();
		SAMRecord right = p.getRightRead();

		Assert.assertEquals(recs.get(0), left);
		Assert.assertEquals(recs.get(1), right);

		Assert.assertEquals(recs.get(0), p.getRead1());
		Assert.assertEquals(recs.get(1), p.getRead2());



	}

	@Test
	public void testFlipped () {
		List<SAMRecord> recs = getPairedRead("test", 1, 10, 20);

		Assert.assertNotNull(recs);
		// put in the reads the other way.  results still the same.
		ReadPair p = new ReadPair (recs.get(1), recs.get(0));
		SAMRecord r1 = p.getFirstRead();
		Assert.assertTrue(r1.getFirstOfPairFlag());
		SAMRecord r2 = p.getSecondRead();
		Assert.assertFalse(r2.getFirstOfPairFlag());

		SAMRecord left = p.getLeftRead();
		SAMRecord right = p.getRightRead();

		Assert.assertEquals(recs.get(0), left);
		Assert.assertEquals(recs.get(1), right);
		Assert.assertTrue(p.testProperlyPaired());

		Assert.assertNotNull(p.toString());
	}

	@Test
	public void testSet () {
		List<SAMRecord> recs = getPairedRead("test", 1, 10, 20);

		Assert.assertNotNull(recs);
		// put in the reads the other way.  results still the same.
		ReadPair p = new ReadPair();
		// set the reads in the wrong order.
		p.setRead1(recs.get(1));
		p.setRead2(recs.get(0));

		SAMRecord r1 = p.getFirstRead();
		SAMRecord r2 = p.getSecondRead();
		Assert.assertEquals(recs.get(0), r1);
		Assert.assertEquals(recs.get(1), r2);

		SAMRecord left = p.getLeftRead();
		SAMRecord right = p.getRightRead();

		Assert.assertEquals(recs.get(0), left);
		Assert.assertEquals(recs.get(1), right);



	}

	private List<SAMRecord> getPairedRead (final String name, final int contig, final int start1, final int start2) {
		List<SAMRecord> result = new ArrayList<> ();

		SAMRecordSetBuilder builder = new SAMRecordSetBuilder();
		builder.addPair(name, contig, start1, start2);

		Collection<SAMRecord> recs = builder.getRecords();

		for (SAMRecord r: recs) {
			if (r.getFirstOfPairFlag()) result.add(0, r);
			if (r.getSecondOfPairFlag()) result.add(1, r);
			r.setMappingQuality(10);
		}

		return (result);

	}

}
