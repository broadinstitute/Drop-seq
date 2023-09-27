package org.broadinstitute.dropseqrna.utils;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class BaseQualityFilterTest {

	@Test ()
	public void testScoreBaseQualityWrongLength() {
		SAMRecordSetBuilder b = new SAMRecordSetBuilder();
		// This read is plenty long and high quality.
		b.addFrag("read1", 0, 1, false, false, "20M", null, 30);
		// this read is too short.
		b.addFrag("read2", 0, 1, false, false, "5M", null, 30);
		Iterator<SAMRecord> iter = b.getRecords().iterator();
		SAMRecord r = iter.next();

		List <BaseRange> baseRangeList = new ArrayList<>();
		baseRangeList.add(new BaseRange(1, 12));

		BaseQualityFilter bqf = new BaseQualityFilter(baseRangeList, 20);
		int quality = bqf.scoreBaseQuality(r);
		Assert.assertEquals(quality, 0);

		// give a read that's too short.
		// this causes a null pointer exception, let's make it a nicer reported exception.
		r = iter.next();

		try {
			quality = bqf.scoreBaseQuality(r);
		} catch (TranscriptomeException e) {
			System.out.println(e.getMessage());
		}


		Assert.assertEquals(quality, 0);
	}
}
