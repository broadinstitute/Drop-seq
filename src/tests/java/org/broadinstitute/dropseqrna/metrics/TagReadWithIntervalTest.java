package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.dropseqrna.utils.CompareBAMTagValues;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.Interval;

public class TagReadWithIntervalTest {
	// any BAM will do.
	private static final File IN_BAM = new File("testdata/org/broadinstitute/dropseq/metrics/NucBYReg4Reg.MOUSE.GCTAAGTAAGAT.Elp2.fixed.bam");
	private static final File TAGGED_BAM = new File("testdata/org/broadinstitute/dropseq/metrics/NucBYReg4Reg.MOUSE.GCTAAGTAAGAT.Elp2.tagged.bam");
	private static final File IN_INTERVAL = new File("testdata/org/broadinstitute/dropseq/metrics/NucBYReg4Reg.MOUSE.GCTAAGTAAGAT.Elp2.intervals");

	@Test
	public void testDoWork() {
		TagReadWithInterval t = new TagReadWithInterval();
		File outBAM=null;
		try {
			outBAM = File.createTempFile("TagReadWithIntervalTest.", ".bam");
			outBAM.deleteOnExit();
		} catch (IOException e) {
			e.printStackTrace();
		}

		t.INPUT=IN_BAM;
		t.INTERVALS=IN_INTERVAL;
		// t.OUTPUT=outBAM;
		t.OUTPUT=outBAM;
		int r = t.doWork();
		Assert.assertTrue (r==0);

		CompareBAMTagValues cbtv = new CompareBAMTagValues();
		cbtv.INPUT_1=outBAM;
		cbtv.INPUT_2=TAGGED_BAM;
		List<String> tags = new ArrayList<>();
		tags.add("ZI");
		cbtv.TAGS=tags;
		int result = cbtv.doWork();
		Assert.assertTrue(result==0);

	}

	@Test
	public void testGetIntervalName () {
		TagReadWithInterval t = new TagReadWithInterval();
		List<Interval> intervals = new ArrayList<>();
		intervals.add(new Interval ("1", 1,2,true, "foo"));
		intervals.add(new Interval ("1", 1,2,true, "bar"));
		String result = t.getIntervalName(intervals);
		Assert.assertEquals(result, "foo,bar");

	}
}
