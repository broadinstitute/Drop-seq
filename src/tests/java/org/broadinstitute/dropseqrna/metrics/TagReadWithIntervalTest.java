package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.dropseqrna.utils.CompareBAMTagValues;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.Interval;
import picard.nio.PicardHtsPath;

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

		// compare the tagged BAM to the expected BAM.
		List<String> tags = Collections.singletonList("ZI");
		compareBAMTagValues(outBAM, TAGGED_BAM, tags, 0);

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

	private void compareBAMTagValues(File input1, File input2, List<String> tags, int expectedProgramValue) {
		CompareBAMTagValues cbtv = new CompareBAMTagValues();
		cbtv.INPUT_1 = Collections.singletonList(new PicardHtsPath(input1));
		cbtv.INPUT_2 = Collections.singletonList(new PicardHtsPath(input2));
		cbtv.TAGS_1 = tags;
		cbtv.TAGS_2 = tags;
		cbtv.STRICT = true;
		int result = cbtv.doWork();
		Assert.assertTrue(result == expectedProgramValue);
	}

}
