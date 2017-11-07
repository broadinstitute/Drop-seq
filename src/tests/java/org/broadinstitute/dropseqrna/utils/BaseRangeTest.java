/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils;

import java.util.ArrayList;
import java.util.List;

import junit.framework.Assert;

import org.testng.annotations.Test;

public class BaseRangeTest {

	@Test(enabled=true)
	public void parseBaseRange() {
		String rangeStr = "13-20";
		char [] foo = rangeStr.toCharArray();
		int val = foo[2];
		int val2 = foo[1];
		List<BaseRange> baseRanges = BaseRange.parseBaseRange(rangeStr);
		Assert.assertNotNull(baseRanges);
	}

	@Test(enabled=true)
	public void parseBaseRange3() {
		String rangeStr = "13 -20";

		List<BaseRange> baseRanges = BaseRange.parseBaseRange(rangeStr);
		Assert.assertNotNull(baseRanges);
	}

	@Test
	public void parseBaseRange2() {
		// weird error where ascii character 173 is in here.
		String BASE_RANGE="13\u00ad-20";
		char [] foo = BASE_RANGE.toCharArray();
		int val = foo[2];
		int val2 = foo[1];

		List<BaseRange> baseRanges = BaseRange.parseBaseRange(BASE_RANGE);
		Assert.assertNotNull(baseRanges);

	}

	@Test(expectedExceptions=IllegalArgumentException.class)
	public void parseBaseRangeBAD() {
		// weird error where ascii character 173 is in here.
		String BASE_RANGE="13\u00ad20";
		char [] foo = BASE_RANGE.toCharArray();
		int val = foo[2];
		int val2 = foo[1];

		List<BaseRange> baseRanges = BaseRange.parseBaseRange(BASE_RANGE);
		Assert.assertNotNull(baseRanges);

	}

	@Test(enabled=true, groups = { "dropseq", "transcriptome"})
	public void testParseBaseRange() {
		String test = "1-4:17-22";
		List<BaseRange> ranges = BaseRange.parseBaseRange(test);
		BaseRange b=ranges.get(0);
		Assert.assertEquals(b.getStart(), 1);
		Assert.assertEquals(b.getEnd(), 4);

		BaseRange b2=ranges.get(1);
		Assert.assertEquals(b2.getStart(), 17);
		Assert.assertEquals(b2.getEnd(), 22);

	}

	@Test(enabled=true, groups = { "dropseq", "transcriptome"})
	public void testParseBaseRange2() {
		String test = "17-22";
		List<BaseRange> ranges = BaseRange.parseBaseRange(test);

		BaseRange b2=ranges.get(0);
		Assert.assertEquals(b2.getStart(), 17);
		Assert.assertEquals(b2.getEnd(), 22);

	}

	@Test(enabled=true, groups = { "dropseq","transcriptome"})
	public void testParseSingleBaseRange() {
		String test = "1-4";
		BaseRange b = BaseRange.parseSingleBaseRange(test);
		Assert.assertEquals(b.getStart(), 1);
		Assert.assertEquals(b.getEnd(), 4);
	}

	@Test(enabled=false, groups = { "dropseq","transcriptome" })
	public void testParseSingleBaseRange2() {
		String test = "Foo";
		BaseRange b = BaseRange.parseSingleBaseRange(test);
		Assert.assertEquals(b.getStart(), 1);
		Assert.assertEquals(b.getEnd(), 4);
	}

	@Test(enabled=true, groups = { "dropseq","transcriptome" })
	public void testGetTotalRangeSize() {
		String test = "1-4:17-22";
		int t1=10;
		int t2=BaseRange.getTotalRangeSize(test);
		Assert.assertEquals(t1, t2);

	}

	@Test(enabled=true, groups = { "dropseq","transcriptome" })
	public void testGetSequenceForBaseRange() {
		String seq="NNNNJJJJJJJJJJJJNNNNVV";
		String test = "1-4:17-22";
		List<BaseRange> ranges = BaseRange.parseBaseRange(test);
		String seq1 = BaseRange.getSequenceForBaseRange (ranges, seq);
		String seqTrue="NNNNNNNNVV";
		Assert.assertEquals(seq1, seqTrue);


		String test2 = "5-16";
		List<BaseRange> ranges2 = BaseRange.parseBaseRange(test2);
		String seq2 = BaseRange.getSequenceForBaseRange (ranges2, seq);
		String seq2True="JJJJJJJJJJJJ";
		Assert.assertEquals(seq2, seq2True);

	}

	@Test
	public void testSize () {
		List<BaseRange> list = new ArrayList<BaseRange>();
		list.add(new BaseRange(1,5));
		list.add(new BaseRange(11,15));
		int size = BaseRange.getTotalRangeSize(list);
		Assert.assertEquals(10, size);
	}

	@Test
	public void testInvert () {
		List<BaseRange> list = new ArrayList<BaseRange>();
		list.add(new BaseRange(1,5));
		list.add(new BaseRange(11,16));

		List<BaseRange> inverted = BaseRange.invert(list, 20);
		for (int i=0; i<inverted.size(); i++) {
			BaseRange br = inverted.get(i);
			if (i==0) {
				Assert.assertEquals(6, br.getStart());
				Assert.assertEquals(10, br.getEnd());
			}
			if (i==1) {
				Assert.assertEquals(17, br.getStart());
				Assert.assertEquals(20, br.getEnd());
			}
		}
	}

	@Test
	public void testCopySequenceRange () {
		byte [] data = {0,0,0,1,1,1,2,2,2};
		List<BaseRange> list = new ArrayList<BaseRange>();
		list.add(new BaseRange(1,3));
		list.add(new BaseRange(7,9));
		byte [] out = BaseRange.getBytesForBaseRange(list, data);
		byte [] expected = {0,0,0,2,2,2};
		for (int i=0; i<expected.length; i++)
			Assert.assertEquals(expected[i], out[i]);

	}


}
