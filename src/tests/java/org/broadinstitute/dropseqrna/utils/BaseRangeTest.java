package org.broadinstitute.dropseqrna.utils;

import java.util.List;

import junit.framework.Assert;

import org.testng.annotations.Test;

public class BaseRangeTest {

	@Test(enabled=true)
	public void parseBaseRange() {
		String rangeStr = "13-20";
		char [] foo = rangeStr.toCharArray();
		int val = (int) foo[2];
		int val2 = (int) foo[1];
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
		int val = (int) foo[2];
		int val2 = (int) foo[1];
		
		List<BaseRange> baseRanges = BaseRange.parseBaseRange(BASE_RANGE);
		Assert.assertNotNull(baseRanges);
		
	}
}
