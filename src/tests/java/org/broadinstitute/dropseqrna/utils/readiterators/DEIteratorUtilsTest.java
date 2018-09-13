package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

public class DEIteratorUtilsTest {
	@Test
	public void testRoundTrip() {
		List<String> tags = Arrays.asList("XC", "XM");
		List<Short> shortTags = DEIteratorUtils.getShortBAMTags(tags);
		List<String> fromShortToStringTags = DEIteratorUtils.getStringBAMTags(shortTags);
		Assert.assertEquals(tags, fromShortToStringTags);

	}
}
