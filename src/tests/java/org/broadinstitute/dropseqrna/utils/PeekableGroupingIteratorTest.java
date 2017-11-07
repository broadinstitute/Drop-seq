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
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.CloserUtil;

public class PeekableGroupingIteratorTest {


	@Test(enabled=true)
	public void testVsGrouping () {
		Comparator<String> cmp = String.CASE_INSENSITIVE_ORDER;

		String [] data = {"A", "A", "A", "B", "B", "C", "C", "C", "D"};

		GroupingIterator<String> gi = new GroupingIterator<>(Arrays.stream(data).iterator(), cmp);

		PeekableGroupingIterator<String> pi = new PeekableGroupingIterator<>(Arrays.stream(data).iterator(), cmp);

		while (gi.hasNext()){
			List<String> giResult = gi.next();
			List<String> piResult = new ArrayList<>();
			// load up the first result if available.
			if (pi.hasNext()) piResult.add(pi.next());
			// deal with subsequent entries in the group.
			while (pi.hasNextInGroup())
				piResult.add(pi.next());
			Assert.assertEquals(giResult.size(), piResult.size());
			Assert.assertEquals(giResult, piResult);
		}
		CloserUtil.close(gi);
		CloserUtil.close(pi);


	}
}
