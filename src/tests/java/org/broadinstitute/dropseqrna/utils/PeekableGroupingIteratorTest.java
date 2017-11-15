/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
