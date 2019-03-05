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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SummarizeUMIBaseQualities;
import org.testng.annotations.Test;

import htsjdk.samtools.util.StringUtil;
import junit.framework.Assert;

public class SummarizeUMIBaseQualitiesTest {

	@Test(enabled=true)
	public void test1() {

		char [] b = {'A', 'A', 'A'};
		byte [] bases = getBases(b);
		byte [] qualities = {10,20,30};

		SummarizeUMIBaseQualities s = new SummarizeUMIBaseQualities(bases, qualities);
		byte r = s.getMostCommonBase();
		Assert.assertEquals(r, (byte) 'A');

		int score = s.getSummarizedPhredScore();
		Assert.assertEquals(14, score);
	}

	@Test
	public void test2() {

		char [] b = {'A', 'A', 'A', 'A', 'A'};
		byte [] bases = getBases(b);
		byte [] qualities = {30,30,30,30,30};

		SummarizeUMIBaseQualities s = new SummarizeUMIBaseQualities(bases, qualities);
		byte r = s.getMostCommonBase();
		Assert.assertEquals(r, (byte) 'A');

		int score = s.getSummarizedPhredScore();
		Assert.assertEquals(30, score);
	}

	@Test
	public void test3() {

		char [] b = {'A', 'A', 'T'};
		byte [] bases = getBases(b);
		byte [] qualities = {30, 30, 10};

		SummarizeUMIBaseQualities s = new SummarizeUMIBaseQualities(bases, qualities);
		byte r = s.getMostCommonBase();
		Assert.assertEquals(r, (byte) 'A');

		int score = s.getSummarizedPhredScore();
		Assert.assertEquals(5, score);
	}

	@Test
	public void testMode () {
		char [] b = {'A', 'A', 'T', 'A', 'A', 'A', 'A'};
		byte [] bases = getBases(b);
		byte [] qualities = {20, 30, 10, 20, 25, 25, 25};

		SummarizeUMIBaseQualities s = new SummarizeUMIBaseQualities(bases, qualities);
		int value = s.getSummarizedPhreadScoreByMode();
		Assert.assertSame(25, value);
	}




	private byte [] getBases (final char [] bases) {
		byte [] result = new byte [bases.length];
		StringUtil.charsToBytes(bases, 0, bases.length, result, 0);
		return (result);
	}
}
