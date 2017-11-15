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

import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.List;

import junit.framework.Assert;

import org.testng.annotations.Test;

public class TagBamWithReadSequenceExtendedTest {

	@Test (enabled=true)
	public void testHardClip () {
		SAMRecord r = new SAMRecord(null);
		r.setReadString("AAAAAAAAAATTTTTT"); //10A6T
		byte [] quals = {10,10,10,0,0,0,10,10,10,10,0,0,0,0,10,10}; // 0 qual at 4-6, 11-14
		r.setBaseQualities(quals);
		List<BaseRange> ranges = new ArrayList<BaseRange>();
		ranges.add(new BaseRange (4,6));
		ranges.add(new BaseRange (11,14));

		r=TagBamWithReadSequenceExtended.hardClipBasesFromRead(r, ranges);

		String newSeq = r.getReadString();
		byte [] newQuals = r.getBaseQualities();
		String expectedSeq ="AAAAAAATT";
		byte [] expectedQuals = {10,10,10,10,10,10,10,10,10};

		Assert.assertEquals(expectedSeq.length(), newSeq.length());
		Assert.assertEquals(expectedSeq, newSeq);
		for (int i=0; i<expectedQuals.length;i++)
			Assert.assertEquals(expectedQuals[i], newQuals[i]);
	}
}
