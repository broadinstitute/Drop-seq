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
package org.broadinstitute.dropseqrna.readtrimming;

import org.broadinstitute.dropseqrna.readtrimming.SimplePolyAFinder;
import org.testng.Assert;
import org.testng.annotations.Test;


public class PolyAFinderTest {

	@Test(enabled=false, groups = { "dropseq", "transcriptome" })
	public void testTrimExample () {
        SimplePolyAFinder p = new SimplePolyAFinder(6, 2);
		String seq = "GGGTATCGCCGATTACAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACCCCCACCCCAATCAAAAACAAAAAAACAAAAAAAAACGCAAAAACACTCTCGGGCGGCCACACCCCACAACTACAACAAAA";
		int expected=14;
		int obs = p.getPolyAStart(seq).startPos;
		Assert.assertEquals(expected, obs);
	}

	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testTrimExample2 () {
        SimplePolyAFinder p = new SimplePolyAFinder(6, 1);
		String seq = "GTGAATGGGGGTATCGACGATTACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCGGAACCCCCCCCACCTGAAACAAAAACCTTAAACAAACAAACAAAAAAACCAAACAGCCCCCCCCCACAATAAAAACAAAAACACAGCACAC";
		int expected=22;
		int obs = p.getPolyAStart(seq).startPos;
		Assert.assertEquals(expected, obs);
	}

	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testTrimExample3 () {
        SimplePolyAFinder p = new SimplePolyAFinder(6, 0);
		String seq = "GAAAAGAAACTAAATATGAAAAGACTTCGAACTGACAATGTTTCAGACTTTTCTGAGAGCAGTGACTCAGAAAATTCAAATAAGAGAATAATAGATAATTCCTCAGAACAGAAGCCAGAGAATGAATTGAAAAAAAAAAAAAAAAAAAAA";
		int expected=129;
		int obs = p.getPolyAStart(seq).startPos;
		Assert.assertEquals(expected, obs);
	}

	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testTrimExample4 () {
        SimplePolyAFinder p = new SimplePolyAFinder(6, 0);
		String seq = "TGTAAAAAAAAAATGGGCCTAAAAAAAAAAAAAAAAAAAA";
		int expected=20;
		int obs = p.getPolyAStart(seq).startPos;
		Assert.assertEquals(expected, obs);
	}

	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testTrimExample5 () {
        SimplePolyAFinder p = new SimplePolyAFinder(6, 0);
		String seq = "GCGTGGAGGGGTAAAGCAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAA";
		int expected=38;
		int obs = p.getPolyAStart(seq).startPos;
		Assert.assertEquals(expected, obs);
	}

	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testTrimExample6 () {
        SimplePolyAFinder p = new SimplePolyAFinder(6, 1);
		String seq = "GCGTGGAGGGGTAAAGCAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAA";
		int expected=17;
		int obs = p.getPolyAStart(seq).startPos;
		Assert.assertEquals(expected, obs);
	}


	//

}
