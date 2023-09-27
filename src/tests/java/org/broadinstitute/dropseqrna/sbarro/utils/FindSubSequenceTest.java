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
package org.broadinstitute.dropseqrna.sbarro.utils;

import htsjdk.samtools.util.SequenceUtil;
import org.testng.Assert;

import org.broadinstitute.dropseqrna.sbarro.utils.FindSubSequence;
import org.broadinstitute.dropseqrna.sbarro.utils.SubSequenceResultI;
import org.testng.annotations.Test;

public class FindSubSequenceTest {

	private String anchorSequence = "CCGGTGGCGCCACTGC";


	/**
	 * GLOBAL ALIGNMENT TESTS
	 */
	@Test(enabled=true)
	public void findSequenceGlobal() {
		String sequence = "AAAGGGCCCTTTCCGGTGGCGCCACTGCTTTGGGCCCAAA";
		FindSubSequence f = new FindSubSequence(anchorSequence);

		SubSequenceResultI r = f.findSequenceGlobalAlignment(sequence);
		int ed = r.getEditDistance().getEditDistance();
		Assert.assertEquals(0, ed);

		Assert.assertEquals(anchorSequence, r.getQuerySequence());
		Assert.assertEquals(sequence, r.getTargetSequence());
		Assert.assertEquals("CCGGTGGCGCCACTGC", r.getSubSequence());

		int start = 13;
		int end = 28;
		int matchLen=16;

		Assert.assertEquals(start, r.getStart());
		Assert.assertEquals(end, r.getEnd());
		Assert.assertEquals(matchLen, r.getMatchLength());
	}

	@Test(enabled=true)
	// we don't expect this to work yet.
	public void findSequenceReverseComplimentedGlobal() {
		String sequence = "AAAGGGCCCTTTCCGGTGGCGCCACTGCTTTGGGCCCAAA";
		sequence = SequenceUtil.reverseComplement(sequence);
		FindSubSequence f = new FindSubSequence(anchorSequence);

		SubSequenceResultI r = f.findSequenceGlobalAlignment(sequence);
		int ed = r.getEditDistance().getEditDistance();
		Assert.assertEquals(14, ed);

		Assert.assertEquals(anchorSequence, r.getQuerySequence());
		Assert.assertEquals(sequence, r.getTargetSequence());
		Assert.assertEquals("CCCAAAGCAGTGGCGCCACCGGAAAGGGC", r.getSubSequence());

		int start = 7;
		int end = 35;
		int matchLen=29;

		Assert.assertEquals(start, r.getStart());
		Assert.assertEquals(end, r.getEnd());
		Assert.assertEquals(matchLen, r.getMatchLength());
	}


	@Test(enabled=true)
	// has a 1 BP insertion in the anchor sequence.
	public void findSequenceWithInsertionGlobal() {
		String sequence = "AAAGGGCCCTTTCCGGTGGCCGCCACTGCTTTGGGCCCAAA";
		FindSubSequence f = new FindSubSequence(anchorSequence);

		SubSequenceResultI r = f.findSequenceGlobalAlignment(sequence);
		int ed = r.getEditDistance().getEditDistance();
		Assert.assertEquals(1, ed);

		Assert.assertEquals(anchorSequence, r.getQuerySequence());
		Assert.assertEquals(sequence, r.getTargetSequence());
		Assert.assertEquals("CCGGTGGCCGCCACTGC", r.getSubSequence());

		int start = 13;
		int end = 29;
		int matchLen=17;

		Assert.assertEquals(start, r.getStart());
		Assert.assertEquals(end, r.getEnd());
		Assert.assertEquals(matchLen, r.getMatchLength());
	}

	@Test(enabled=true)
	// has a 1 BP insertion in the anchor sequence.
	public void findSequenceWithDeletionGlobal() {
		String sequence = "AAAGGGCCCTTTCCGGTGGGCCACTGCTTTGGGCCCAAA";
		FindSubSequence f = new FindSubSequence(anchorSequence);

		SubSequenceResultI r = f.findSequenceGlobalAlignment(sequence);

		Assert.assertEquals(anchorSequence, r.getQuerySequence());
		Assert.assertEquals(sequence, r.getTargetSequence());
		Assert.assertEquals("CCGGTGGGCCACTGC", r.getSubSequence());

		int ed = r.getEditDistance().getEditDistance();
		Assert.assertEquals(1, ed);

		int start = 13;
		int end = 27;
		int matchLen=15;

		Assert.assertEquals(start, r.getStart());
		Assert.assertEquals(end, r.getEnd());
		Assert.assertEquals(matchLen, r.getMatchLength());
	}

	@Test(enabled=true)
	public void findSequenceSimpleGlobal() {
		String querySeq = "AG";
		String targetSeq = "TAGT";

		FindSubSequence f = new FindSubSequence(querySeq);

		SubSequenceResultI r = f.findSequenceGlobalAlignment(targetSeq);
		int ed = r.getEditDistance().getEditDistance();
		Assert.assertEquals(0, ed);

		Assert.assertEquals(querySeq, r.getQuerySequence());
		Assert.assertEquals(targetSeq, r.getTargetSequence());
		Assert.assertEquals("AG", r.getSubSequence());

		int start = 2;
		int end = 3;

		Assert.assertEquals(start, r.getStart());
		Assert.assertEquals(end, r.getEnd());

	}

	@Test(enabled=true)
	public void findSequenceMissingGlobal() {
		String querySeq = "GGTTT";
		String targetSeq = "AAAAA";

		FindSubSequence f = new FindSubSequence(querySeq);

		SubSequenceResultI r = f.findSequenceGlobalAlignment(targetSeq);
		Assert.assertEquals(5, r.getMatchLength());

		Assert.assertEquals(querySeq, r.getQuerySequence());
		Assert.assertEquals(targetSeq, r.getTargetSequence());
		Assert.assertEquals(targetSeq, r.getSubSequence());

		int ed = r.getEditDistance().getEditDistance();
		Assert.assertEquals(5, ed);

		int start = 1;
		int end = 5;

		Assert.assertEquals(start, r.getStart());
		Assert.assertEquals(end, r.getEnd());

	}


	/**
	 * LOCAL ALIGNMENT TESTS
	 */
	@Test(enabled=true)
	public void findSequenceLocal() {
		String sequence = "AAAGGGCCCTTTCCGGTGGCGCCACTGCTTTGGGCCCAAA";
		FindSubSequence f = new FindSubSequence(anchorSequence);

		SubSequenceResultI r = f.findSequenceLocalAlignment(sequence);
		int ed = r.getEditDistance().getEditDistance();
		Assert.assertEquals(0, ed);

		Assert.assertEquals(anchorSequence, r.getTargetSequence());
		Assert.assertEquals(sequence, r.getQuerySequence());
		Assert.assertEquals("CCGGTGGCGCCACTGC", r.getSubSequence());

		int start = 13;
		int end = 28;
		int matchLen=16;

		Assert.assertEquals(start, r.getStart());
		Assert.assertEquals(end, r.getEnd());
		Assert.assertEquals(matchLen, r.getMatchLength());
	}

	@Test(enabled=true)
	// we don't expect this to work yet.
	public void findSequenceReverseComplimentedLocal() {
		String sequence = "AAAGGGCCCTTTCCGGTGGCGCCACTGCTTTGGGCCCAAA";
		sequence = SequenceUtil.reverseComplement(sequence);
		FindSubSequence f = new FindSubSequence(anchorSequence);

		SubSequenceResultI r = f.findSequenceLocalAlignment(sequence);
		Assert.assertEquals(anchorSequence, r.getTargetSequence());
		Assert.assertEquals(sequence, r.getQuerySequence());
		Assert.assertEquals("GCAGTGGCGCCACCGG", r.getSubSequence());


		int ed = r.getEditDistance().getEditDistance();
		Assert.assertEquals(4, ed);

		int start = 13;
		int end = 28;
		int matchLen=16;

		Assert.assertEquals(start, r.getStart());
		Assert.assertEquals(end, r.getEnd());
		Assert.assertEquals(matchLen, r.getMatchLength());
	}


	@Test(enabled=true)
	// has a 1 BP insertion in the anchor sequence.
	public void findSequenceWithInsertionLocal() {
		String sequence = "AAAGGGCCCTTTCCGGTGGCCGCCACTGCTTTGGGCCCAAA";
		FindSubSequence f = new FindSubSequence(anchorSequence);

		SubSequenceResultI r = f.findSequenceLocalAlignment(sequence);
		int ed = r.getEditDistance().getEditDistance();
		Assert.assertEquals(1, ed);

		Assert.assertEquals(anchorSequence, r.getTargetSequence());
		Assert.assertEquals(sequence, r.getQuerySequence());
		Assert.assertEquals("CCGGTGGCCGCCACTGC", r.getSubSequence());

		int start = 13;
		int end = 29;
		int matchLen=17;

		Assert.assertEquals(start, r.getStart());
		Assert.assertEquals(end, r.getEnd());
		Assert.assertEquals(matchLen, r.getMatchLength());
	}

	@Test(enabled=true)
	// has a 1 BP insertion in the anchor sequence.
	public void findSequenceWithDeletionLocal() {
		String sequence = "AAAGGGCCCTTTCCGGTGGGCCACTGCTTTGGGCCCAAA";
		FindSubSequence f = new FindSubSequence(anchorSequence);

		SubSequenceResultI r = f.findSequenceLocalAlignment(sequence);

		Assert.assertEquals(anchorSequence, r.getTargetSequence());
		Assert.assertEquals(sequence, r.getQuerySequence());
		Assert.assertEquals("CCGGTGGGCCACTGC", r.getSubSequence());

		int ed = r.getEditDistance().getEditDistance();
		Assert.assertEquals(1, ed);

		int start = 13;
		int end = 27;
		int matchLen=15;

		Assert.assertEquals(start, r.getStart());
		Assert.assertEquals(end, r.getEnd());
		Assert.assertEquals(matchLen, r.getMatchLength());
	}

	@Test(enabled=true)
	public void findSequenceSimpleLocal() {
		String querySeq = "AG";
		String targetSeq = "TAGT";

		FindSubSequence f = new FindSubSequence(querySeq);

		SubSequenceResultI r = f.findSequenceLocalAlignment(targetSeq);
		int ed = r.getEditDistance().getEditDistance();
		Assert.assertEquals(0, ed);

		Assert.assertEquals(querySeq, r.getTargetSequence());
		Assert.assertEquals(targetSeq, r.getQuerySequence());
		Assert.assertEquals("AG", r.getSubSequence());

		int start = 2;
		int end = 3;

		Assert.assertEquals(start, r.getStart());
		Assert.assertEquals(end, r.getEnd());

	}

	@Test(enabled=true)
	public void findSequenceMissingLocal() {
		String querySeq = "GGTTT";
		String targetSeq = "AAAAA";
		String result = "";

		FindSubSequence f = new FindSubSequence(targetSeq);

		SubSequenceResultI r = f.findSequenceLocalAlignment(querySeq);


		Assert.assertEquals(querySeq, r.getQuerySequence());
		Assert.assertEquals(targetSeq, r.getTargetSequence());
		Assert.assertEquals(result, r.getSubSequence());
		Assert.assertEquals(0, r.getMatchLength());

		int ed = r.getEditDistance().getEditDistance();
		Assert.assertEquals(5, ed);

		int start = -1;
		int end = -1;

		Assert.assertEquals(start, r.getStart());
		Assert.assertEquals(end, r.getEnd());

	}



}
