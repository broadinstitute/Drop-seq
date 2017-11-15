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
package org.broadinstitute.dropseqrna.utils.editdistance;

import org.apache.commons.lang.StringUtils;
import org.junit.Assert;
import org.testng.annotations.Test;

public class ModifiedLevenshteinDistanceTest {

	@Test(enabled=true, groups={"transcriptome", "dropseq"})
	public void websiteTest () {
		// see http://stackoverflow.com/questions/5849139/levenshtein-distance-inferring-the-edit-operations-from-the-matrix

		LevenshteinDistanceResult r =  LevenshteinDistance.computeLevenshteinDistanceResult("democrat", "republican");
		int ed = r.getEditDistance();
		Assert.assertEquals(8, ed);
		String [] operations = r.getOperations();
		String result = StringUtils.join(operations);
		Assert.assertEquals(result, "SMIISSSSMS");

	}

	@Test(enabled=true, groups={"transcriptome", "dropseq"})
	public void simpleTests1 () {
		LevenshteinDistanceResult r =  LevenshteinDistance.computeLevenshteinDistanceResult("AAA", "AAA");
		int ed = r.getEditDistance();
		Assert.assertEquals(0, ed);

		r =  LevenshteinDistance.computeLevenshteinDistanceResult("AAA", "ATA");
		ed = r.getEditDistance();
		Assert.assertEquals(1, ed);

		r =  LevenshteinDistance.computeLevenshteinDistanceResult("AAA", "ATA", 1,1,2);
		ed = r.getEditDistance();
		Assert.assertEquals(2, ed);

		// test deletion
		r =  LevenshteinDistance.computeLevenshteinDistanceResult("AATA", "AAA");
		ed = r.getEditDistance();
		Assert.assertEquals(1, ed);

		// test insertion
		r =  LevenshteinDistance.computeLevenshteinDistanceResult("AAA", "AATA");
		ed = r.getEditDistance();
		Assert.assertEquals(1, ed);

		// test deletion - modify deletion weight
		r =  LevenshteinDistance.computeLevenshteinDistanceResult("AATA", "AAA", 2,1,1);
		ed = r.getEditDistance();
		Assert.assertEquals(2, ed);

		// test insertion - modify insertion weight
		r =  LevenshteinDistance.computeLevenshteinDistanceResult("AAA", "AATA", 1,2,1);
		ed = r.getEditDistance();
		Assert.assertEquals(2, ed);


	}

	@Test(enabled=true, groups={"transcriptome", "dropseq"})
	public void testEditDistanceOperations () {
		LevenshteinDistanceResult r =  LevenshteinDistance.computeLevenshteinDistanceResult("AAA", "AATA");
		int ed = r.getEditDistance();
		Assert.assertEquals(1, ed);
		String [] operations = r.getOperations();
		String result = StringUtils.join(operations);
		Assert.assertEquals("MMIM", result);

	}


	@Test(enabled=true, groups={"transcriptome", "dropseq"})
	public void testEditDistanceOperationsFull () {

		testEditDistance("AGGCTAGGTGTT","AGGCTAAGGTGT", 1);
		testEditDistance("AGGCTAGGTGTT","AGCTAGGTGTTT",1);
		testEditDistance("AGGCTAGGTGTT","AGCTAGCTGTTT",2);

		// 1 indel, 2 subs
		testEditDistance("AGGCTAGGTGTT","ACCTAGCTGTTT",3);

		// 1 indel, 0 subs
		testEditDistance("AAAGAGGACGG","AAGAGGACGCG",1);


	}

	@Test(enabled=true, groups={"transcriptome", "dropseq"})
	public void testEditDistanceOperationsFull2 () {

		// 2 bp insertion
		testEditDistance("GCT", "TAG", 2);

		// 2 bp substitution
		testEditDistance("AATG","ZZTG",2);

		// 2 bp insertion of a "GC"
		testEditDistance("TCCCAGCAATAG","TCCCAGCAAGCT",2);

		// has an indel of length 1 and a substitution of the first base.
		testEditDistance("TATTTTACGCAG","AATTTTAACGCA",2);

		// 2 bp insertion, 1 substiution
		testEditDistance("ACCCAGCAATAG","TCCCAGCAAGCT",3);

		// go bigger.

		// 4 bp insertion, 2 substiution.  Turns out this is ED=4, 'cause long runs of a base "lose" info.
		testEditDistance("TCCCAGCAAGCT","TCCAAAACATTA",4);

		// 3 bp insertion, 3bp substitution.
		testEditDistance("AAACCATTTGGG","GAGCCCCCCTTT",6);



	}

	private void testEditDistance (final String a, final String b, final int expectedED) {
		LevenshteinDistanceResult r = LevenshteinDistance.computeLevenshteinDistanceResult(a,b, 1,1,2);

		int ed2 = r.getEditDistanceIndelCorrected();
		Assert.assertEquals(expectedED, ed2);

	}

	@Test(enabled=true)
	public void testIndelOnly() {
		String bc1="AAAAGTGAGGAC";
		String [] bc2 ={"AAAAGTGAGGCA", "AAAAGTGAGGCC", "AAAAGTGAGGCG", "AAAAGTGAGGCT"};
		String [] opsExpected={"M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "D", "M", "I"};
		for (String element : bc2) {
			LevenshteinDistanceResult r = LevenshteinDistance.computeLevenshteinDistanceResult(bc1,element, 1,1,2);
			int ed1 = r.getEditDistance();
			int ed11= r.getEditDistanceIndelCorrected();

			Assert.assertEquals(ed1, 2);
			Assert.assertEquals(ed11, 1);
			String []  ops = r.getOperations();
			for (int i=0; i<ops.length; i++)
				Assert.assertEquals(opsExpected[i], ops[i]);
		}
	}


}
