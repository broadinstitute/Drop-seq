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

import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistance;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistanceResult;
import org.testng.Assert;
import org.testng.annotations.Test;


public class TrimSequenceTemplateTest {

	/**
	@Test(enabled=false)
	public void smithWatermanAlignmentTest () throws CompoundNotFoundException {
		String targetSeq = "GGGGGGGGGGGG";

		DNASequence target = new DNASequence(targetSeq,
				AmbiguityDNACompoundSet.getDNACompoundSet());

		//String querySeq = "GTACTCTGCGTTGATACCACTGCTTAAAA";
		String querySeq = "GGGGGAAGGGGG";

		DNASequence query = new DNASequence(querySeq,
				AmbiguityDNACompoundSet.getDNACompoundSet());

		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
		SimpleGapPenalty gapP = new SimpleGapPenalty();
		gapP.setOpenPenalty((short)5);
		gapP.setExtensionPenalty((short)2);
		SequencePair<DNASequence, NucleotideCompound> psa =
				Alignments.getPairwiseAlignment(query, target,
						PairwiseSequenceAlignerType.LOCAL, gapP, matrix);

		System.out.println(psa);

		Iterator<AlignedSequence<DNASequence, NucleotideCompound>> itr = psa.iterator();
		while (itr.hasNext()) {
			AlignedSequence<DNASequence, NucleotideCompound> a = itr.next();
			System.out.println(a);
		}



	}
	*/
	@Test(enabled=true)
	public void testPosInRead () {
		String read = "GGTCTTATGTGTCTACGTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGA";
		TrimSequenceTemplate t = new TrimSequenceTemplate("GTACTCTGCGTTGATACCACTGCTT");
		int pos = t.getPositionInRead(read, 8, 0);
		Assert.assertEquals(pos, 16);

	}

	@Test(enabled=true)
	public void testPosInRead2 () {
		String read = "AGTACTCTG";
		TrimSequenceTemplate t = new TrimSequenceTemplate("GTACTCTGCGTTGATACCACTGCTT");
		int pos = t.getPositionInRead(read, 8, 0);
		Assert.assertEquals(pos, 1);
	}

	//CCTATAAGCAGTGTATCATAAATTTAAATCACTACAAACACTCAACCACTACATACAAAC
	@Test(enabled=true)
	public void testPosInRead3 () {
		String read = "CCTATAAGCAGTGTATCATAAATTTAAATCACTACAAACACTCAACCACTACATACAAAC";
		TrimSequenceTemplate t = new TrimSequenceTemplate("AAGCAGTGGTATCAACGCAGAGTAC");
		int pos = t.getPositionInRead(read, 8, 0);
		Assert.assertEquals(pos, 5);
	}

	//AGCAGTGGATGTGAGTGCTTCAAGTTATGGACAGACAAGATCAGAAGCACACATAATCTA
	@Test(enabled=true)
	public void testPosInRead4 () {
		String read = "AGCAGTGGATGTGAGTGCTTCAAGTTATGGACAGACAAGATCAGAAGCACACATAATCTA";
		//TrimSequenceTemplate t = new TrimSequenceTemplate("AAGCAGTGGTATCAACGCAGAGTAC");
		//int pos = t.getPositionInRead(read, 8, 0);
		//Assert.assertEquals(pos, 1);
	}

	@Test(enabled=true)
	public void testPosInRead5 () {
		String read = "AGCAGTGGATGTGAGTGCTTCAAGTTATGGACAGACAAGATCAGAAGCACACATAATCTA";
		String template="AAGCAGTGGTATCAACGCAGAGTAC";

		LevenshteinDistanceResult r = LevenshteinDistance.computeLevenshteinDistanceResult(read, template, 1, 1, 1);

		String [] rr = r.getOperations();

		Assert.assertNotNull(rr);
		//TrimSequenceTemplate t = new TrimSequenceTemplate("AAGCAGTGGTATCAACGCAGAGTAC");
		//int pos = t.getPositionInRead(read, 8, 0);
		//Assert.assertEquals(pos, 1);
	}

	@Test(enabled=true)
	public void testTrimExample () {

		String primerSeq="ACGGCGATCGCAAGCAGTGGTATCAACGCAGAGTGAATGGG";
		String readSeq=null;

		TrimSequenceTemplate t = new TrimSequenceTemplate(primerSeq);

		readSeq="TATCAACGCAGAGTGAATGGGCCGTATGTACAAAGGGCAGGGACTTAA";
		int pos = t.getPositionInTemplate(readSeq, 4, 0);
		Assert.assertEquals(pos, 20);

		readSeq="ACGGCGATCGCAAGCAGT";
		pos = t.getPositionInTemplate(readSeq, 1, 0);
		Assert.assertEquals(pos, 0);

		readSeq="TTTTTTTTTTTT";
		pos = t.getPositionInTemplate(readSeq, 1, 0);
		Assert.assertEquals(pos, -1);



	}

	@Test
	public void testHasForwardMatch() {
		String barcode="TATCAACNCAGAGTGA";
		String readSeq=null;

		TrimSequenceTemplate t = new TrimSequenceTemplate(barcode, "N");

		readSeq="TATCAACGCAGAGTGAATGGGCCGTATGTACAAAGGGCAGGGACTTAA";
		boolean test =t.hasForwardMatch(readSeq);
		Assert.assertTrue(test);

		// not at the start of the read.
		readSeq="CAGAGTGAATGGGCCGTATGTACAAAGGGCAGGGACTTAA";
		test =t.hasForwardMatch(readSeq);
		Assert.assertFalse(test);

	}



	@Test(enabled=true)
	public void test2 () {
		String primerSeq="ACGGCGATCGCAAGCAGTGGTATCAACGCAGAGTGAATGGG";
		TrimSequenceTemplate t = new TrimSequenceTemplate(primerSeq);

		String readSeq = "CTTGTGCAGCAATGGCCAAGATCAAGGCTCGAGATCTTCGCGGGAAGAAG";
		int pos = t.getPositionInTemplate(readSeq, 4, 1);
		Assert.assertEquals(pos, -1);
	}

	@Test(enabled=true)
	public void test3() {
		String primerSeq="AAGCAGTGGTATCAACGCAGAGTGAATGGG";
		TrimSequenceTemplate t = new TrimSequenceTemplate(primerSeq);

		String readSeq = "CAGTGGTATCAACGCAGAGTGAATGGGTATGGTGAAACCCCGTCTCCACTAAAAATACAAACATTAGCTGGGCGTGGTGGTGCGCGCCTGTAATCCCAGCTACTCAAGAGGCTGAGGCAGGAGCATCGCTTGAACCTGGGAGGCGGAGGT";
		int pos = t.getPositionInTemplate(readSeq, 4, 1);
		Assert.assertEquals(pos, 3);

	}

	@Test(enabled=true)
	public void test4() {
		String primerSeq="AAGCAGTGGTATCAACGCAGAGTGAATGGG";
		TrimSequenceTemplate t = new TrimSequenceTemplate(primerSeq);

		String readSeq = "GTGAATGGGAGCAGAGGGCGGCGGCGGTGCGGGCGGACCCGGGTCCCTAA";
		int pos = t.getPositionInTemplate(readSeq, 4, 1);
		Assert.assertEquals(pos, 21);

	}

	@Test(enabled=true)
	public void test5() {
		String primerSeq="AAGCAGTGGTATCAACGCAGAGTGAATGGG";
		TrimSequenceTemplate t = new TrimSequenceTemplate(primerSeq);

		String readSeq = "GAATGGGAGATCGCACGTGGCGCCCGAGAAGTAGTGGTGATCCCGAGACC";
		int pos = t.getPositionInTemplate(readSeq, 4, 1);
		Assert.assertEquals(pos, 23);

	}


}
