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

import org.broadinstitute.dropseqrna.sbarro.utils.ConsensusSequence;
import org.broadinstitute.dropseqrna.sbarro.utils.ConsensusSequenceFactory;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistance;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ConsensusSequenceFactoryTest {


	@Test(enabled=true)
	public void test3WithBQ () {
		String seq1 = "GATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCTACGGATCAGCCCCGGTGGCGCCACTGCTGTAATTAGGCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";
		String bq1 =  "BCCCCFFFFFBCGGGGGGGGGGGGHHHHHHHHHHHHHH3B2BE2B333BEEGGGGHGGGGGHHHH3434BF4B3BFGHHHHHHGGFGHGHHHGGGHHHHHHHHHHHEHHHHHHHHHGHHHHGGHHHHHHHHHHHGHHHHHHHHHHHGGGGHHHFHHHGG";
		String seq2 = "AAGTCTCAGGTTGAGAAGTGTTGCCAGTTGTTCTTCTGATCTAATGTTTTTTTCTCGACTGAAAAGCCCAATCGCAGCAGTGGCGCCACCGGAGCTGAGTCGTAGCTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCAGATCGGAAGAGCGTCGTGTA";
		String bq2  = "3AABBFFFFFFFGGGGGGGGGGFHHHHHHHHHHHHHHHHHHHHHHHGHHHGGGHHHGHGGGHHHHH1G1AG231??EEHHHHHGGGGGGHGG///03?04B>?FFGHHHHHHHHHHHHHHHGFDGGHEGHGGGCCGHHHGHEFHHHG?FG0>FCECD<EGG";


		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(seq1, seq2, false);
		cs.addReadBaseQualities(bq1, bq2);
		String result = cs.getConsensusSequence();

		String consensus = "TACACGACGCTCTTCCGATCTGATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCTACGAATCAGCCCCGGTGGCGCCACTGCTGTAATTAGGCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";
		Assert.assertEquals(result, consensus);

		int ed = cs.getLocalAlignmentEditDistance().getEditDistance();
		Assert.assertEquals(ed, 6);
	}

	/**
	 * The first read has leading bases, the second read has trailing bases.
	 *   CCGATTTTTGATCACTCTTTCCCC
              CCCCGATCACTCTAGGGGGATCGGAAGAG
	 *
	 */
	@Test(enabled=true)
	public void testConsensusSequenceIndex1 () {
		String seq1 = "CCGATTTTTGATCACTCTTTCCCC";
		String bq1 =  "BCCCCFFFFFBCGGGGGGGGGGGG";
		String seq2 = "CCCCGATCACTCTAGGGGGATCGGAAGAG";
		String bq2  = "3AABBFFFFFFF3333333333FFFFFFF";
		String consensus="CCGATTTTTGATCACTCTTTCCCCATCGGAAGAG";

		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(seq1, seq2, false);
		cs.addReadBaseQualities(bq1, bq2);
		String result = cs.getConsensusSequence();
		Assert.assertEquals(result, consensus);

	}

	/**
	 * Same as the above but the reads are flipped for which is first vs second.
	 *
     *        CCCCGATCACTCTAGGGGGATCGGAAGAG
	 *	 CCGATTTTTGATCACTCTTTCCCC
	 */
	@Test(enabled=true)
	public void testConsensusSequenceIndex2 () {
		String seq2 = "CCGATTTTTGATCACTCTTTCCCC";
		String bq2 =  "BCCCCFFFFFBCGGGGGGGGGGGG";
		String seq1 = "CCCCGATCACTCTAGGGGGATCGGAAGAG";
		String bq1  = "3AABBFFFFFFF3333333333FFFFFFF";

		String consensus="CCGATTTTTGATCACTCTTTCCCCATCGGAAGAG";

		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(seq1, seq2, false);
		cs.addReadBaseQualities(bq1, bq2);
		String result = cs.getConsensusSequence();
		Assert.assertEquals(result, consensus);

	}

	/**
	 * Add insertion into first read.
	 *
     *        CCCCGATCATCTCTAGGGGGATCGGAAGAG
	 *	 CCGATTTTTGATCA-CTCTTTCCCC
	 */
	@Test(enabled=true)
	public void testConsensusSequenceIndex3 () {
		String seq2 = "CCGATTTTTGATCATCTCTTTCCCC";
		String bq2 =  "BCCCCFFFFFBCGGGGGGGGGGGGG";
		String seq1 = "CCCCGATCACTCTAGGGGGATCGGAAGAG";
		String bq1  = "3AABBFFFFFFF3333333333FFFFFFF";


		String consensus="CCGATTTTTGATCATCTCTTTCCCCATCGGAAGAG";

		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(seq1, seq2, false);
		String seq2RC = cs.getOriginalReadTwo();
		cs.addReadBaseQualities(bq1, bq2);
		String result = cs.getConsensusSequence();
		Assert.assertEquals(result, consensus);

	}

	@Test(enabled=true)
	/**
	 * Test one sequence being a subsequence of the other
	 * CCGATTTTTGATCACTCTAATTCCCC
	 *      CCCCGATCACTCTAAAGGG
	 */
	public void testConsensusSequenceIndex4 () {
		String seq1 = "CCGATTTTTGATCACTCTAATTCCCC";
		String bq1 =  "BCCCCFFFFFBCGGGGGGGGGGGGGG";
		String seq2 = "CCCCGATCACTCTAAAGGG";
		String bq2  = "3AABBFFFFFFFGGGGGGG";

		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(seq1, seq2, false);
		cs.addReadBaseQualities(bq1, bq2);
		String result = cs.getConsensusSequence();
		String consensus = "CCGATTTTTGATCACTCTAATTCCCC";
		Assert.assertEquals(result, consensus);
	}


	@Test(enabled=true)
	/**
	 * Set up with trailing disagreeing bases on both sides to test local alignment.
	 */
	public void test3WithBQSecondReadLeft () {
		String seq1 = "AAGATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCTACGGATCAGCCCCGGTGGCGCCACTGCTGTAATTAGGCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";
		String bq1 =  "CCBCCCCFFFFFBCGGGGGGGGGGGGHHHHHHHHHHHHHH3B2BE2B333BEEGGGGHGGGGGHHHH3434BF4B3BFGHHHHHHGGFGHGHHHGGGHHHHHHHHHHHEHHHHHHHHHGHHHHGGHHHHHHHHHHHGHHHHHHHHHHHGGGGHHHFHHHGG";
		String seq2 = "AAAAGTCTCAGGTTGAGAAGTGTTGCCAGTTGTTCTTCTGATCTAATGTTTTTTTCTCGACTGAAAAGCCCAATCGCAGCAGTGGCGCCACCGGAGCTGAGTCGTAGCTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCAGATCGGAAGAGCGTCGTGTA";
		String bq2  = "FF3AABBFFFFFFFGGGGGGGGGGFHHHHHHHHHHHHHHHHHHHHHHHGHHHGGGHHHGHGGGHHHHH1G1AG231??EEHHHHHGGGGGGHGG///03?04B>?FFGHHHHHHHHHHHHHHHGFDGGHEGHGGGCCGHHHGHEFHHHG?FG0>FCECD<EGG";


		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(seq1, seq2, false);
		cs.addReadBaseQualities(bq1, bq2);
		String result = cs.getConsensusSequence();

		String consensus = "TACACGACGCTCTTCCGATCTGATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCTACGAATCAGCCCCGGTGGCGCCACTGCTGTAATTAGGCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";
		Assert.assertEquals(result, consensus);

		int ed = cs.getLocalAlignmentEditDistance().getEditDistance();
		Assert.assertEquals(ed, 6);
	}


	// A weird read from real data.
	// M00282:312:000000000-AUK2C:1:1101:10000:14973 from GFPBC_unmapped_adapter_marked.bam
	@Test(enabled=true)
	public void test4() {
		String seq1 = "TCGCTTGGTCAACCCCTCAGCGGCAAAAATTAAAATTTTTACCGCTTCGGCGTTATAACCTCACACTCAATCTTTTATCACGAAGTCATGATTGAATCGCGAGTGGTCGGCAGATTGCGATAAACGGTCACATTAAATTTAACCTGACTATTCCACTGCAA";
		String bq1 =  "AAAABABAFFFFGGEECEGGGGGCEEHBEFGGHHFGHHDE2F32AAEAEEEGGGGHHHHHF?GG3FGFCFHGGBBGEHFHHGC1E?GHDGHFGHFHHF?CEEC@HEFHGCCGCEHHFHDCCFAF1DGAEDH1FFGHHHHHHGF1FGGFDDGBGB0DD<DDD";
		String seq2 = "ACTTGCCGCCGCGTGAAATTTCTATGAAGGATGTTTTCCGTTCTGGTGATTCGTCTAAGAAGTTTAAGATTGCTGAGGGTCAGTGGTATCGTTATGCGCCTTCGTATGTTTCTCCTGCTTATCACCTTCTTGAAGGCTTCCCATTCATTCAGGAACCGC";
		String bq2  = "BBBBB5CDDDBDGEGCFGFGGGHFFGHHCHHH3BG5GDFA3BGDHHGFFHHHGAA1AFFFBFGFGHF5@GBEGGHFGEEEHFH4GB4B4FGHHEEBEFFAG/F?FGHFBFG4F?G33?3BDGBFGGGHH1112FFHAFGH10GG1@@GFFHEBBCA/@/";

		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(seq1, seq2, false);
		cs.addReadBaseQualities(bq1, bq2);
		String result = cs.getConsensusSequence();
		int ed = cs.getLocalAlignmentEditDistance().getEditDistance();
		Assert.assertEquals(ed, 75);

	}

	@Test(enabled=true)
	public void testExtremeTrimmedRead1 () {
		String seq1 = "A";
		String bq1 =  "A";
		String seq2 = "A";
		String bq2  = "B";

		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(seq1, seq2, true);
		cs.addReadBaseQualities(bq1, bq2);
		String result = cs.getConsensusSequence();
		int ed = cs.getLocalAlignmentEditDistance().getEditDistance();
		Assert.assertEquals(ed, 0);
	}

	@Test(enabled=true, expectedExceptions=IllegalArgumentException.class)
	public void testExtremeTrimmedRead2 () {
		String seq1 = "";
		String bq1 =  "";
		String seq2 = "";
		String bq2  = "";

		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(seq1, seq2, true);
		cs.addReadBaseQualities(bq1, bq2);
		String result = cs.getConsensusSequence();
		int ed = cs.getLocalAlignmentEditDistance().getEditDistance();
		Assert.assertEquals(ed, 0);
	}



	@Test(enabled=true)
	public void testExtractIndividualReads () {
		String seq1 = "GATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCTACGGATCAGCCCCGGTGGCGCCACTGCTGTAATTAGGCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";
		String bq1 =  "BCCCCFFFFFBCGGGGGGGGGGGGHHHHHHHHHHHHHH3B2BE2B333BEEGGGGHGGGGGHHHH3434BF4B3BFGHHHHHHGGFGHGHHHGGGHHHHHHHHHHHEHHHHHHHHHGHHHHGGHHHHHHHHHHHGHHHHHHHHHHHGGGGHHHFHHHGG";
		String seq2 = "AAGTCTCAGGTTGAGAAGTGTTGCCAGTTGTTCTTCTGATCTAATGTTTTTTTCTCGACTGAAAAGCCCAATCGCAGCAGTGGCGCCACCGGAGCTGAGTCGTAGCTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCAGATCGGAAGAGCGTCGTGTA";
		String bq2  = "3AABBFFFFFFFGGGGGGGGGGFHHHHHHHHHHHHHHHHHHHHHHHGHHHGGGHHHGHGGGHHHHH1G1AG231??EEHHHHHGGGGGGHGG///03?04B>?FFGHHHHHHHHHHHHHHHGFDGGHEGHGGGCCGHHHGHEFHHHG?FG0>FCECD<EGG";

		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(seq1, seq2, false);
		cs.addReadBaseQualities(bq1, bq2);
		String result = cs.getConsensusSequence();

		String consensus = "TACACGACGCTCTTCCGATCTGATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCTACGAATCAGCCCCGGTGGCGCCACTGCTGTAATTAGGCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";
		Assert.assertEquals(result, consensus);

		// OK, now let's extract the sequence where some of the sequence differences are, from 62 to 88.  Those should match the original sequence.
		String seq1Fragment = cs.getOriginalSequenceAtConsensusLocation(1, 62, 95);
		String seq2Fragment = cs.getOriginalSequenceAtConsensusLocation(2, 62, 95);

		Assert.assertNotEquals(seq1Fragment, seq2Fragment);

		// each should be a substring of the original sequence.
		int indexSeq1=seq1.indexOf(seq1Fragment);
		String seq2Fixed = seq2;
		if (cs.isReadTwoReverseComplimented())
			seq2Fixed = SequenceUtil.reverseComplement(seq2Fixed);

		int indexSeq2=seq2Fixed.indexOf(seq2Fragment);

		Assert.assertEquals(indexSeq1, 40);
		Assert.assertEquals(indexSeq2, 61);

		int ed = LevenshteinDistance.computeLevenshteinDistanceResult(seq1Fragment, seq2Fragment).getEditDistance();
		Assert.assertEquals(ed, 6);
	}

	@Test(enabled=true)
	public void testJankyRead1 () {
		// read trimmed from original M00282:312:000000000-AUK2C:1:1101:10000:14973
		String seq1 = "TCGCTTGGTCAACCCCTCAGCGGCAAAAATTAAAATTTTTACCGCTTCGGCGTTATAACCTCACACTCAATCTTTTATCACGAAGTCATGATTGAATCGCGAGTGGTCGGCAGATTGCGATAAACGGTCACATTAAATTTAACCTGACTATTCCACTGCAA";
		String bq1 =  "AAAABABAFFFFGGEECEGGGGGCEEHBEFGGHHFGHHDE2F32AAEAEEEGGGGHHHHHF?GG3FGFCFHGGBBGEHFHHGC1E?GHDGHFGHFHHF?CEEC@HEFHGCCGCEHHFHDCCFAF1DGAEDH1FFGHHHHHHGF1FGGFDDGBGB0DD<DDD";
		String seq2 = "ACTTGCCGCCGCGTGAAATTTCTATGAAGGATGTTTTCCGTTCTGGTGATTCGTCTAAGAAGTTTAAGATTGCTGAGGGTCAGTGGTATCGTTATGCGCCTTCGTATGTTTCTCCTGCTTATCACCTTCTTGAAGGCTTCCCATTCATTCAGGAACCGC";
		String bq2  = "BBBBB5CDDDBDGEGCFGFGGGHFFGHHCHHH3BG5GDFA3BGDHHGFFHHHGAA1AFFFBFGFGHF5@GBEGGHFGEEEHFH4GB4B4FGHHEEBEFFAG/F?FGHFBFG4F?G33?3BDGBFGGGHH1112FFHAFGH10GG1@@GFFHEBBCA/@/";

		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(seq1, seq2, true);
		cs.addReadBaseQualities(bq1, bq2);
		seq2 = cs.getOriginalReadTwo();
		String result = cs.getConsensusSequence();


	}

	/*
	@Test(enabled=false)
	public void test1NoBQ() {
		String seq1 = "GATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCTACGGATCAGCCCCGG";
		String seq2 = "TACACGACGCTCTTCCGATCTGATCACTCTCGGCATGGACGAGCTGTACAAG";
		String consensus="TACACGACGCTCTTCCGATCTGATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCTACGGATCAGCCCCGG";


		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(seq1, seq2);
		Assert.assertNotNull(cs);
		String result = cs.getConsensusSequence();
		Assert.assertEquals(result, consensus);

		int ed = cs.getLocalAlignmentEditDistance().getEditDistance();
		Assert.assertEquals(ed, 0);

	}
	*/
	// @Test(enabled=false)
	/**
	 * This is read @M00282:312:000000000-AUK2C:1:1101:13302:1744 1:N:0:1 from file
	 * /broad/mccarroll/arpy/sequence/rabies_barcoding/160919_plasmidBC_pSPBN-GFP_Rep1/GFPBC_S1_L001_R1_001.fastq.gz
	 */
	/*
	public void test2NoBQ () {
		String seq1 = "GATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCTACGGATCAGCCCCGGTGGCGCCACTGCTGTAATTAGGCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";
		String seq2 = "AAGTCTCAGGTTGAGAAGTGTTGCCAGTTGTTCTTCTGATCTAATGTTTTTTTCTCGACTGAAAAGCCCAATCGCAGCAGTGGCGCCACCGGAGCTGAGTCGTAGCTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCAGATCGGAAGAGCGTCGTGTA";
		String consensus = "TACACGACGCTCTTCCGATCTGATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCTACGGATCAGCCCCGGTGGCGCCACTGCTGTAATTAGGCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";
		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(seq1, seq2);
		Assert.assertNotNull(cs);
		String result = cs.getConsensusSequence();
		Assert.assertEquals(result, consensus);

		int ed = cs.getLocalAlignmentEditDistance().getEditDistance();
		Assert.assertEquals(ed, 6);
	}
	*/


}
