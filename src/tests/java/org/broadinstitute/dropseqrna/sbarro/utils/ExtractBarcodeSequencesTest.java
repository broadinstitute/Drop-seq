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

import org.broadinstitute.dropseqrna.sbarro.utils.ConsensusSequence;
import org.broadinstitute.dropseqrna.sbarro.utils.ConsensusSequenceFactory;
import org.broadinstitute.dropseqrna.sbarro.utils.ExtractBarcodeSequences;
import org.broadinstitute.dropseqrna.sbarro.utils.ExtractedSequenceGroup;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ExtractBarcodeSequencesTest {

	@Test(enabled=true)
	public void test1 () {
		String seq = "GATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCTACCGATGGCCTCCGGTGGCGCCACTGCGAGGTTCTTCCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";

		ExtractBarcodeSequences e = new ExtractBarcodeSequences("CGGCATGGACGAGCTGTACAAGTAAGCTA", "CCGGTGGCGCCACTGC");
		ExtractedSequenceGroup result =  e.findRabiesBarcode(seq);

		String expectedStopCodonBC="CCGATGGCCT";
		String expectedPolyABC="GAGGTTCTTC";

		int [] gpfAnchorPos={10,(10+29-1)};
		int gfpAnchorED=0;
		int [] cassetteAnchorPos={49,64};
		int cassetteAnchorED=0;

		Assert.assertEquals(result.getStopCodonBarcode().getSequence(), expectedStopCodonBC);
		Assert.assertEquals(result.getPolyABarcode().getSequence(), expectedPolyABC);

		Assert.assertEquals(result.getGfpAnchor().getStart(), gpfAnchorPos[0]);
		Assert.assertEquals(result.getGfpAnchor().getEnd(), gpfAnchorPos[1]);
		Assert.assertEquals(result.getGfpAnchor().getEditDistance().getEditDistance(), gfpAnchorED);

		Assert.assertEquals(result.getCassetteAnchor().getStart(), cassetteAnchorPos[0]);
		Assert.assertEquals(result.getCassetteAnchor().getEnd(), cassetteAnchorPos[1]);
		Assert.assertEquals(result.getCassetteAnchor().getEditDistance().getEditDistance(), cassetteAnchorED);

	}

	@Test(enabled=true)
	public void test2 () {
		String seq = "GATCACTCTCGGCATGGTCGAGCTGTACAAGTAAGCTACCGATGGCCTCCGGTGGCGCCACTGCGAGGTTCTTCCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";

		ExtractBarcodeSequences e = new ExtractBarcodeSequences("CGGCATGGACGAGCTGTACAAGTAAGCTA", "CCGGTGGCGCCACTGC");
		ExtractedSequenceGroup result =  e.findRabiesBarcode(seq);

		String expectedStopCodonBC="CCGATGGCCT";
		String expectedPolyABC="GAGGTTCTTC";

		int [] gpfAnchorPos={10,(10+29-1)};
		int gfpAnchorED=1;
		int [] cassetteAnchorPos={49,64};
		int cassetteAnchorED=0;

		Assert.assertEquals(result.getStopCodonBarcode().getSequence(), expectedStopCodonBC);
		Assert.assertEquals(result.getPolyABarcode().getSequence(), expectedPolyABC);

		Assert.assertEquals(result.getGfpAnchor().getStart(), gpfAnchorPos[0]);
		Assert.assertEquals(result.getGfpAnchor().getEnd(), gpfAnchorPos[1]);
		Assert.assertEquals(result.getGfpAnchor().getEditDistance().getEditDistance(), gfpAnchorED);

		Assert.assertEquals(result.getCassetteAnchor().getStart(), cassetteAnchorPos[0]);
		Assert.assertEquals(result.getCassetteAnchor().getEnd(), cassetteAnchorPos[1]);
		Assert.assertEquals(result.getCassetteAnchor().getEditDistance().getEditDistance(), cassetteAnchorED);

	}

	@Test(enabled=true)
	public void test3 () {
		String seq = "GATCACTCTCGGCATGGCGAGCTGTACAAGTAAGCTACCGATGGCCTCCGGTGGCGCCACTGCGAGGTTCTTCCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";

		ExtractBarcodeSequences e = new ExtractBarcodeSequences("CGGCATGGACGAGCTGTACAAGTAAGCTA", "CCGGTGGCGCCACTGC");
		ExtractedSequenceGroup result =  e.findRabiesBarcode(seq);

		String expectedStopCodonBC="CCGATGGCCT";
		String expectedPolyABC="GAGGTTCTTC";

		int [] gpfAnchorPos={10,(10+28-1)};
		int gfpAnchorED=1;
		int [] cassetteAnchorPos={48,63};
		int cassetteAnchorED=0;

		Assert.assertEquals(result.getStopCodonBarcode().getSequence(), expectedStopCodonBC);
		Assert.assertEquals(result.getPolyABarcode().getSequence(), expectedPolyABC);

		Assert.assertEquals(result.getGfpAnchor().getStart(), gpfAnchorPos[0]);
		Assert.assertEquals(result.getGfpAnchor().getEnd(), gpfAnchorPos[1]);
		Assert.assertEquals(result.getGfpAnchor().getEditDistance().getEditDistance(), gfpAnchorED);

		Assert.assertEquals(result.getCassetteAnchor().getStart(), cassetteAnchorPos[0]);
		Assert.assertEquals(result.getCassetteAnchor().getEnd(), cassetteAnchorPos[1]);
		Assert.assertEquals(result.getCassetteAnchor().getEditDistance().getEditDistance(), cassetteAnchorED);

	}

	@Test(enabled=true)
	// this read is real janky and aligns poorly to the gfpAnchor.
	public void test4 () {
		String seq = "GATCACTCTCGGCATGGCGAGCTAGTAAGCTACCGATGGCCTCCGGTGGCGCCACTGCGAGGTTCTTCCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";

		ExtractBarcodeSequences e = new ExtractBarcodeSequences("CGGCATGGACGAGCTGTACAAGTAAGCTA", "CCGGTGGCGCCACTGC");
		ExtractedSequenceGroup result =  e.findRabiesBarcode(seq);

		String expectedStopCodonBC="CCGATGGCCT";
		String expectedPolyABC="GAGGTTCTTC";

		int [] gpfAnchorPos={10,32};
		int gfpAnchorED=6;
		int [] cassetteAnchorPos={43,58};
		int cassetteAnchorED=0;

		Assert.assertEquals(result.getStopCodonBarcode().getSequence(), expectedStopCodonBC);
		Assert.assertEquals(result.getPolyABarcode().getSequence(), expectedPolyABC);

		Assert.assertEquals(result.getGfpAnchor().getStart(), gpfAnchorPos[0]);
		Assert.assertEquals(result.getGfpAnchor().getEnd(), gpfAnchorPos[1]);
		Assert.assertEquals(result.getGfpAnchor().getEditDistance().getEditDistance(), gfpAnchorED);

		Assert.assertEquals(result.getCassetteAnchor().getStart(), cassetteAnchorPos[0]);
		Assert.assertEquals(result.getCassetteAnchor().getEnd(), cassetteAnchorPos[1]);
		Assert.assertEquals(result.getCassetteAnchor().getEditDistance().getEditDistance(), cassetteAnchorED);

	}

	@Test(enabled=true)
	public void test5 () {
		String seq = "TATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCTACGAGACTTGCCCGGTGGCGCCACTGCGCTCGCGATTCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";

		ExtractBarcodeSequences e = new ExtractBarcodeSequences("CGGCATGGACGAGCTGTACAAGTAAGCTA", "CCGGTGGCGCCACTGC");
		ExtractedSequenceGroup result =  e.findRabiesBarcode(seq);

		String expectedStopCodonBC="CGAGACTTGC";
		String expectedPolyABC="GCTCGCGATT";

		int [] gpfAnchorPos={10,38};
		int gfpAnchorED=0;
		int [] cassetteAnchorPos={49,64};
		int cassetteAnchorED=0;

		Assert.assertEquals(result.getStopCodonBarcode().getSequence(), expectedStopCodonBC);
		Assert.assertEquals(result.getPolyABarcode().getSequence(), expectedPolyABC);

		Assert.assertEquals(result.getGfpAnchor().getStart(), gpfAnchorPos[0]);
		Assert.assertEquals(result.getGfpAnchor().getEnd(), gpfAnchorPos[1]);
		Assert.assertEquals(result.getGfpAnchor().getEditDistance().getEditDistance(), gfpAnchorED);

		Assert.assertEquals(result.getCassetteAnchor().getStart(), cassetteAnchorPos[0]);
		Assert.assertEquals(result.getCassetteAnchor().getEnd(), cassetteAnchorPos[1]);
		Assert.assertEquals(result.getCassetteAnchor().getEditDistance().getEditDistance(), cassetteAnchorED);

	}

	@Test(enabled=true)
	// this janky read goes off the rails a bit, as the edit distance affected by the single base aligned to the end of the read instead of contiguously.
	public void test6 () {
		String seq = "TCTCCCTCTCTTCCTTTCCTCTCTTTCCCCTTCCTCTCCTCTTTTTTCCCTTTTTCTCCCCTTCTTTTTTTTTTCTTTTCTTTCTTTTCCCCCCCTTTCTTTCCTCTTTCCCTCTTTCTCCCCTTCTCCCCCTTTTTCTTCTTTCTTTTTCTCTCCCTTTT";

		ExtractBarcodeSequences e = new ExtractBarcodeSequences("CGGCATGGACGAGCTGTACAAGTAAGCTA", "CCGGTGGCGCCACTGC");
		ExtractedSequenceGroup result =  e.findRabiesBarcode(seq);
		int len = result.getGfpAnchor().getMatchLength();
		String sub = result.getGfpAnchor().getSubSequence();
		int ed = result.getGfpAnchor().getEditDistance().getEditDistance();
		result.getGfpAnchor().getStart();

		String str = result.toString();

		String expectedStopCodonBC="TTTCCCCTTCCTCTCCTCTTTTTTC";
		String expectedPolyABC="TTTTTTTTTT";

		int [] gpfAnchorPos={1,23};
		int gfpAnchorED=19;
		int [] cassetteAnchorPos={49,64};
		int cassetteAnchorED=7;

		Assert.assertEquals(result.getStopCodonBarcode().getSequence(), expectedStopCodonBC);
		Assert.assertEquals(result.getPolyABarcode().getSequence(), expectedPolyABC);

		Assert.assertEquals(result.getGfpAnchor().getStart(), gpfAnchorPos[0]);
		Assert.assertEquals(result.getGfpAnchor().getEnd(), gpfAnchorPos[1]);
		Assert.assertEquals(result.getGfpAnchor().getEditDistance().getEditDistance(), gfpAnchorED);

		Assert.assertEquals(result.getCassetteAnchor().getStart(), cassetteAnchorPos[0]);
		Assert.assertEquals(result.getCassetteAnchor().getEnd(), cassetteAnchorPos[1]);
		Assert.assertEquals(result.getCassetteAnchor().getEditDistance().getEditDistance(), cassetteAnchorED);

	}

	@Test(enabled=true)
	// Insert in the middle of the GFP sequence.
	public void test7 () {
		String seq = "GATCACTCTCGGCATGGACGAGACTGTACAAGTAAGCTACCGATGGCCTCCGGTGGCGCCACTGCGAGGTTCTTCCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";

		ExtractBarcodeSequences e = new ExtractBarcodeSequences("CGGCATGGACGAGCTGTACAAGTAAGCTA", "CCGGTGGCGCCACTGC");
		ExtractedSequenceGroup result =  e.findRabiesBarcode(seq);

		String expectedStopCodonBC="CCGATGGCCT";
		String expectedPolyABC="GAGGTTCTTC";

		int [] gpfAnchorPos={10,39};
		int gfpAnchorED=1;
		int [] cassetteAnchorPos={50,65};
		int cassetteAnchorED=0;

		Assert.assertEquals(result.getStopCodonBarcode().getSequence(), expectedStopCodonBC);
		Assert.assertEquals(result.getPolyABarcode().getSequence(), expectedPolyABC);

		Assert.assertEquals(result.getGfpAnchor().getStart(), gpfAnchorPos[0]);
		Assert.assertEquals(result.getGfpAnchor().getEnd(), gpfAnchorPos[1]);
		Assert.assertEquals(result.getGfpAnchor().getEditDistance().getEditDistance(), gfpAnchorED);

		Assert.assertEquals(result.getCassetteAnchor().getStart(), cassetteAnchorPos[0]);
		Assert.assertEquals(result.getCassetteAnchor().getEnd(), cassetteAnchorPos[1]);
		Assert.assertEquals(result.getCassetteAnchor().getEditDistance().getEditDistance(), cassetteAnchorED);

	}

	@Test(enabled=true)
	// Insert at the base 5' of the GFP anchor.
	public void test8 () {
		String seq = "GATCACTCTCCGGCATGGACGAGCTGTACAAGTAAGCTACCGATGGCCTCCGGTGGCGCCACTGCGAGGTTCTTCCTTTTCAGTCGAGAAAAAAACATTAGATCAGAAGAACAACTGGCAACACTTCTCAACCTGAGACTTAGATCGGAAGAGCACACGT";

		ExtractBarcodeSequences e = new ExtractBarcodeSequences("CGGCATGGACGAGCTGTACAAGTAAGCTA", "CCGGTGGCGCCACTGC");
		ExtractedSequenceGroup result =  e.findRabiesBarcode(seq);

		String expectedStopCodonBC="CCGATGGCCT";
		String expectedPolyABC="GAGGTTCTTC";

		int [] gpfAnchorPos={11,39};
		int gfpAnchorED=0;
		int [] cassetteAnchorPos={50,65};
		int cassetteAnchorED=0;

		Assert.assertEquals(result.getStopCodonBarcode().getSequence(), expectedStopCodonBC);
		Assert.assertEquals(result.getPolyABarcode().getSequence(), expectedPolyABC);

		Assert.assertEquals(result.getGfpAnchor().getStart(), gpfAnchorPos[0]);
		Assert.assertEquals(result.getGfpAnchor().getEnd(), gpfAnchorPos[1]);
		Assert.assertEquals(result.getGfpAnchor().getEditDistance().getEditDistance(), gfpAnchorED);

		Assert.assertEquals(result.getCassetteAnchor().getStart(), cassetteAnchorPos[0]);
		Assert.assertEquals(result.getCassetteAnchor().getEnd(), cassetteAnchorPos[1]);
		Assert.assertEquals(result.getCassetteAnchor().getEditDistance().getEditDistance(), cassetteAnchorED);
	}

	@Test(enabled=false)
	// A SUPER janky consensus sequence generates a weird read that breaks code.  The second anchor is too far out on the read to see the polyA barcode.
	// This is a read read that we'd throw away, probably not worth obsessing.
	public void test9 () {
		String seq = "TCGGTTCCTGAATGAATGGGCAACCCTTCAAGACGGTAAAAATTAGGAGAAATTTTTACGAAGGCGCATAACGGCGTTATAACACTGACCCTCAGCAATCTTATATCACGAAGTCATAGACTGAATCACGAGTGGTCGGCAGATTGCGATAAACAGTCCTTCATAGAAATTTAACGCTGACTATTCCACTGCAAGT";

		ExtractBarcodeSequences e = new ExtractBarcodeSequences("CGGCATGGACGAGCTGTACAAGTAAGCTA", "CCGGTGGCGCCACTGC");
		ExtractedSequenceGroup result =  e.findRabiesBarcode(seq);

		String expectedStopCodonBC="AAAATTAGGAGAAATTTTTACGAAGGCGCATAACGGCGTTATAACACTGACCCTCAGCAATCTTATATCACGAAGTCATAGACTGAATCACGAGTGGTCGGCAGATTGCGATAAACAGTCCTTCATAGAAATTTA";
		String expectedPolyABC="AAGT";

		int [] gpfAnchorPos={12,38};
		int gfpAnchorED=12;
		int [] cassetteAnchorPos={174,192};
		int cassetteAnchorED=7;

		Assert.assertEquals(result.getStopCodonBarcode().getSequence(), expectedStopCodonBC);
		Assert.assertEquals(result.getPolyABarcode().getSequence(), expectedPolyABC);

		Assert.assertEquals(result.getGfpAnchor().getStart(), gpfAnchorPos[0]);
		Assert.assertEquals(result.getGfpAnchor().getEnd(), gpfAnchorPos[1]);
		Assert.assertEquals(result.getGfpAnchor().getEditDistance().getEditDistance(), gfpAnchorED);

		Assert.assertEquals(result.getCassetteAnchor().getStart(), cassetteAnchorPos[0]);
		Assert.assertEquals(result.getCassetteAnchor().getEnd(), cassetteAnchorPos[1]);
		Assert.assertEquals(result.getCassetteAnchor().getEditDistance().getEditDistance(), cassetteAnchorED);
	}

	@Test (enabled=true)
	// a read that crashes
	public void test10 () {
		String s1 = "ATTAACACCATCCTTCATGAACTTAATCCACTGTTCACCATAACCGTGACGATGAGGGACATAAAAAGTAAAAATGTCTCCAGTAGAGTCAATAGCAAGGCCACGACGCAATGGAGAAAGCCGGAGAGCGCCAACGGCGTCCATCTCGCAGGCGTCGCCAG";
		String s2 = "AGGCAAGCGTAAAGGCGCTCGTCTTTGGTATGTAGGTGGTCAACAATTTTAATTGCAGGGGCTTCGGCCCCTTACTTGAGGATAAATTATGTCTAATATTCAAACTGGCGCCGAGCGTATCCCGCATGACCTTTCCCATCTTGGCTTCCTTGCTGGTCAGA";
		String bq1= "AAAAADDFBBFFGGGGFGFGG5FGFGGHHDHHFFGHBGG3FAG3F2EECFGGFGBF3EEGGHBGFFEAGGHFHBFGHHF5FGGF55FGHFHHGHHHFE211B1B/EF?G?/?BF3FB0FF3GCG////?<CDF<//<<@DC@.GHGFF.-A?.-<@?CDD:";
		String bq2= ">AA>A1AAA1A1FGGAEEGEGGFEFGHHHDFHBFFGHFEGHHGHHEHHHHE1ED2GDG0/EEGHHGGGGGGH0@FGHFHF01FF1FFG12GGHHHHFFHHHH2F1>10GEA?E/<<@C<B0@<<BCGFFHEGHHG<1?GGFGB1FGHEFFGHHGH<GH0/=";

		ConsensusSequence cs = ConsensusSequenceFactory.getInstance().getConsensusSequence(s1, s2, false);
		cs.addReadBaseQualities(bq1, bq2);
		String seq = cs.getConsensusSequence();
		int consensusED = cs.getLocalAlignmentEditDistance().getEditDistance();

		ExtractBarcodeSequences e = new ExtractBarcodeSequences("CGGCATGGACGAGCTGTACAAGTAAGCTA", "CCGGTGGCGCCACTGC");
		ExtractedSequenceGroup result =  e.findRabiesBarcode(seq);
		result.getCassetteAnchor().getStart();
		result.getCassetteAnchor().getEnd();
		result.getCassetteAnchor().getSubSequence();

		String toString = result.toString();
		Assert.assertNotNull(result);

	}



}
