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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileup;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.SampleGenotypeProbabilities;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import picard.annotation.LocusFunction;

public class SampleGenotypeProbabilitiesTest {

	// Test against R implementation
	// dropseqrna/transcriptome/R/DAC/SampleCellBarcodeAssignment.R

	// See smallTest_snpUMIPileUp.summary.txt for a listing of the
	// bases/qualities/cell/molecular barcodes for all 15 reads
	
	private final static String rootDir="testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts";
	
	private final File smallBAMFile = new File(rootDir+"/smallTest_snpUMIPileUp.retagged.bam");

	@Test(enabled = true)

	public void testSummarziedBQFromFile() {

		String cellBC = "ATCAGGGACAGA";
		String molBC = "CGGGGCTC";
		/*
		 * Read name =NS500217:67:H14GMBGXX:1:22104:14231:13073 Base = G Base phred
		 * quality = 37 XC = ATCAGGGACAGA XM = CGGGGCTC Read name
		 * =NS500217:67:H14GMBGXX:1:13109:1238:14496 Base = G Base phred quality = 37 XC
		 * = ATCAGGGACAGA XM = CGGGGCTC Read name
		 * =NS500217:67:H14GMBGXX:1:13110:26687:14501 Base = G Base phred quality = 27
		 * XC = ATCAGGGACAGA XM = CGGGGCTC
		 */

		int snpPos = 76227022;
		Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
		SampleGenotypeProbabilities result = new SampleGenotypeProbabilities(snpInterval, cellBC);
		SNPUMIBasePileup p = getPileUpFromFile(cellBC, molBC);
		result.add(p);
		List<Byte> quals = result.getQualities();
		Byte summarized = quals.get(0);
		Assert.assertTrue((byte) 31 == summarized);

		// let's add another pileup in that's within ed=1 and incorporates an error.
		String molBC2 = "AGGGGCTC";
		SNPUMIBasePileup p2 = getPileUpFromFile(cellBC, molBC2);
		char[] b2 = { 'G', 'G', 'T' };
		byte[] bases2 = new byte[b2.length];
		StringUtil.charsToBytes(b2, 0, b2.length, bases2, 0);
		byte[] quals2 = { 37, 37, 10 };
		p2.setBasesAndQualities(bases2, quals2);
		result.add(p2);
		result.collapseUMIs(1);
		List<Byte> qualResult2 = result.getQualities();
		Assert.assertEquals(qualResult2.size(), 1);
		Assert.assertEquals(Byte.valueOf((byte) 8), qualResult2.get(0));

	}

	@Test
	public void testCollapseByEditDistance() {
		/*
		 * Read name =NS500217:67:H14GMBGXX:1:22104:14231:13073 Base = G Base phred
		 * quality = 37 XC = ATCAGGGACAGA XM = CGGGGCTC Read name
		 * =NS500217:67:H14GMBGXX:1:13109:1238:14496 Base = G Base phred quality = 37 XC
		 * = ATCAGGGACAGA XM = CGGGGCTC Read name
		 * =NS500217:67:H14GMBGXX:1:13110:26687:14501 Base = G Base phred quality = 27
		 * XC = ATCAGGGACAGA XM = CGGGGCTC
		 */
		String cellBC = "ATCAGGGACAGA";
		String molBC = "CGGGGCTC";
		int snpPos = 76227022;
		Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
		SampleGenotypeProbabilities result = new SampleGenotypeProbabilities(snpInterval, cellBC);
		SNPUMIBasePileup p = getPileUpFromFile(cellBC, molBC);
		result.add(p);

		// let's add another pileup in that's within ed=1 and incorporates an error.
		// result has only a single base because two UMIs are merged.
		String molBC2 = "AGGGGCTC";
		SNPUMIBasePileup p2 = getPileUpFromFile(cellBC, molBC2);
		char[] b2 = { 'G', 'G', 'T' };
		byte[] bases2 = new byte[b2.length];
		StringUtil.charsToBytes(b2, 0, b2.length, bases2, 0);
		byte[] quals2 = { 37, 37, 10 };
		p2.setBasesAndQualities(bases2, quals2);
		result.add(p2);
		result.collapseUMIs(1);
		List<Byte> qualResult2 = result.getQualities();
		Assert.assertEquals(qualResult2.size(), 1);
		Assert.assertEquals(Byte.valueOf((byte) 8), qualResult2.get(0));
	}

	@Test
	public void testNoCollapseByEditDistance() {
		/*
		 * Read name =NS500217:67:H14GMBGXX:1:22104:14231:13073 Base = G Base phred
		 * quality = 37 XC = ATCAGGGACAGA XM = CGGGGCTC Read name
		 * =NS500217:67:H14GMBGXX:1:13109:1238:14496 Base = G Base phred quality = 37 XC
		 * = ATCAGGGACAGA XM = CGGGGCTC Read name
		 * =NS500217:67:H14GMBGXX:1:13110:26687:14501 Base = G Base phred quality = 27
		 * XC = ATCAGGGACAGA XM = CGGGGCTC
		 */
		String cellBC = "ATCAGGGACAGA";
		String molBC = "CGGGGCTC";
		int snpPos = 76227022;
		Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
		SampleGenotypeProbabilities result = new SampleGenotypeProbabilities(snpInterval, cellBC);
		SNPUMIBasePileup p = getPileUpFromFile(cellBC, molBC);
		result.add(p);

		// let's add another pileup in nowhere near the ED
		// result has 2 different bases for the two pileups.
		String molBC2 = "ZZZZZZZZ";
		SNPUMIBasePileup p2 = getPileUpFromFile(cellBC, molBC2);
		char[] b2 = { 'G', 'G', 'T' };
		byte[] bases2 = new byte[b2.length];
		StringUtil.charsToBytes(b2, 0, b2.length, bases2, 0);
		byte[] quals2 = { 37, 37, 10 };
		p2.setBasesAndQualities(bases2, quals2);
		result.add(p2);
		result.collapseUMIs(1);

		List<Byte> qualResult2 = result.getQualities();
		Assert.assertEquals(qualResult2.size(), 2);
		Assert.assertEquals(Byte.valueOf((byte) 31), qualResult2.get(0));
		Assert.assertEquals(Byte.valueOf((byte) 5), qualResult2.get(1));
	}

	/**
	 * Data added. Read name =NS500217:67:H14GMBGXX:1:22107:20436:18875 Base = NA
	 * Base phred quality = NA XC = ATCAGGGACAGA XM = CGCGGGAC // this read doesn't
	 * overlap the SNP so it isn't added!
	 *
	 * Read name =NS500217:67:H14GMBGXX:1:22104:14231:13073 Base = G Base phred
	 * quality = 37 XC = ATCAGGGACAGA XM = CGGGGCTC Read name
	 * =NS500217:67:H14GMBGXX:1:13109:1238:14496 Base = G Base phred quality = 37 XC
	 * = ATCAGGGACAGA XM = CGGGGCTC Read name
	 * =NS500217:67:H14GMBGXX:1:13110:26687:14501 Base = G Base phred quality = 27
	 * XC = ATCAGGGACAGA XM = CGGGGCTC
	 *
	 * Read name =NS500217:67:H14GMBGXX:3:23612:5963:17486 Base = G Base phred
	 * quality = 37 XC = ATCAGGGACAGA XM = CTGGGCTC
	 *
	 * Read name =NS500217:67:H14GMBGXX:4:13611:8735:3829 Base = A Base phred
	 * quality = 8 XC = ATCAGGGACAGA XM = GGAATGTG Read name
	 * =NS500217:67:H14GMBGXX:1:13202:10555:15929 Base = A Base phred quality = 8 XC
	 * = ATCAGGGACAGA XM = GTGCCGTG Read name
	 * =NS500217:67:H14GMBGXX:1:21110:17335:17576 Base = G Base phred quality = 37
	 * XC = ATCAGGGACAGA XM = TCGCAGAC
	 *
	 */
	@Test
	public void testGetCounts() {
		String cellBC = "ATCAGGGACAGA";

		int snpPos = 76227022;
		Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
		SampleGenotypeProbabilities result = new SampleGenotypeProbabilities(snpInterval, cellBC);
		SNPUMIBasePileup p = getPileUpFromFile(cellBC, "CGGGGCTC");
		p.addLocusFunction(LocusFunction.CODING);
		result.add(p);
		result.add(getPileUpFromFile(cellBC, "CGCGGGAC")); // this read doesn't overlap the SNP so it isn't added!
		result.add(getPileUpFromFile(cellBC, "CTGGGCTC"));
		result.add(getPileUpFromFile(cellBC, "GGAATGTG"));
		result.add(getPileUpFromFile(cellBC, "GTGCCGTG"));
		result.add(getPileUpFromFile(cellBC, "TCGCAGAC"));

		// the code coverage gods are pleased.
		Assert.assertEquals(result.getCell(), cellBC);
		Assert.assertEquals(result.getSNPInterval(), snpInterval);

		Assert.assertNotNull(result.toString());

		// get the unfiltered results.
		ObjectCounter<Character> umiCounts = result.getUMIBaseCounts();
		ObjectCounter<Character> umiCountsExpected = new ObjectCounter<>();
		umiCountsExpected.incrementByCount('A', 2);
		umiCountsExpected.incrementByCount('G', 3);
		Assert.assertEquals(umiCounts, umiCountsExpected);

		ObjectCounter<Character> readCounts = result.getReadBaseCounts();
		ObjectCounter<Character> readCountsExpected = new ObjectCounter<>();
		readCountsExpected.incrementByCount('A', 2);
		readCountsExpected.incrementByCount('G', 5);
		Assert.assertEquals(readCounts, readCountsExpected);

		Set<LocusFunction> expected = new HashSet<>(Arrays.asList(LocusFunction.CODING));
		Assert.assertEquals(result.getLocusFunctions(), expected);

	}

	@Test(expectedExceptions=IllegalArgumentException.class)
	public void testIllegalToAddPileupInterval() {
		String cellBC = "ATCAGGGACAGA";
		int snpPos = 76227022;
		Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
		Interval snpInterval2 = new Interval("HUMAN_2", snpPos, snpPos, true, "test2");
		SampleGenotypeProbabilities result = new SampleGenotypeProbabilities(snpInterval, cellBC);
		result.add(getPileUpFromFile(cellBC, "CTGGGCTC"));
		SNPUMIBasePileup p = new SNPUMIBasePileup(snpInterval2,"TESTGENE", cellBC, "AAAAAAA");
		p.addBaseAndQuality((byte)1, (byte) 1);
		result.add(p);
	}

	@Test(expectedExceptions=IllegalArgumentException.class)
	public void testIllegalToAddPileupCellBC() {
		String cellBC = "ATCAGGGACAGA";
		int snpPos = 76227022;
		Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
		SampleGenotypeProbabilities result = new SampleGenotypeProbabilities(snpInterval, cellBC);
		result.add(getPileUpFromFile(cellBC, "CTGGGCTC"));
		SNPUMIBasePileup p = new SNPUMIBasePileup(snpInterval,"TESTGENE", "ZZZZZZ", "AAAAAAA");
		p.addBaseAndQuality((byte)1, (byte) 1);
		result.add(p);
	}

	@Test
	// See LikelihoodUtilsTest for more in-depth testing, since this class defers to
	// that one. 
	public void testGetLikelihoods() {
		String cellBC = "ATCAGGGACAGA";

		int snpPos = 76227022;
		Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
		SampleGenotypeProbabilities result = new SampleGenotypeProbabilities(snpInterval, cellBC);
		result.add(getPileUpFromFile(cellBC, "CGCGGGAC")); // this read doesn't overlap the SNP so it isn't added!
		result.add(getPileUpFromFile(cellBC, "CGGGGCTC"));
		result.add(getPileUpFromFile(cellBC, "CTGGGCTC"));
		result.add(getPileUpFromFile(cellBC, "GGAATGTG"));
		result.add(getPileUpFromFile(cellBC, "GTGCCGTG"));
		result.add(getPileUpFromFile(cellBC, "TCGCAGAC"));

		// the A bases are low confidence, so the GG likelihood is closer to the GA
		// likelihood even though this is a 3G2A SNP pileup.
		double likeHet = result.getLogLikelihood('G', 'A', null, null, null, null, null, null);
		double likeHetExpected = -1.505;
		Assert.assertEquals(likeHet, likeHetExpected, 0.001);
		double likeRef = result.getLogLikelihood('G', 'G', null, null, null,null, null, null);
		double likeRefExpected = -1.600;
		Assert.assertEquals(likeRef, likeRefExpected, 0.001);
		double likeAlt = result.getLogLikelihood('A', 'A', null, null, null,null, null, null);
		double likeAltExpected = -10.6498;
		Assert.assertEquals(likeAlt, likeAltExpected, 0.001);

		// test with a fixed error rate. More testing takes place in LikelihoodUtilsTest
		double likeAltFixed = result.getLogLikelihood('A', 'A', 1e-1, null, null,null, null, null);
		double likeAltExpectedFixed = -3.0915;
		Assert.assertEquals(likeAltFixed, likeAltExpectedFixed, 0.001);

		// test with a capped error rate. More testing takes place in
		// LikelihoodUtilsTest
		double likeAltCapped = result.getLogLikelihood('A', 'A', null, null, 1e-2,null, null, null);
		double likeAltExpectedCapped = -6.149;
		Assert.assertEquals(likeAltCapped, likeAltExpectedCapped, 0.001);

	}

	private SNPUMIBasePileup getPileUpFromFile(final String cellBC, final String molBC) {
		SamReader reader = SamReaderFactory.makeDefault().open(smallBAMFile);
		int snpPos = 76227022;
		Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
		SNPUMIBasePileup p = new SNPUMIBasePileup(snpInterval, "ACADM", cellBC, molBC);
		// cheap and cheerful way to not use the Iterator (yet).
		for (SAMRecord r : reader) {
			String currentCell = r.getStringAttribute("XC");
			String currentMolBC = r.getStringAttribute("XM");
			if (currentCell.equals(cellBC) && currentMolBC.equals(molBC))
				p.addRead(r);
		}

		return p;

	}

}
