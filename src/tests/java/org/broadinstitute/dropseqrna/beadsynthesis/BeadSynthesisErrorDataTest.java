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
package org.broadinstitute.dropseqrna.beadsynthesis;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Assert;
import org.testng.annotations.Test;

public class BeadSynthesisErrorDataTest {
	// from gunzip -c Spermidine3_1000.molBC.txt.gz |grep ATACAGTCTACA |cut -f 3
	// |head -n 25
	private final String[] umis = { "GCCTCACT", "CGTCTCAT", "TCGCTTTT",
			"TATGACTT", "TGATGTTT", "GCTTAATT", "TTAATGTT", "CCGCGGGT",
			"TGCCTTGT", "TATGCTGT", "CAGACGCT", "GAAGATTT", "CCTCGTCT",
			"TATTAAGT", "GTGTTTGT", "TCATATCT", "ATTGGAAT", "AGCCGGTT",
			"TCCTTTTT", "ATCCGTCT", "GGGGTGCT", "GATGTTTT", "GGGCGTAT",
			"TAGCTCGT", "GATCTATT" };

	@Test
	public void test1() {
		double threshold = 0.9;

		BeadSynthesisErrorData d = new BeadSynthesisErrorData("ATACAGTCTACA");
		// to test adding a collection of UMIs.
		List<String> umiCollection = Arrays.stream(umis).collect(Collectors.toList());
		d.addUMI(umiCollection);

		d.getErrorType(threshold);
		d.finalize();
		boolean hasError = d.hasSynthesisError(threshold);
		Assert.assertTrue(hasError);

		int errorPositionExpected = 8;
		int errorPosition = d.getErrorBase(threshold);
		Assert.assertEquals(errorPositionExpected, errorPosition);

		// base 1: 12 16 32 40
		// base 2: 32 28 24 16
		// base 3: 16 20 28 36
		// base 4: 8 40 24 28
		// base 5: 20 12 28 40
		// base 6: 20 12 20 48
		// base 7: 12 24 24 40
		// base 8: 0 0 0 100
		double[] metric = d.synthesisErrorMetric();
		Assert.assertEquals(0.4, metric[0], 0.001);
		Assert.assertEquals(0.32, metric[1], 0.001);
		Assert.assertEquals(0.36, metric[2], 0.001);
		Assert.assertEquals(0.40, metric[3], 0.001);
		Assert.assertEquals(0.40, metric[4], 0.001);
		Assert.assertEquals(0.48, metric[5], 0.001);
		Assert.assertEquals(0.40, metric[6], 0.001);
		Assert.assertEquals(1.0, metric[7], 0.001);

		int errorBase = d.getPolyTErrorPosition(0.9);
		Assert.assertEquals(8, errorBase);
		// for debugging to hit cache.
		errorBase = d.getPolyTErrorPosition(0.9);
		Assert.assertEquals(8, errorBase);

		String s = d.toString();
		Assert.assertNotNull(s);

		Assert.assertTrue(d.hasSynthesisError(threshold));

		boolean biasedBase8 = d.isPolyTBiasedAtPosition(8,  threshold);
		Assert.assertTrue(biasedBase8);

		boolean biasedBase7 = d.isPolyTBiasedAtPosition(7,  threshold);
		Assert.assertFalse(biasedBase7);

		double lastBaseBias = d.getPolyTFrequencyLastBase();
		Assert.assertEquals(1, lastBaseBias, 0.01);

		d.incrementReads(this.umis.length);
		int numReads = d.getNumReads();
		Assert.assertEquals(this.umis.length, numReads);

	}

	@Test(expectedExceptions=IllegalArgumentException.class)
	public void testOutOfBounds () {
		BeadSynthesisErrorData d = new BeadSynthesisErrorData("ATACAGTCTACA");
		// to test adding a collection of UMIs.
		List<String> umiCollection = Arrays.stream(umis).collect(Collectors.toList());
		d.addUMI(umiCollection);
		d.isPolyTBiasedAtPosition(12, 0.9);
	}

	@Test
	public void testSingleUMIError () {
		double threshold = 0.9;

		BeadSynthesisErrorData d = new BeadSynthesisErrorData("ATACAGTCTACA");
		for (String umi : this.umis)
			d.addUMI(umis[0]);
		BeadSynthesisErrorType t= d.getErrorType(threshold);
		Assert.assertEquals(BeadSynthesisErrorType.SINGLE_UMI, t);


	}

	@Test
	public void testFixedFirstBase () {
		double threshold = 0.9;

		BeadSynthesisErrorData d = new BeadSynthesisErrorData("ATACAGTCTACA");
		for (String umi : this.umis) {
			// fix the first base.
			String newUMI = "A"+ umi.substring(0,7);
			d.addUMI(newUMI);
		}
		BeadSynthesisErrorType t= d.getErrorType(threshold);
		Assert.assertEquals(BeadSynthesisErrorType.FIXED_FIRST_BASE, t);
	}

	@Test
	public void testOtherError () {
		double threshold = 0.9;

		BeadSynthesisErrorData d = new BeadSynthesisErrorData("ATACAGTCTACA");
		for (String umi : this.umis) {
			// fix the first base, last base already fixed.
			String newUMI = "A"+ umi;
			d.addUMI(newUMI);
		}
		BeadSynthesisErrorType t= d.getErrorType(threshold);
		Assert.assertEquals(BeadSynthesisErrorType.OTHER_ERROR, t);
	}

	@Test
	public void testNoError () {
		double threshold = 0.9;

		BeadSynthesisErrorData d = new BeadSynthesisErrorData("ATACAGTCTACA");
		for (String umi : this.umis) {
			// fix the first base, last base already fixed.
			String newUMI = umi.substring(0,7);
			d.addUMI(newUMI);
		}
		BeadSynthesisErrorType t= d.getErrorType(threshold);
		Assert.assertEquals(BeadSynthesisErrorType.NO_ERROR, t);

	}

	@Test
	public void testDetectPrimer () {
		double threshold = 0.9;
		String UMI = this.umis[0];
		String primer = UMI;
		DetectPrimerInUMI dp = new DetectPrimerInUMI(primer);
		BeadSynthesisErrorData d = new BeadSynthesisErrorData("ATACAGTCTACA");
		d.addUMI(UMI);
		BeadSynthesisErrorType t= d.getErrorType(threshold, dp, new Integer (1));
		Assert.assertEquals(BeadSynthesisErrorType.PRIMER, t);

		d = new BeadSynthesisErrorData("ATACAGTCTACA");
		for (String umi : this.umis) {
			// fix the first base, last base already fixed.
			String newUMI = umi.substring(0,7);
			d.addUMI(newUMI);
		}
		t= d.getErrorType(threshold, dp, new Integer (1));
		Assert.assertEquals(BeadSynthesisErrorType.NO_ERROR, t);

	}








}
