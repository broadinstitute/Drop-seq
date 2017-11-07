/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.beadsynthesis;

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
		for (String umi : umis) {
			d.addUMI(umi);
		}
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
	}
}
