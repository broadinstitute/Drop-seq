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

import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

public class DetectPrimerTest {

	private final String primer = "AAGCAGTGGTATCAACGCAGAGT";
	@Test (enabled=true)
	public void testOne() {
		DetectPrimerInUMI dpu = new DetectPrimerInUMI(this.primer);
		String test = "AAGCAGTG";
		boolean flag = dpu.isStringInPrimer(test, 1);
		Assert.assertTrue(flag);
	}
	
	@Test (enabled=true)
	public void testGetSubstrings() {
		DetectPrimerInUMI dpu = new DetectPrimerInUMI(this.primer);
		List<String> substrings = dpu.getSubstrings(8);
		String [] expected = {"AAGCAGTG","AGCAGTGG","GCAGTGGT","CAGTGGTA","AGTGGTAT","GTGGTATC","TGGTATCA","GGTATCAA","GTATCAAC","TATCAACG","ATCAACGC","TCAACGCA","CAACGCAG","AACGCAGA","ACGCAGAG","CGCAGAGT"};
		Assert.assertTrue(substrings.size()==expected.length);
		for (int i=0; i<expected.length; i++) {
			Assert.assertEquals(substrings.get(i), expected[i]);
		}		
	}
}
