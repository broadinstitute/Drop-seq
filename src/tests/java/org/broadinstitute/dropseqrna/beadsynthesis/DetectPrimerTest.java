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
