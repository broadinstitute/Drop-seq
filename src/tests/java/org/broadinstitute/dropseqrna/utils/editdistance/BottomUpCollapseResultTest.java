/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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

import java.util.List;
import java.util.Set;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.google.common.collect.ImmutableList;

public class BottomUpCollapseResultTest {

	@Test
	public void testFilterToCommonModes () {
		BottomUpCollapseResult c = new BottomUpCollapseResult();
		// A->T base 1 (position 0).
		c.addPair("TCGTA", "ACGTA");
		c.addPair("TCTTT", "ACTTT");
		c.addPair("TCTAA", "ACTAA");
		c.addPair("TGTAA", "AGTAA");
		// A->C position 0. (should be ambiguous)
		c.addPair("CCGTA","ACGTA");



		// G->A base 3 (position 2)
		c.addPair("ACATA","ACGTA");
		c.addPair("TGATA","TGGTA");
		c.addPair("ATATA","ATGTA");
		c.addPair("TCATA","TCGTA");

		// add a little noise (G->C) (should be ambiguous)
		c.addPair("TCCTA","TCGTA");


		// add in some noise at other positions.
		// position 1 changes, expected freq =0.33.  All should be ambiguous
		c.addPair("GATCT", "GGTCT");
		c.addPair("GCTCT", "GGTCT");
		c.addPair("GTTCT", "GGTCT");

		// position 3 changes, expected freq =0.33.  All should be ambiguous
		c.addPair("CTGAA", "CTGTA");
		c.addPair("CTGCA", "CTGTA" );
		c.addPair("CTGGA", "CTGTA");


		BottomUpCollapseResult result = c.makeNonCommonChangesAmbiguous(0.5);

		// there should be lots of ambiguous barcodes I didn't add directly.
		List<String> expectedAmbiguous = ImmutableList.of("CCGTA", "TCCTA", "GATCT", "GCTCT", "GTTCT", "CTGAA", "CTGCA", "CTGGA");
		Set<String> ambiguousBC = result.getAmbiguousBarcodes();

		for (String ea : expectedAmbiguous)
			Assert.assertTrue(ambiguousBC.contains(ea));

		// test that the smaller->larger barcodes are preserved for the expected patterns, but not the noise
		Assert.assertEquals(result.getLargerRelatedBarcode("TCGTA"), "ACGTA");
		Assert.assertEquals(result.getLargerRelatedBarcode("TCTTT"), "ACTTT");
		Assert.assertEquals(result.getLargerRelatedBarcode("TCTAA"), "ACTAA");
		Assert.assertEquals(result.getLargerRelatedBarcode("TGTAA"), "AGTAA");
		Assert.assertNull(result.getLargerRelatedBarcode("CCGTA"));

		Assert.assertEquals(result.getLargerRelatedBarcode("ACATA"), "ACGTA");
		Assert.assertEquals(result.getLargerRelatedBarcode("TGATA"), "TGGTA");
		Assert.assertEquals(result.getLargerRelatedBarcode("ATATA"), "ATGTA");
		Assert.assertEquals(result.getLargerRelatedBarcode("TCATA"), "TCGTA");
		Assert.assertNull(result.getLargerRelatedBarcode("TCCTA"));

	}
}
