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

import org.broadinstitute.dropseqrna.utils.editdistance.BarcodeSubstitutionCollection.BarcodeSubstitutionElement;
import org.testng.Assert;
import org.testng.annotations.Test;

public class BarcodeSubstitutionCollectionTest {

	@Test
	public void testBiasBase1And3 () {
		BarcodeSubstitutionCollection c = new BarcodeSubstitutionCollection();
		// A->T base 1 (position 0).
		for (int i=0; i<5; i++) {
			String intended = "ACGTA";
			String neighbor = "TCGTA";
			c.add(intended, neighbor);
		}
		// G->A base 3 (position 2)
		for (int i=0; i<5; i++) {
			String intended = "ACGTA";
			String neighbor = "ACATA";
			c.add(intended, neighbor);
		}

		List<Integer> positions = c.getPositions();
		Assert.assertEquals(positions.size(), 2);
		int pos = 0;

		BarcodeSubstitutionElement e = c.getMostCommonSubsitution(pos);
		double freq = c.getSubsitutionFrequency(e, pos);
		Assert.assertTrue(e.getIntendedBase().equals("A") && e.getNeighborbase().equals("T"));
		Assert.assertEquals(freq, 1, 0.001);

		pos=2;
		e = c.getMostCommonSubsitution(pos);
		freq = c.getSubsitutionFrequency(e, pos);
		Assert.assertTrue(e.getIntendedBase().equals("G") && e.getNeighborbase().equals("A") );
		Assert.assertEquals(freq, 1, 0.001);

		// for another position, there is no dominant base
		pos=3;
		e = c.getMostCommonSubsitution(pos);
		Assert.assertNull(e);
		freq = c.getSubsitutionFrequency(e, pos);
		Assert.assertEquals(freq, 0, 0.001);
	}

	@Test
	public void testFilterToCommonSubstitutionPatterns() {
		BarcodeSubstitutionCollection c = new BarcodeSubstitutionCollection();
		// A->T base 1 (position 0).
		for (int i=0; i<5; i++) {
			String intended = "ACGTA";
			String neighbor = "TCGTA";
			c.add(intended, neighbor);
		}
		// add a little noise
		c.add("ACGTA", "CCGTA");

		// G->A base 3 (position 2)
		for (int i=0; i<5; i++) {
			String intended = "ACGTA";
			String neighbor = "ACATA";
			c.add(intended, neighbor);
		}
		// add a little noise
		c.add("ACGTA", "ACCTA");

		//validate that the freqs are not one because of the noise.
		int pos = 0;

		BarcodeSubstitutionElement e = c.getMostCommonSubsitution(pos);
		double freq = c.getSubsitutionFrequency(e, pos);
		Assert.assertTrue(e.getIntendedBase().equals("A") && e.getNeighborbase().equals("T"));
		Assert.assertEquals(freq, 0.8333333, 0.001);

		pos=2;
		e = c.getMostCommonSubsitution(pos);
		freq = c.getSubsitutionFrequency(e, pos);
		Assert.assertTrue(e.getIntendedBase().equals("G") && e.getNeighborbase().equals("A") );
		Assert.assertEquals(freq, 0.8333333, 0.001);

		// add in some noise at other positions.
		// position 1 changes, expected freq =0.33
		c.add("GGTCT", "GATCT");
		c.add("GGTCT", "GCTCT");
		c.add("GGTCT", "GTTCT");

		// position 3 changes, expected freq =0.33
		c.add("CTGTA", "CTGAA");
		c.add("CTGTA", "CTGCA");
		c.add("CTGTA", "CTGGA");

		BarcodeSubstitutionCollection result= c.filterToCommonSubstitutionPatterns(0.5);
		// there should only be 2 positions with >50% the same substitution.
		List<Integer> positions = result.getPositions();
		Assert.assertEquals(positions.size(), 2);

		// after cleanup!
		pos = 0;

		e = result.getMostCommonSubsitution(pos);
		freq = result.getSubsitutionFrequency(e, pos);
		Assert.assertTrue(e.getIntendedBase().equals("A") && e.getNeighborbase().equals("T"));
		Assert.assertEquals(freq, 1, 0.001);

		pos=2;
		e = result.getMostCommonSubsitution(pos);
		freq = result.getSubsitutionFrequency(e, pos);
		Assert.assertTrue(e.getIntendedBase().equals("G") && e.getNeighborbase().equals("A") );
		Assert.assertEquals(freq, 1, 0.001);

	}
}
