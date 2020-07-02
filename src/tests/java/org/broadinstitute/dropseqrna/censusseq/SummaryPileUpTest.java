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
package org.broadinstitute.dropseqrna.censusseq;

import org.broadinstitute.dropseqrna.censusseq.SummaryPileUp;
import org.testng.Assert;
import org.testng.annotations.Test;


public class SummaryPileUpTest {

	@Test
	public void testOne () {
		SummaryPileUp p = new SummaryPileUp("TEST");
		p.incrementAltCount(5932);
		p.incrementRefCount(70519);
		double result = p.getRatioByAltAlleleFraction();
		Assert.assertEquals(result, 0.1551844, 0.0001);
	}

	@Test
	public void testBoringStuff() {
		SummaryPileUp p = new SummaryPileUp("TEST");
		p.incrementAltCount(5932);
		p.incrementRefCount(70519);
		p.incrementNumSNPs();

		Assert.assertEquals(p.getDonor(), "TEST");
		Assert.assertEquals(p.getAltCount(), 5932);
		Assert.assertEquals(p.getRefCount(), 70519);
		Assert.assertEquals(p.getNumSNPs(), 1);
		Assert.assertNotNull(p.toString());


	}
}
