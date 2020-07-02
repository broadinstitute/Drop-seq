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
package org.broadinstitute.dropseqrna.utils.statistics;

import junit.framework.Assert;

import org.testng.annotations.Test;

public class DiversityTest {

	@Test
	//
	public void testDiversity () {
		// diversity(c(22, 69, 72, 26, 41, 13, 63, 83, 21, 46))
		// 2.161481
		int [] intData = {22, 69, 72, 26, 41, 13, 63, 83, 21, 46};
		double r1 = Diversity.diversity(intData);
		Assert.assertEquals(2.161481d, r1, 0.00001);

		// diversity (c(0.596662266, 0.047729539, 0.014175908, 0.070480358, 0.002488427, 0.600280966, 0.325990001, 0.204724048, 0.213749792, 0.429259571))
		// [1] 1.878304
		double [] doubleData = {0.596662266, 0.047729539, 0.014175908, 0.070480358, 0.002488427, 0.600280966, 0.325990001, 0.204724048, 0.213749792, 0.429259571};
		double r2 = Diversity.diversity(doubleData);
		Assert.assertEquals(1.878304d, r2, 0.00001);
	}

	@Test
	public void testEquitability () {
		// diversity(c(22, 69, 72, 26, 41, 13, 63, 83, 21, 46))/log(10)
		// 0.9387194
		int [] intData = {22, 69, 72, 26, 41, 13, 63, 83, 21, 46};
		double r1 = Diversity.equitability(intData);
		Assert.assertEquals(0.9387194d, r1, 0.00001);

		// diversity (c(0.596662266, 0.047729539, 0.014175908, 0.070480358, 0.002488427, 0.600280966, 0.325990001, 0.204724048, 0.213749792, 0.429259571))/log(10)
		// [1] 0.8157371
		double [] doubleData = {0.596662266, 0.047729539, 0.014175908, 0.070480358, 0.002488427, 0.600280966, 0.325990001, 0.204724048, 0.213749792, 0.429259571};
		double r2 = Diversity.equitability(doubleData);
		Assert.assertEquals(0.8157371d, r2, 0.00001);
	}
}
