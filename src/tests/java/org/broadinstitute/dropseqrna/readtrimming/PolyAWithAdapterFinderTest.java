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
package org.broadinstitute.dropseqrna.readtrimming;

import org.broadinstitute.dropseqrna.readtrimming.AdapterDescriptor;
import org.broadinstitute.dropseqrna.readtrimming.PolyAFinder;
import org.broadinstitute.dropseqrna.readtrimming.PolyAWithAdapterFinder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.util.ClippingUtility;

public class PolyAWithAdapterFinderTest {


    @Test(dataProvider = "testWithAdapterTestCases")
    public void testWithAdapter(final String testName, final String read, final String adapter, final int expectedResult) {
        PolyAWithAdapterFinder polyAFinder = new PolyAWithAdapterFinder(
                new AdapterDescriptor(AdapterDescriptor.DEFAULT_ADAPTER),
                4,
                ClippingUtility.MAX_ERROR_RATE,
                20,
                6,
                0.1,
                6);
        PolyAFinder.PolyARun ret = polyAFinder.getPolyAStart(read, adapter);
        Assert.assertEquals(ret.startPos, expectedResult, testName);
    }

    @DataProvider(name="testWithAdapterTestCases")
    public Object[][] testWithAdapterTestCases() {
        return new Object[][]{
                {"trailingPolyAWithoutAdapter", "AAAACATGTTGATTATTATTTTTATTAAATTAATTACAATAAAAAATAAA",
                        "CGCCTAACTGGCAATTCAATACGTACTCTGCGTTGCTACCACTG", 40},
                {"completeTrimWithPolyA", "AAACCTTACACCCTTCTAATTCCACGTACTCTGCGTTGATACCACTGCTTCC",
                        "CCTTACACCCTTCTAATTCCACGTACTCTGCGTTGCTACCACTG", 0},
                {"completeTrimWithoutPolyA", "CCTTACACCCTTCTAATTCCACGTACTCTGCGTTGATACCACTGCTTCC",
                        "CCTTACACCCTTCTAATTCCACGTACTCTGCGTTGCTACCACTG", 0},
                {"trailingPolyAWithAdapter", "CGTCTATGTAGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCAGG",
                        "CCCCCCAGGACCGACCCGAGACGTACTCTGCGTTGCTACCACTG", 9},
                {"noTrim", "GGTTTGATAAAGAACTTAGAAAAAAAAAAAAAAAAAAAATTTTCTAACGT",
                        "CCCCCGTCGTCGCATAGGTAACGTACTCTGCGTTGCTACCACTG", -1}
        };
    }
}
