/*
 * MIT License
 *
 * Copyright 2026 Broad Institute
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

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Random;

public class HammingDistanceTest {

    private static String randomBases(final Random r, final int length) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < length; i++) {
            sb.append("ACGT".charAt(r.nextInt(4)));
        }
        return sb.toString();
    }

    /** withinHammingDistance must give the same answer as comparing getHammingDistance to the threshold. */
    @Test
    public void testWithinHammingDistanceMatchesHammingDistance() {
        final Random r = new Random(1);
        for (int i = 0; i < 20_000; i++) {
            final String a = randomBases(r, r.nextInt(8));
            // a mutated copy is often close to a, so that small distances are well covered
            final String b = r.nextBoolean() ? randomBases(r, r.nextInt(8)) : mutate(r, a);
            final int distance = HammingDistance.getHammingDistance(a, b);
            for (int max = -1; max <= 8; max++) {
                Assert.assertEquals(HammingDistance.withinHammingDistance(a, b, max), distance <= max,
                        a + " " + b + " max " + max + " distance " + distance);
            }
        }
    }

    private static String mutate(final Random r, final String s) {
        final StringBuilder sb = new StringBuilder(s);
        for (int i = 0; i < sb.length(); i++) {
            if (r.nextInt(4) == 0) {
                sb.setCharAt(i, "ACGT".charAt(r.nextInt(4)));
            }
        }
        return sb.toString();
    }

    @Test
    public void testWithinHammingDistanceExamples() {
        Assert.assertTrue(HammingDistance.withinHammingDistance("ACGT", "ACGT", 0));
        Assert.assertFalse(HammingDistance.withinHammingDistance("ACGT", "ACGA", 0));
        Assert.assertTrue(HammingDistance.withinHammingDistance("ACGT", "ACGA", 1));
        Assert.assertFalse(HammingDistance.withinHammingDistance("ACGT", "AGGA", 1));
        // each position of the longer string past the shorter one is a difference
        Assert.assertTrue(HammingDistance.withinHammingDistance("ACG", "ACGT", 1));
        Assert.assertFalse(HammingDistance.withinHammingDistance("AC", "ACGT", 1));
        Assert.assertTrue(HammingDistance.withinHammingDistance("", "", 0));
    }
}
