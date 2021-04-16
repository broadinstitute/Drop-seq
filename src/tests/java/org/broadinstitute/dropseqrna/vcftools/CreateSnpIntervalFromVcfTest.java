/*
 * MIT License
 *
 * Copyright 2019 Broad Institute
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
package org.broadinstitute.dropseqrna.vcftools;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

public class CreateSnpIntervalFromVcfTest {

    private static final File TEST_FILE = new File("testdata/org/broadinstitute/vcftools/test.vcf");
    private static final String EXPECTED_CONTIG = "1";
    private static final int[] EXPECTED_START_POS = {
            243901, 76227022, 150199123, 150199124
    };
    
    private static final int[] EXPECTED_START_POS_HET_ONLY = {
            76227022
    };
    

    @Test
    public void testBasic() throws IOException {
        final CreateSnpIntervalFromVcf clp = new CreateSnpIntervalFromVcf();
        clp.INPUT = TEST_FILE;
        clp.OUTPUT = File.createTempFile("CreateSnpIntervalFromVcfTest.", ".intervals");
        clp.OUTPUT.deleteOnExit();
        Assert.assertEquals(clp.doWork(), 0);
        final IntervalList intervals = IntervalList.fromFile(clp.OUTPUT);
        // 6 variants in input, but one is an indel and one is filtered
        Assert.assertEquals(intervals.size(), EXPECTED_START_POS.length);
        final Iterator<Interval> it = intervals.iterator();
        for (final int startPos : EXPECTED_START_POS) {
            Assert.assertTrue(it.hasNext());
            final Interval interval = it.next();
            Assert.assertEquals(interval.getContig(), EXPECTED_CONTIG);
            Assert.assertEquals(interval.getStart(), startPos);
            Assert.assertEquals(interval.length(), 1);
        }
        Assert.assertFalse(it.hasNext());
    }
    
    @Test
    public void testHetOnly() throws IOException {
        final CreateSnpIntervalFromVcf clp = new CreateSnpIntervalFromVcf();
        clp.INPUT = TEST_FILE;
        clp.HET_SNPS_ONLY=true;
        clp.OUTPUT = File.createTempFile("CreateSnpIntervalFromVcfTest.", ".intervals");
        clp.OUTPUT.deleteOnExit();
        Assert.assertEquals(clp.doWork(), 0);
        final IntervalList intervals = IntervalList.fromFile(clp.OUTPUT);
        // 6 variants in input, but one is an indel and one is filtered
        Assert.assertEquals(intervals.size(), EXPECTED_START_POS_HET_ONLY.length);
        final Iterator<Interval> it = intervals.iterator();
        for (final int startPos : EXPECTED_START_POS_HET_ONLY) {
            Assert.assertTrue(it.hasNext());
            final Interval interval = it.next();
            Assert.assertEquals(interval.getContig(), EXPECTED_CONTIG);
            Assert.assertEquals(interval.getStart(), startPos);
            Assert.assertEquals(interval.length(), 1);
        }
        Assert.assertFalse(it.hasNext());
    }
    
    
    
    
}
