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
package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class EnhanceGTFRecordsTest {

	File GTF_FILE1 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_ISG15.gtf.gz");
	File GTF_FILE2 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_ISG15_FAM41C.gtf.gz");
	File SD = new File("testdata/org/broadinstitute/transcriptome/annotation/human_g1k_v37_decoy_50.dict");
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testIntron() {
		List<Interval> exons = new ArrayList<Interval>();
		Interval e1= new Interval(null, 948803, 948956);
		Interval e2= new Interval(null, 949364, 949920);
		exons.add(e1);
		exons.add(e2);
		
		EnhanceGTFRecords e = new EnhanceGTFRecords();
		List<Interval> introns = e.getIntronIntervals(exons);
		
		Assert.assertEquals(exons.size()-1, introns.size());
		Interval intron = introns.get(0);
		Assert.assertEquals(intron.getStart(), 948957);
		Assert.assertEquals(intron.getEnd(), 949363);
		
	}
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void test1Enhanced() {
		EnhanceGTFRecords e = new EnhanceGTFRecords();
		GTFParser parser = new GTFParser(GTF_FILE1, ValidationStringency.STRICT);
        List<GTFRecord> records;
        try {
            records = e.enhanceGTFRecords(parser);
        } finally {
            CloserUtil.close(parser);
        }
        Assert.assertNotNull(records);
		
	}
	
}
