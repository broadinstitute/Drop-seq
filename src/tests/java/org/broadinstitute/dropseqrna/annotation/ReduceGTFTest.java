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

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.CollectionUtil;


public class ReduceGTFTest {

	File GTF_FILE1 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_ISG15.gtf.gz");
	File GTF_FILE2 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_ISG15_FAM41C.gtf.gz");
	File GTF_FILE3 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_AL592188.5.gtf.gz");
	// contains just the APITD1 gene
	File GTF_FILE4 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_APITD1.gtf.gz");
	// contains both the APITD1 gene and the APITD1-CORT genes, that occupy the same bounds, but have different sets of exons.
	File GTF_FILE5 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_APITD1_both.gtf.gz");
	File GTF_FILE5_REDUCED= new File("testdata/org/broadinstitute/transcriptome/annotation/human_APITD1_both.gtf.reduced.gz");

	File SD = new File("testdata/org/broadinstitute/transcriptome/annotation/human_g1k_v37_decoy_50.dict");



    private Iterator<GTFRecord> parseGtf(final File file) {
        ReduceGTF r = new ReduceGTF();
        r.SEQUENCE_DICTIONARY = SD;
        r.GTF = file;
        return r.parseGTF();



    }

    @Test
    public void testDoWork () {
    	ReduceGTF r = new ReduceGTF();
    	r.GTF=GTF_FILE5;
    	r.SEQUENCE_DICTIONARY=SD;
    	File o = getTempReportFile("ReduceGTF", ".reduced.gtf.gz");
    	// File o2 = getTempReportFile("ReduceGTF", ".reduced.gtf");
    	o.deleteOnExit();
    	// o2.deleteOnExit();
    	//r.OUTPUT=GTF_FILE5_REDUCED;
    	r.OUTPUT=o;
    	r.ENHANCE_GTF=true;
    	int ret = r.doWork();
    	Assert.assertTrue(ret==0);
//    	r.OUTPUT=o2;
//    	int ret2 = r.doWork();
//    	Assert.assertTrue(ret2==0);


    	try {
			boolean t1 = FileUtils.contentEquals(o, GTF_FILE5_REDUCED);
			Assert.assertTrue(t1);
		} catch (IOException e) {
			e.printStackTrace();
		}


    }

    private File getTempReportFile (final String prefix, final String suffix) {
		File tempFile=null;

		try {
			tempFile = File.createTempFile(prefix, suffix);
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		return tempFile;
	}


	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void test1() {
		Assert.assertNotNull(parseGtf(GTF_FILE1));

	}

	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void test2() {
		final Iterator<GTFRecord> gtfIterator = parseGtf(GTF_FILE3);
		Assert.assertNotNull(gtfIterator);
		EnhanceGTFRecords e = new EnhanceGTFRecords();
		List<GTFRecord> records = e.enhanceGTFRecords(gtfIterator);
		Assert.assertNotNull(records);

		for (GTFRecord a: records) {
			Assert.assertEquals(a.getStart(), 152252);
			Assert.assertEquals(a.getEnd(), 152312);
			Assert.assertEquals(a.getStrandAsString(), "+");

			if (a.getFeatureType().equals("gene")) {
				Assert.assertEquals(a.getTranscriptID(), null);
				Assert.assertEquals(a.getTranscriptName(), null);
				Assert.assertEquals(a.getTranscriptType(), "miRNA");
			} else {
				Assert.assertEquals(a.getTranscriptID(), "ENST00000577630");
				Assert.assertEquals(a.getTranscriptName(), "AL592188.5-201");
			}

		}

	}

	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testAPITD1() {
		Iterator<GTFRecord> gtfIterator = parseGtf(GTF_FILE4);
		Assert.assertNotNull(gtfIterator);
        Collection<GTFRecord> records = CollectionUtil.makeCollection(gtfIterator);
		// gunzip -c human_APITD1.gtf.gz | grep -v CDS |grep -v start_codon |grep -v stop_codon |wc -l
		Assert.assertEquals(records.size(),26);

        final GeneFromGTFBuilder geneBuilder = new GeneFromGTFBuilder(records.iterator());
        Collection<GeneFromGTF> genes = CollectionUtil.makeCollection(geneBuilder);
		Assert.assertEquals(genes.size(),1);

        final EnhanceGTFRecords enhancer = new EnhanceGTFRecords();
        for (final GeneFromGTF gene : genes)
			Assert.assertNotNull(enhancer.enhanceGene(gene));
	}

	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testAPITD1Complex() {
        Iterator<GTFRecord> gtfIterator = parseGtf(GTF_FILE5);
		Assert.assertNotNull(gtfIterator);
        Collection<GTFRecord> records = CollectionUtil.makeCollection(gtfIterator);
		// gunzip -c human_APITD1.gtf.gz | grep -v CDS |grep -v start_codon |grep -v stop_codon |wc -l
		Assert.assertEquals(records.size(),42);

        final GeneFromGTFBuilder geneBuilder = new GeneFromGTFBuilder(records.iterator());
        Collection<GeneFromGTF> genes = CollectionUtil.makeCollection(geneBuilder);
        Assert.assertEquals(genes.size(),2);

        final EnhanceGTFRecords enhancer = new EnhanceGTFRecords();
        for (final GeneFromGTF gene : genes)
			Assert.assertNotNull(enhancer.enhanceGene(gene));
	}
}
