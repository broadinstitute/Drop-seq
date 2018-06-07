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

import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.dropseqrna.annotation.GeneAnnotationReader;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.annotation.Gene;

import java.io.File;


public class GeneAnnotationReaderTest {

	File GTF_FILE1 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_ISG15.gtf.gz");
	File refFlatCompressed = new File("testdata/org/broadinstitute/transcriptome/annotation/Homo_sapiens.GRCh37.74.refFlat.gz");
	File refFlatUncompressed = new File("testdata/org/broadinstitute/transcriptome/annotation/Homo_sapiens.GRCh37.74.refFlat");
	File SD = new File("testdata/org/broadinstitute/transcriptome/annotation/human_g1k_v37_decoy_50.dict");
	
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testGTF() {
		OverlapDetector<Gene> od = GeneAnnotationReader.loadAnnotationsFile(GTF_FILE1, SD);
        Assert.assertNotNull(od);
	}
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testRefFlatCompressed() {
		OverlapDetector<Gene> od = GeneAnnotationReader.loadAnnotationsFile(refFlatCompressed, SD);
        Assert.assertNotNull(od);
	}
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testRefFlatUncompressed() {
		OverlapDetector<Gene> od = GeneAnnotationReader.loadAnnotationsFile(refFlatUncompressed, SD);
        Assert.assertNotNull(od);
	}
}
