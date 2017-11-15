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

import java.io.File;
import java.util.Collection;

import org.broadinstitute.dropseqrna.annotation.ConvertToRefFlat;
import org.broadinstitute.dropseqrna.annotation.GeneAnnotationReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import picard.annotation.Gene;


public class ConvertToRefFlatTest {

	File GTF_FILE1 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_ISG15.gtf.gz");
	// cat Homo_sapiens_assembly19.refFlat |grep ISG15 > ~/1chip/trunk/transcriptome_java/testdata/org/broadinstitute/transcriptome/annotation/human_ISG15_refFlat
	File REF_FLAT1 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_ISG15.refFlat.gz");
	File SD = new File("testdata/org/broadinstitute/transcriptome/annotation/human_g1k_v37_decoy_50.dict");
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testConvert() {
		ConvertToRefFlat crf = new ConvertToRefFlat();
		
		// GTFReader r = new GTFReader(GTF_FILE1, SD);
		OverlapDetector<Gene> od = GeneAnnotationReader.loadAnnotationsFile(GTF_FILE1, SD);
		Assert.assertNotNull(od);
		
		Collection<Gene> genes = od.getAll();
		Gene g1 = genes.iterator().next();
		Gene.Transcript t = g1.iterator().next();
		
		OverlapDetector<Gene> od2 = GeneAnnotationReader.loadAnnotationsFile(REF_FLAT1, SD);
		Gene.Transcript t2  = od2.getAll().iterator().next().iterator().next();
		
		//assertTwoTranscriptsSame(t, t2);
		
		RefFlatRecord  result = crf.convertTranscript(t);
		Assert.assertNotNull(result);
		
		RefFlatRecord expected = new RefFlatRecord("ISG15", "ISG15-001", "1", false, 948803, 949920, 948954, 949855);
		expected.addExonStart(948802+1);
		expected.addExonStart(949363+1);
		expected.addExonEnd(948956);
		expected.addExonEnd(949920);
		
		Assert.assertEquals(result, expected);
		
		RefFlatRecord  result2 = crf.convertTranscript(t2);
		Assert.assertNotNull(result2);
		// I add data in 1 space, the file is in 0 space.  Thus the very obvious +1's.
		RefFlatRecord expected2 = new RefFlatRecord("ISG15", "NM_005101", "1", false, 948846+1, 949919, 948953+1, 949858);
		expected2.addExonStart(948846+1);
		expected2.addExonStart(949363+1);
		expected2.addExonEnd(948956);
		expected2.addExonEnd(949919);
		
		Assert.assertEquals(result2, expected2);
		
		
	}
	
	/*
	private void assertTwoTranscriptsSame (Transcript t, Transcript t2) {
		// the transcript names don't match up for these guys, but the other crap does.
		
		//Assert.assertEquals(t.getGene().getName(), t2.getGene().getName());
		Assert.assertEquals(t.getGene().getSequence(), t2.getGene().getSequence());
		Assert.assertEquals(t.getGene().getStart(), t2.getGene().getStart());
		Assert.assertEquals(t.getGene().getEnd(), t2.getGene().getEnd());
		Assert.assertEquals(t.codingStart, t2.codingStart);
		Assert.assertEquals(t.codingEnd, t2.codingEnd);
		Assert.assertEquals(t.transcriptionStart, t2.transcriptionStart);
		Assert.assertEquals(t.transcriptionEnd, t2.transcriptionEnd);

		Exon [] e1 = t.exons;
		Exon [] e2 = t2.exons;

		Assert.assertEquals(e1.length, e2.length);
		for (int i=0; i<e1.length; i++) {
			Assert.assertEquals(e1[i].start, e2[i].start);
			Assert.assertEquals(e1[i].end, e2[i].end);
		}
	}
	*/
}
