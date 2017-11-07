/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
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
