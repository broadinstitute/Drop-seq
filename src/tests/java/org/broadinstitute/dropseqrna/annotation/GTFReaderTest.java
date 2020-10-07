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
import htsjdk.samtools.util.OverlapDetector;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.annotation.AnnotationException;
import picard.annotation.Gene;
import picard.annotation.Gene.Transcript;

import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;



public class GTFReaderTest {
	File TEST_DATA_DIR = new File("testdata/org/broadinstitute/transcriptome/annotation");

	private final File GTF_FILE1 = new File(TEST_DATA_DIR, "human_ISG15.gtf.gz");
	private final File GTF_FILE2 = new File(TEST_DATA_DIR, "human_ISG15_FAM41C.gtf.gz");
	private final File GTF_FILE3 = new File(TEST_DATA_DIR, "human_SNORD18.gtf.gz");
	private final File SD = new File(TEST_DATA_DIR, "human_g1k_v37_decoy_50.dict");
	File PSEUDOGENE_GTF = new File(TEST_DATA_DIR, "pseudogene.gtf");



	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public void test1() {
		GTFReader r = new GTFReader(GTF_FILE1, SD);
		OverlapDetector<Gene> od = (OverlapDetector) r.load();
		Assert.assertNotNull(od);
		Collection<Gene> genes = od.getAll();
		Assert.assertEquals(genes.size(), 1);
		Gene g = genes.iterator().next();
		Assert.assertEquals(g.getStart(), 948803);
		Assert.assertEquals(g.getEnd(), 949920);
		Assert.assertTrue(g.isPositiveStrand());
		Gene.Transcript t = g.iterator().next();
		Assert.assertEquals(t.transcriptionStart, 948803);
		Assert.assertEquals(t.transcriptionEnd, 949920);
	}



	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	// I like the negative strand genes da best, so I put one in this set.
	public void test2() {
		GTFReader r = new GTFReader(GTF_FILE2, SD);

		@SuppressWarnings({ "rawtypes", "unchecked" })
		OverlapDetector<Gene> od = (OverlapDetector) r.load();
		Assert.assertNotNull(od);
		Collection<Gene> genes = od.getAll();
		Map<String, Gene> geneMap = new HashMap<>(genes.size());
		for (Gene g: genes)
			geneMap.put(g.getName(), g);

		Gene g = geneMap.get("FAM41C");

		Assert.assertEquals(genes.size(), 2);

		Assert.assertEquals(g.getStart(), 803451);
		Assert.assertEquals(g.getEnd(), 812283);
		Assert.assertTrue(g.isNegativeStrand());
        for (Transcript t : g)
			if (t.name.equals("FAM41C-001")) {
                Assert.assertEquals(t.name, "FAM41C-001");
                Assert.assertEquals(t.transcriptionStart, 803451);
                Assert.assertEquals(t.transcriptionEnd, 812283);
            }
	}

	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	// I like the negative strand genes da best, so I put one in this set.
	public void testMixedChromosomeGene() {
		GTFReader r = new GTFReader(GTF_FILE3, SD);

		@SuppressWarnings({ "rawtypes", "unchecked" })
		OverlapDetector<Gene> od = (OverlapDetector) r.load();
		Assert.assertNotNull(od);
		Collection<Gene> genes = od.getAll();
		Map<String, Gene> geneMap = new HashMap<>(genes.size());
		for (Gene g: genes)
			geneMap.put(g.getName(), g);

		Gene g = geneMap.get("SNORD18");

		/*
		Assert.assertEquals(genes.size(), 2);

		Assert.assertEquals(g.getStart(), 803451);
		Assert.assertEquals(g.getEnd(), 812283);
		Assert.assertTrue(g.isNegativeStrand());
		Iterator<Transcript> iter = g.iterator();
		while (iter.hasNext()) {
			Transcript t = iter.next();
			if (t.name.equals("FAM41C-001")) {
				Assert.assertEquals(t.name, "FAM41C-001");
				Assert.assertEquals(t.transcriptionStart, 803451);
				Assert.assertEquals(t.transcriptionEnd, 812283);
			}
		}
		*/
	}



	@Test(enabled=false, groups={"dropseq", "transcriptome"})
	// this tests loading in a full GTF, which I'm not going to store in the code base itself, so I'll leave the test disabled.
	// this is more to see what errors are thrown in a full parse, rather than to debug individual elements.
	public void fullLoad () {
		File fullD = new File ("/humgen/cnp04/sandbox/data/Evan/common/individual_references/Human/Homo_sapiens.GRCh37.74.gtf");
		GTFReader r = new GTFReader(fullD, SD);
		@SuppressWarnings({ "rawtypes", "unchecked" })
		OverlapDetector<Gene> od = (OverlapDetector) r.load();
		Assert.assertNotNull(od);
	}

	@Test
	public void testGeneWithNoExonTranscriptRejectedLenient() {
		final OverlapDetector<GeneFromGTF> od = load(PSEUDOGENE_GTF, SD, ValidationStringency.LENIENT);
		Assert.assertTrue(od.getAll().isEmpty());
	}

	@Test(expectedExceptions = AnnotationException.class)
	public void testGeneWithNoExonTranscriptRejectedStrict() {
		final OverlapDetector<GeneFromGTF> od = load(PSEUDOGENE_GTF, SD, ValidationStringency.STRICT);
		Assert.assertTrue(od.getAll().isEmpty());
	}

	private OverlapDetector<GeneFromGTF> load(final File gtf, final File sd, final ValidationStringency stringency) {
		final GTFReader reader = new GTFReader(gtf, sd);
		reader.setValidationStringency(stringency);
		return reader.load();
	}

}
