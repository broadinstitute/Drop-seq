package org.broadinstitute.dropseqrna.utils.alignmentcomparison;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class CompareDropSeqAlignmentsTest {

	// modify reads by hand, then test for changes.
	// NS500217:67:H14GMBGXX:1:22207:3769:12483 HUMAN_3 -> HUMAN_4
	// NS500217:67:H14GMBGXX:3:21507:14155:11331 HUMAN_3 -> HUMAN_22 also NKTR->NKTR2
	// NS500217:67:H14GMBGXX:2:23305:7629:7836 16 SNRPA1 -> SNRPA2
	// NS500217:67:H14GMBGXX:2:22312:14606:19968 HUMAN_10 -> HUMAN_GL000224.1
	// NS500217:67:H14GMBGXX:3:11507:4353:14582 HUMAN_15 -> HUMAN_16, SNRPA1 -> SNRPA2

	private File OLD = new File ("testdata/org/broadinstitute/dropseq/utils/alignmentcomparison/old_alignment.bam");
	private File NEW = new File ("testdata/org/broadinstitute/dropseq/utils/alignmentcomparison/new_alignment.bam");

	private File CONTIG_REPORT = new File ("testdata/org/broadinstitute/dropseq/utils/alignmentcomparison/contig_report.txt");
	private File GENE_REPORT = new File ("testdata/org/broadinstitute/dropseq/utils/alignmentcomparison/gene_report.txt");

	@Test
	public void testDoWork() {
		CompareDropSeqAlignments c = new CompareDropSeqAlignments();
		File outGeneReport=null;
		File outContigReport=null;
		try {
			outGeneReport = File.createTempFile("CompareDropSeqAlignmentsTest.", ".gene_report.txt");
			outContigReport = File.createTempFile("CompareDropSeqAlignmentsTest.", ".contig_report.txt");
		} catch (IOException e) {
			e.printStackTrace();
		}


		c.INPUT_1=OLD;
		c.INPUT_2=NEW;
		c.GENE_REPORT=outGeneReport;
		c.CONTIG_REPORT=outContigReport;
		int ret = c.doWork();
		Assert.assertTrue(ret==0);

		boolean t1 = TestUtils.testFilesSame(this.CONTIG_REPORT, outContigReport);
		boolean t2 = TestUtils.testFilesSame(this.GENE_REPORT, outGeneReport);
		Assert.assertTrue(t1);
		Assert.assertTrue(t2);

	}

	@Test
	public void testGeneResultAddMapping() {
		GeneResult r = new GeneResult("GeneA", "1", CompareDropSeqAlignments.noGeneTag);
		r.addMapping(Arrays.asList("GeneA"), Arrays.asList("1"), 1);
		Assert.assertSame(r.getCountSameMapping(), 1);

		r = new GeneResult("GeneA", "1", CompareDropSeqAlignments.noGeneTag);
		r.addMapping(Arrays.asList("GeneA"), Arrays.asList("1"), 2);
		Assert.assertSame(r.getCountSameGeneMapsNonUniqueCount(), 1);

		r = new GeneResult("GeneA", "1", CompareDropSeqAlignments.noGeneTag);
		r.addMapping(Arrays.asList("GeneA", CompareDropSeqAlignments.noGeneTag), Arrays.asList("1"), 2);
		Assert.assertSame(r.getCountSameGeneMapsNonUniqueCount(), 1);

		r = new GeneResult("GeneA", "1", CompareDropSeqAlignments.noGeneTag);
		r.addMapping(Arrays.asList(CompareDropSeqAlignments.noGeneTag), Arrays.asList("1"), 2);
		Assert.assertSame(r.getCountIntronicOrIntergenic(), 1);

		r = new GeneResult("GeneA", "1", CompareDropSeqAlignments.noGeneTag);
		r.addMapping(Arrays.asList(CompareDropSeqAlignments.noGeneTag), Arrays.asList("1"), 1);
		Assert.assertSame(r.getCountIntronicOrIntergenic(), 1);

		r = new GeneResult("GeneA", "1", CompareDropSeqAlignments.noGeneTag);
		r.addMapping(Arrays.asList("GeneB"), Arrays.asList("1"), 1);
		Assert.assertSame(r.getCountDifferentUniqueGene(), 1);

		r = new GeneResult("GeneA", "1", CompareDropSeqAlignments.noGeneTag);
		r.addMapping(Arrays.asList("GeneB"), Arrays.asList("1"), 2);
		Assert.assertSame(r.getCountDifferentGeneNonUniqueCount(), 1);

		r = new GeneResult("GeneA", "1", CompareDropSeqAlignments.noGeneTag);
		r.addMapping(Arrays.asList("GeneA", "GeneB"), Arrays.asList("1"), 2);
		Assert.assertSame(r.getCountMultiGeneMappingCount(), 1);

		Assert.assertNotNull(r.toString());

	}

	@Test
	public void testContigResult () {
		ContigResult r = new ContigResult("1", Arrays.asList("2"), false);
		ContigResult r2 = new ContigResult("1", Arrays.asList("2"), true);

		Assert.assertTrue(r.compareTo(r)==0);
		Assert.assertTrue(r.compareTo(r2)>0);
		Assert.assertTrue(r2.compareTo(r)<0);

		Assert.assertFalse(r.equals(r2));
		Assert.assertTrue(r.equals(r));

		Assert.assertNotNull(r.toString());

	}




}
