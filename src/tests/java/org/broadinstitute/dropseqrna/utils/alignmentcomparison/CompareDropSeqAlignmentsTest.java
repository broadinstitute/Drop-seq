package org.broadinstitute.dropseqrna.utils.alignmentcomparison;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import picard.util.TabbedInputParser;

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

		boolean t1 = testFilesSame(this.CONTIG_REPORT, outContigReport);
		boolean t2 = testFilesSame(this.GENE_REPORT, outGeneReport);
		Assert.assertTrue(t1);
		Assert.assertTrue(t2);

	}

	// Like FileUtils.contentEquals(file1, file2), but ignores the header lines which may be different due to absolute file paths.
	private boolean testFilesSame (final File expected, final File actual) {
		TabbedInputParser e = new TabbedInputParser(true, expected);
		TabbedInputParser a = new TabbedInputParser(true, expected);

		while (e.hasNext() && a.hasNext()) {
			e.next();
			a.next();
			String le = e.getCurrentLine();
			String la = a.getCurrentLine();
			if (!le.equals(la)) return false;
		}

		return true;
	}

}
