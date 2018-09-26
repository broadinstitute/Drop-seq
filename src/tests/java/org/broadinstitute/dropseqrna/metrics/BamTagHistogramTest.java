package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.io.IOException;

import org.testng.annotations.Test;

import junit.framework.Assert;
import picard.util.TabbedInputParser;

public class BamTagHistogramTest {

	private static final File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.bam");
	private static final File EXPECTED_RESULT = new File ("testdata/org/broadinstitute/dropseq/metrics/5cell3gene.counts_per_XC.txt");
	private static final File EXPECTED_RESULT_INT = new File ("testdata/org/broadinstitute/dropseq/metrics/5cell3gene.counts_per_NM.txt");

	@Test
	public void doWorkStringTag() throws IOException {
		File outFile = File.createTempFile("BamTagHistogramTest.", ".counts_XC.txt");
		outFile.deleteOnExit();
		BamTagHistogram bth = new BamTagHistogram();
		bth.INPUT=IN_FILE;
		bth.OUTPUT=outFile;
		bth.MINIMUM_MAPPING_QUALITY=10;
		bth.TAG="XC";

		int r = bth.doWork();
		Assert.assertTrue(r==0);
		boolean t1 = testFilesSame(EXPECTED_RESULT, outFile);
		Assert.assertTrue(t1);
	}

	@Test
	public void doWorkIntegerTag() throws IOException {
		File outFile=null;
		outFile = File.createTempFile("BamTagHistogramTest.", ".counts_NM.txt");
		outFile.deleteOnExit();
		BamTagHistogram bth = new BamTagHistogram();
		bth.INPUT=IN_FILE;
		bth.OUTPUT=outFile;
		bth.MINIMUM_MAPPING_QUALITY=10;
		bth.TAG="NM";

		int r = bth.doWork();
		Assert.assertTrue(r==0);
		boolean t1 = testFilesSame(EXPECTED_RESULT_INT, outFile);
		Assert.assertTrue(t1);
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
