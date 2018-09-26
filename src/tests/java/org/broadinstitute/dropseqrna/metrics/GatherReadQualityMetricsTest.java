package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.io.IOException;

import org.testng.annotations.Test;

import junit.framework.Assert;
import picard.util.TabbedInputParser;

public class GatherReadQualityMetricsTest {

	private static final File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.bam");
	private static final File EXPECTED_RESULT = new File ("testdata/org/broadinstitute/dropseq/metrics/5cell3gene.read_quality_metrics.txt");

	@Test
	public void testDoWork() throws IOException {
		File outFile = File.createTempFile("GatherReadQualityMetricsTest.", ".read_quality_metrics.txt");
		GatherReadQualityMetrics g = new GatherReadQualityMetrics();
		g.INPUT=IN_FILE;
		g.MINIMUM_MAPPING_QUALITY=10;
		g.OUTPUT=outFile;
		g.TAG="XC";

		int r = g.doWork();
		Assert.assertTrue(r==0);
		boolean t1 = testFilesSame(EXPECTED_RESULT, outFile);
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
