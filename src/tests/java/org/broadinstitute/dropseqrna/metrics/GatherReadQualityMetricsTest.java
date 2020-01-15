package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.io.IOException;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.annotations.Test;

import junit.framework.Assert;

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
		//g.OUTPUT.deleteOnExit();
		g.TAG="XC";

		int r = g.doWork();
		Assert.assertTrue(r==0);
		boolean t1 = TestUtils.testFilesSame(EXPECTED_RESULT, outFile);
		Assert.assertTrue(t1);

	}



}
