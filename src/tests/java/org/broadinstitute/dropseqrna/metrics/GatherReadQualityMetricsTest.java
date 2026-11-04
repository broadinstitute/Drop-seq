package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.annotations.Test;

import org.testng.Assert;
import picard.nio.PicardHtsPath;

public class GatherReadQualityMetricsTest {

	private static final Path IN_FILE = Paths.get("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.bam");
	private static final File EXPECTED_RESULT = new File ("testdata/org/broadinstitute/dropseq/metrics/5cell3gene.read_quality_metrics.txt");

	@Test
	public void testDoWork() throws IOException {
		File outFile = File.createTempFile("GatherReadQualityMetricsTest.", ".read_quality_metrics.txt");
		GatherReadQualityMetrics g = new GatherReadQualityMetrics();
		g.INPUT = Collections.singletonList(PicardHtsPath.fromPath(IN_FILE));
		g.MINIMUM_MAPPING_QUALITY=10;
		g.OUTPUT=outFile;
		g.OUTPUT.deleteOnExit();
		g.TAG="XC";

		int r = g.doWork();
		Assert.assertTrue(r==0);
		boolean t1 = TestUtils.testFilesSame(EXPECTED_RESULT, outFile);
		Assert.assertTrue(t1);

	}

	private static File run(final List<PicardHtsPath> inputs, final int numThreads) throws IOException {
		File outFile = File.createTempFile("GatherReadQualityMetricsTest.", ".read_quality_metrics.txt");
		outFile.deleteOnExit();
		GatherReadQualityMetrics g = new GatherReadQualityMetrics();
		g.INPUT = inputs;
		g.MINIMUM_MAPPING_QUALITY = 10;
		g.OUTPUT = outFile;
		g.TAG = "XC";
		g.NUM_THREADS = numThreads;
		Assert.assertEquals(g.doWork(), 0);
		return outFile;
	}

	@Test
	public void testNumThreadsDoesNotChangeOutput() throws IOException {
		for (final int numThreads : new int[]{1, 2, 4}) {
			File outFile = run(Collections.singletonList(PicardHtsPath.fromPath(IN_FILE)), numThreads);
			Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_RESULT, outFile), "NUM_THREADS=" + numThreads);
		}
	}

	@Test
	public void testNumThreadsDoesNotChangeOutputWithMultipleInputs() throws IOException {
		final List<PicardHtsPath> inputs = Arrays.asList(
				PicardHtsPath.fromPath(IN_FILE), PicardHtsPath.fromPath(IN_FILE), PicardHtsPath.fromPath(IN_FILE));
		final File serial = run(inputs, 1);
		final File threaded = run(inputs, 4);
		Assert.assertTrue(TestUtils.testFilesSame(serial, threaded));
		// three copies of the input must not give the same output as one, or the comparison above proves nothing
		Assert.assertFalse(TestUtils.testFilesSame(EXPECTED_RESULT, threaded));
	}

	@Test
	public void testNumThreadsMustBePositive() {
		GatherReadQualityMetrics g = new GatherReadQualityMetrics();
		g.NUM_THREADS = 0;
		Assert.assertTrue(g.customCommandLineValidation() != null);
	}



}
