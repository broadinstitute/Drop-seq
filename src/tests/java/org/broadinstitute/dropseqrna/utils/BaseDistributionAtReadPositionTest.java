package org.broadinstitute.dropseqrna.utils;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class BaseDistributionAtReadPositionTest {

	private static final File TAGGED_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.bam");
	private static final File UNMAPPED_BAM_FILE = new File("testdata/org/broadinstitute/dropseq/utils/unmapped_paired_reads.bam");
	private static final File OUTPUT_EXPECTED_RESULT=new File("testdata/org/broadinstitute/dropseq/utils/BaseDistributionAtReadPosition.expected_output.txt");

	@Test
	public void gatherBaseQualitiesFile() {
		BaseDistributionAtReadPosition test= new BaseDistributionAtReadPosition();
		BaseDistributionMetricCollection resultR1 = test.gatherBaseQualities (UNMAPPED_BAM_FILE, 1);
		// {1={A=1, C=1, T=4, G=4, N=0}, 2={A=4, C=2, T=1, G=3, N=0}, 3={A=1, C=4, T=3, G=2, N=0}, 4={A=4, C=3, T=1, G=2, N=0}, 5={A=1, C=6, T=2, G=1, N=0}, 6={A=2, C=2, T=1, G=5, N=0}, 7={A=4, C=1, T=3, G=2, N=0}, 8={A=2, C=3, T=1, G=4, N=0}, 9={A=3, C=1, T=4, G=2, N=0}, 10={A=3, C=2, T=2, G=3, N=0}}
		Assert.assertEquals(new BaseDistributionMetric(1,1,4,4,0), resultR1.getDistributionAtPosition(1));
		Assert.assertEquals(new BaseDistributionMetric(4,2,3,1,0), resultR1.getDistributionAtPosition(2));
		Assert.assertEquals(new BaseDistributionMetric(4,3,2,1,0), resultR1.getDistributionAtPosition(4));
		Assert.assertEquals(new BaseDistributionMetric(2,3,4,1,0), resultR1.getDistributionAtPosition(8));

		BaseDistributionMetricCollection resultR2 = test.gatherBaseQualities (UNMAPPED_BAM_FILE, 2);
		// {1={A=1, C=3, T=2, G=4, N=0}, 2={A=3, C=3, T=1, G=3, N=0}, 3={A=1, C=3, T=0, G=6, N=0}, 4={A=2, C=4, T=2, G=2, N=0}, 5={A=4, C=3, T=2, G=1, N=0}, 6={A=2, C=2, T=3, G=3, N=0}, 7={A=1, C=4, T=0, G=5, N=0}, 8={A=1, C=2, T=6, G=1, N=0}, 9={A=2, C=2, T=3, G=3, N=0}, 10={A=2, C=2, T=4, G=2, N=0}}
		Assert.assertEquals(new BaseDistributionMetric(1,3,4,2,0), resultR2.getDistributionAtPosition(1));
		Assert.assertEquals(new BaseDistributionMetric(4,3,1,2,0), resultR2.getDistributionAtPosition(5));
		Assert.assertEquals(new BaseDistributionMetric(2,2,3,3,0), resultR2.getDistributionAtPosition(9));
	}

	@Test
	public void gatherBaseQualitiesFileString() {
		BaseDistributionAtReadPosition test= new BaseDistributionAtReadPosition();
		BaseDistributionMetricCollection result = test.gatherBaseQualities (TAGGED_FILE, "XC");
		// {1={A=2619, C=0, T=3747, G=4, N=0}, 2={A=1097, C=3, T=3209, G=2061, N=0}, 3={A=24, C=2557, T=0, G=3789, N=0}, 4={A=2592, C=2709, T=3, G=1066, N=0}, 5={A=2148, C=1746, T=23, G=2453, N=0}, 6={A=2072, C=5, T=2806, G=1487, N=0}, 7={A=2051, C=19, T=2802, G=1498, N=0}, 8={A=5405, C=15, T=18, G=932, N=0}, 9={A=2012, C=3216, T=1076, G=66, N=0}, 10={A=1504, C=27, T=1078, G=3761, N=0}, 11={A=961, C=1764, T=39, G=3606, N=0}, 12={A=2524, C=1097, T=968, G=1781, N=0}}
		Assert.assertEquals(new BaseDistributionMetric(2619,0,4,3747,0), result.getDistributionAtPosition(1));
		Assert.assertEquals(new BaseDistributionMetric(24,2557,3789,0,0), result.getDistributionAtPosition(3));
		Assert.assertEquals(new BaseDistributionMetric(2072,5,1487,2806,0), result.getDistributionAtPosition(6));
	}

	@Test
	public void writeOutput() {
		BaseDistributionAtReadPosition test= new BaseDistributionAtReadPosition();
		BaseDistributionMetricCollection result = test.gatherBaseQualities (TAGGED_FILE, "XC");
		File tempReportFile=getTempReportFile();

        result.writeOutput(tempReportFile);

		boolean testSuccess=false;
		try {
			testSuccess=FileUtils.contentEquals(tempReportFile, OUTPUT_EXPECTED_RESULT);
		} catch (IOException e) {
			e.printStackTrace();
		}

		Assert.assertTrue(testSuccess);
		tempReportFile.delete();

	}

	private File getTempReportFile () {
		File tempFile=null;

		try {
			tempFile = File.createTempFile("BaseDistributionAtReadPosition", "test_output");
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		return tempFile;
	}

	@Test
	public void testDoWork () {
		BaseDistributionAtReadPosition test= new BaseDistributionAtReadPosition();
		File tempReportFile=getTempReportFile();
		test.INPUT=TAGGED_FILE;
		test.TAG="XC";
		test.OUTPUT=tempReportFile;
		test.doWork();

		boolean testSuccess=false;
		try {
			testSuccess=FileUtils.contentEquals(tempReportFile, OUTPUT_EXPECTED_RESULT);
		} catch (IOException e) {
			e.printStackTrace();
		}

		Assert.assertTrue(testSuccess);
		tempReportFile.delete();

	}
}
