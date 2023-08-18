package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMFileHeader;
import junit.framework.Assert;

public class FilterReadsByUMISupportTest {
	
	private final static String rootDir="testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/";
	
	private final File INPUT_BAM = new File(rootDir, "filter_umi_test.sam");
	
	
	@Test
	public void testRetain2To4Reads() throws IOException {
		FilterReadsByUMISupport f = new FilterReadsByUMISupport();
		f.INPUT=Collections.singletonList(INPUT_BAM);
		f.MIN_READ_SUPPORT=2;
		f.MAX_READ_SUPPORT=4;
		f.OUTPUT=File.createTempFile("FilterReadsByUMISupportTest.", ".bam");
		f.OUTPUT.deleteOnExit();
		f.METRICS=File.createTempFile("FilterReadsByUMISupportTest.", ".metrics");
		f.METRICS.deleteOnExit();
		f.SORT_ORDER=SAMFileHeader.SortOrder.coordinate;
		f.customCommandLineValidation();
		
		int ret = f.doWork();
		Assert.assertTrue(ret==0);
		
		File EXPECTED_BAM=new File(rootDir, "filter_umi_test_min2_max4.bam");		
		TestUtils.assertSamFilesSame(f.OUTPUT, EXPECTED_BAM);

		File EXPECTED_METRICS=new File(rootDir, "filter_umi_test_min2_max4.metrics");
		TestUtils.testMetricsFilesEqual(EXPECTED_METRICS, f.METRICS);
		
	}
	
	@Test
	public void testRejectSingleReadUMIs() throws IOException {
		FilterReadsByUMISupport f = new FilterReadsByUMISupport();
		f.INPUT=Collections.singletonList(INPUT_BAM);
		f.MIN_READ_SUPPORT=2;
		f.OUTPUT=File.createTempFile("FilterReadsByUMISupportTest.", ".bam");
		f.OUTPUT.deleteOnExit();
		f.METRICS=File.createTempFile("FilterReadsByUMISupportTest.", ".metrics");
		f.METRICS.deleteOnExit();
		f.SORT_ORDER=SAMFileHeader.SortOrder.coordinate;
		f.customCommandLineValidation();
		
		int ret = f.doWork();
		Assert.assertTrue(ret==0);
		
		File EXPECTED_BAM=new File(rootDir, "filter_umi_test_min2.bam");		
		TestUtils.assertSamFilesSame(f.OUTPUT, EXPECTED_BAM);

		File EXPECTED_METRICS=new File(rootDir, "filter_umi_test_min2.metrics");
		TestUtils.testMetricsFilesEqual(EXPECTED_METRICS, f.METRICS);		
	}
	
	@Test
	public void testRetainOnlySingleReadUMIs() throws IOException {
		FilterReadsByUMISupport f = new FilterReadsByUMISupport();
		f.INPUT=Collections.singletonList(INPUT_BAM);
		f.MAX_READ_SUPPORT=1;
		f.OUTPUT=File.createTempFile("FilterReadsByUMISupportTest.", ".bam");
		f.OUTPUT.deleteOnExit();
		f.METRICS=File.createTempFile("FilterReadsByUMISupportTest.", ".metrics");
		f.METRICS.deleteOnExit();
		f.SORT_ORDER=SAMFileHeader.SortOrder.coordinate;
		f.customCommandLineValidation();
		
		int ret = f.doWork();
		Assert.assertTrue(ret==0);
		
		File EXPECTED_BAM=new File(rootDir, "filter_umi_test_max1.bam");		
		TestUtils.assertSamFilesSame(f.OUTPUT, EXPECTED_BAM);

		File EXPECTED_METRICS=new File(rootDir, "filter_umi_test_max1.metrics");
		TestUtils.testMetricsFilesEqual(EXPECTED_METRICS, f.METRICS);
		
	}
}
