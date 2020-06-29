package org.broadinstitute.dropseqrna.censusseq;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.annotations.Test;

import junit.framework.Assert;

public class FindDonorRatiosInBulkDataTest {
	private static final File IN_BAM = new File("testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.bam");
	private static final File IN_VCF = new File("testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.vcf.gz");
	private static final File IN_SAMPLE_LIST = new File("testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.donorList.txt");
	private static final File IN_WRONG_SAMPLE_LIST = new File("testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.wrong_donorList.txt");


	private static final File OUT_ROLL_CALL = new File("testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.roll_call.txt");
	private static final File OUT_ROLL_CALL_VERBOSE = new File("testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.roll_call.verbose.txt.gz");

	private static final File OUT_CENSUS = new File("testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.census.txt");
	private static final File OUT_CENSUS_SNP_HISTOGRAM = new File("testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.census.snp_histogram.txt");
	@Test
	// Tests full path and result files.  Barebones but useful.
	// the math is checked more stringently in other unit tests.
	public void testRollCall() throws IOException {
		FindDonorRatiosInBulkData f = new FindDonorRatiosInBulkData();
		f.INPUT_BAM=IN_BAM;
		f.INPUT_VCF=IN_VCF;
		f.KNOWN_DONOR_TAG="ZS";
		f.MIN_BASE_QUALITY=null;
		f.OUTPUT=File.createTempFile("testRollCall.", ".roll_call.txt");
		f.OUTPUT.deleteOnExit();
		f.COMMON_SNP_ANALYSIS=false;
		f.OUTPUT_VERBOSE=File.createTempFile("testRollCall.", ".roll_call.verbose.txt.gz");
		f.OUTPUT_VERBOSE.deleteOnExit();
		String TMP_DIR=f.OUTPUT.getParent();
		//TODO: what's the proper way to get the TMP DIR?
		f.TMP_DIR=Arrays.asList(new File (TMP_DIR));
		int ret = f.doWork();
		Assert.assertTrue(ret==0);
		Assert.assertTrue(TestUtils.testFilesSame(OUT_ROLL_CALL, f.OUTPUT));
		Assert.assertTrue(TestUtils.testFilesSame(OUT_ROLL_CALL_VERBOSE, f.OUTPUT_VERBOSE));

	}

	@Test
	// Tests full path and result files.  Barebones but useful.
	// the math is checked more stringently in other unit tests.
	public void testCensusCorrectSampleFile() throws IOException {
		FindDonorRatiosInBulkData f = new FindDonorRatiosInBulkData();
		f.INPUT_BAM=IN_BAM;
		f.INPUT_VCF=IN_VCF;
		f.KNOWN_DONOR_TAG="ZS";
		f.MIN_BASE_QUALITY=null;
		f.OUTPUT=File.createTempFile("testCensus.", ".census.txt");
		f.OUTPUT.deleteOnExit();
		f.COMMON_SNP_ANALYSIS=true;
		f.SAMPLE_FILE=IN_SAMPLE_LIST;
		f.REPORT_ALLELE_COUNTS=true;
		f.SNP_COVERAGE_HISTOGRAM=File.createTempFile("testCensus.", ".snp_histogram.txt");
		f.SNP_COVERAGE_HISTOGRAM.deleteOnExit();
		String TMP_DIR=f.OUTPUT.getParent();
		//TODO: what's the proper way to get the TMP DIR?
		f.TMP_DIR=Arrays.asList(new File (TMP_DIR));
		int ret = f.doWork();
		Assert.assertTrue(ret==0);
		Assert.assertTrue(TestUtils.testFilesSame(OUT_CENSUS, f.OUTPUT));
		Assert.assertTrue(TestUtils.testFilesSame(OUT_CENSUS_SNP_HISTOGRAM, f.SNP_COVERAGE_HISTOGRAM));
	}

	@Test
	public void testCensusWrongSampleFile() throws IOException {
		FindDonorRatiosInBulkData f = new FindDonorRatiosInBulkData();
		f.INPUT_BAM=IN_BAM;
		f.INPUT_VCF=IN_VCF;
		f.KNOWN_DONOR_TAG="ZS";
		f.MIN_BASE_QUALITY=null;
		f.OUTPUT=File.createTempFile("testCensus.", ".census.txt");
		f.OUTPUT.deleteOnExit();
		f.COMMON_SNP_ANALYSIS=true;
		f.SAMPLE_FILE=IN_WRONG_SAMPLE_LIST;
		String TMP_DIR=f.OUTPUT.getParent();
		//TODO: what's the proper way to get the TMP DIR?
		f.TMP_DIR=Arrays.asList(new File (TMP_DIR));
		int ret = f.doWork();
		// given a donor list that doesn't completely overlap the donors in the VCF, the program errors and exits.
		Assert.assertTrue(ret==1);

	}





}
