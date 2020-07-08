package org.broadinstitute.dropseqrna.censusseq;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.broadinstitute.dropseqrna.censusseq.CsiAnalysis;
import org.broadinstitute.dropseqrna.censusseq.CsiMetrics;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class CsiAnalysisTest {

	private static final File IN_BAM = new File(
			"testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.bam");
	private static final File IN_VCF = new File(
			"testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.vcf.gz");
	private static final File IN_SAMPLE_LIST = new File(
			"testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.donorList.txt");
	private static final File IN_WRONG_SAMPLE_LIST = new File(
			"testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.wrong_donorList.txt");
	private static final File IN_INCOMPLETE_SAMPLE_LIST = new File(
			"testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.donorList.incomplete.txt");

	private static final File OUT_INCOMPLETE_SAMPLE_LIST = new File(
			"testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.FACS_poolAF_result.txt");

	private static final File OUT_INCOMPLETE_SAMPLE_LIST_1KG_AF = new File(
			"testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.FACS_1KG_AF_result.txt");

	@Test (enabled=true)
	// Using all the samples that generated the data, all variants are accounted
	// for, so there are no SNPs to examine.
	public void testCsi() throws IOException {
		CsiAnalysis f = new CsiAnalysis();
		f.INPUT_BAM = IN_BAM;
		f.INPUT_VCF = IN_VCF;
		f.OUTPUT = File.createTempFile("FACS.", ".txt");
		f.OUTPUT.deleteOnExit();
		f.SAMPLE_FILE = IN_SAMPLE_LIST;
		f.VCF_OUTPUT=File.createTempFile("FACS.", ".vcf.gz");
		f.VCF_OUTPUT.deleteOnExit();
		String TMP_DIR = f.OUTPUT.getParent();
		// TODO: what's the proper way to get the TMP DIR?
		f.TMP_DIR = Arrays.asList(new File(TMP_DIR));
		int ret = f.doWork();
		Assert.assertTrue(ret == 0);

	}

	@Test (enabled=true)
	//
	public void testCsiWrongSampleList() throws IOException {
		CsiAnalysis f = new CsiAnalysis();
		f.INPUT_BAM = IN_BAM;
		f.INPUT_VCF = IN_VCF;
		f.OUTPUT = File.createTempFile("FACS.", ".txt");
		f.OUTPUT.deleteOnExit();
		f.SAMPLE_FILE = IN_WRONG_SAMPLE_LIST;
		String TMP_DIR = f.OUTPUT.getParent();
		// TODO: what's the proper way to get the TMP DIR?
		f.TMP_DIR = Arrays.asList(new File(TMP_DIR));
		int ret = f.doWork();
		Assert.assertTrue(ret == 1);


	}

	@Test
	// By omitting samples from the list of donors, the data for those donors
	// becomes "contamination" of the sequencing data.
	public void testCsiIncompleteSampleList() throws IOException {
		CsiAnalysis f = new CsiAnalysis();
		f.INPUT_BAM = IN_BAM;
		f.INPUT_VCF = IN_VCF;
		f.OUTPUT = File.createTempFile("FACS.", ".txt");
		f.OUTPUT.deleteOnExit();
		f.SAMPLE_FILE = IN_INCOMPLETE_SAMPLE_LIST;
		f.KNOWN_DONOR_TAG = "ZS";
		String TMP_DIR = f.OUTPUT.getParent();
		// TODO: what's the proper way to get the TMP DIR?
		f.TMP_DIR = Arrays.asList(new File(TMP_DIR));
		int ret = f.doWork();
		Assert.assertTrue(ret == 0);
		Assert.assertTrue(TestUtils.testFilesSame(f.OUTPUT, OUT_INCOMPLETE_SAMPLE_LIST));

	}

	@Test
	// By omitting samples from the list of donors, the data for those donors
	// becomes "contamination" of the sequencing data.
	// use the thousand genomes allele freqs to estimate contamination.
	// What is the true contamination rate? Use R.
	// b=read.table("10_donors_chr22.selected_sites.readsPerDonor.txt", header=F,
	// stringsAsFactors=F)
	// a=read.table("8_donors_chr22.selected_sites.donorList.txt", header=F,
	// stringsAsFactors=F)
	// idx=match(a$V1, b$V2)
	public void testCsiIncompleteSampleList1KGAF() throws IOException {
		CsiAnalysis f = new CsiAnalysis();
		f.INPUT_BAM = IN_BAM;
		f.INPUT_VCF = IN_VCF;
		f.OUTPUT = File.createTempFile("FACS.", ".txt");
		f.ALLELE_FREQ_TAG = "KG_AF";
		f.KNOWN_DONOR_TAG="ZS";
		f.OUTPUT.deleteOnExit();
		f.SAMPLE_FILE = IN_INCOMPLETE_SAMPLE_LIST;
		String TMP_DIR = f.OUTPUT.getParent();
		// TODO: what's the proper way to get the TMP DIR?
		f.TMP_DIR = Arrays.asList(new File(TMP_DIR));		
		int ret = f.doWork();
		Assert.assertTrue(ret == 0);
		Assert.assertTrue(TestUtils.testFilesSame(f.OUTPUT, OUT_INCOMPLETE_SAMPLE_LIST_1KG_AF));
	}
	
	@Test (enabled=true)
	/**
	 * Instead of having an incomplete sample list, have the complete sample list, but specifically exclude a single sample from analysis.
	 * a=read.table("10_donors_chr22.selected_sites.donorList.txt", header=F)
	 * b=read.table("10_donors_chr22.selected_sites.donorList.incomplete.txt", header=F)
	 * setdiff (a$V1, b$V1)
     * [1] "CHB5_P25_140801"  "HUES74_P7_150201"
	 */
	public void testCsiExcludeSample1KGAF() throws IOException {
		CsiAnalysis f = new CsiAnalysis();
		f.INPUT_BAM = IN_BAM;
		f.INPUT_VCF = IN_VCF;
		f.OUTPUT = File.createTempFile("FACS.", ".txt");
		f.ALLELE_FREQ_TAG = "KG_AF";
		f.OUTPUT.deleteOnExit();
		f.SAMPLE_FILE = IN_INCOMPLETE_SAMPLE_LIST;
		f.KNOWN_DONOR_TAG="ZS";
		String TMP_DIR = f.OUTPUT.getParent();
		f.EXCLUDE_KNOWN_DONOR=Arrays.asList("CHB5_P25_140801", "HUES74_P7_150201");
		// TODO: what's the proper way to get the TMP DIR?
		f.TMP_DIR = Arrays.asList(new File(TMP_DIR));		
		int ret = f.doWork();
		Assert.assertTrue(ret == 0);
		Assert.assertTrue(TestUtils.testFilesSame(f.OUTPUT, OUT_INCOMPLETE_SAMPLE_LIST_1KG_AF));
	}
	
	
	
	/*
	@Test (enabled=false)
	public void testBQSimple () {
		CsiMetrics m = new CsiMetrics();
		m.ALT_COUNT=0;
		m.REF_COUNT=2000;
	
		m.ALT_DOSAGE=101;
		m.REF_DOSAGE=1899;
		for (int i=0; i<1000; i++) {
			m.addBaseErrorProbability(0.1);
		}
		for (int i=0; i<1000; i++)
			m.addBaseErrorProbability(0.001);
		
		double altFreqByCount = m.getAltFrequencyByCount();
		double altFreqByDosage = m.getAltFreqByDosage();
		double obsSeqErrorRate = m.getObservedSequencingErrorRate();
		
		Assert.assertEquals(altFreqByDosage, 0.0505, 0.001);
		Assert.assertEquals(altFreqByCount, 0, 0.001);
		Assert.assertEquals(altFreqByCount, 0, 0.001);
		Assert.assertEquals(obsSeqErrorRate, 0.0505, 0.001);
		
		// 0.03633333
				
	}
	*/

}
