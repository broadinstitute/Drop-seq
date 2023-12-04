package org.broadinstitute.dropseqrna.sbarro;

import java.io.File;
import java.io.IOException;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class TagReadWithRabiesBarcodesTest {

	
	/*
	 * Test data generation:
	 * gunzip -c /broad/mccarroll/arpy/sequence/rabies_barcoding/170915_170330_pSPBN_GFP_v9_v2_B19EnvA_15P_2B_day5_Final_RVg/170915_170330_pSPBN_GFP_v9_v2_B19EnvA_15P_2B_day5_Final_RVg_S1_L001_R1_001.fastq.gz |head -n 400 > test.fastq
	 * java -jar ~/jn_branch/3rdParty/picard/picard.jar FastqToSam F1=test.fastq SM=sample001 RG=rg0013 SO=queryname O=TagReadWithRabiesBarcodes_testdata.bam	
	 */
	
	
	private static final File TEST_BAM = new File("testdata/org/broadinstitute/dropseq/sbarro/TagReadWithRabiesBarcodes_testdata.bam");
	private static final File EXPECTED_REPORT = new File("testdata/org/broadinstitute/dropseq/sbarro/TagReadWithRabiesBarcodes_testdata.report.txt");
	private static final File EXPECTED_BAM = new File("testdata/org/broadinstitute/dropseq/sbarro/TagReadWithRabiesBarcodes_output.bam");
	
	
	/*
	 * Test data generation:
	 * 
	 * java -jar ~/jn_branch/3rdParty/picard/picard.jar FastqToSam 
	 * F1=/broad/mccarroll/arpy/sequence/rabies_barcoding/190717_SCC06_1e-2_A/190717_SCC06_1e-2_A_S1_L001_R1_001.fastq.gz
	 * F2=/broad/mccarroll/arpy/sequence/rabies_barcoding/190717_SCC06_1e-2_A/190717_SCC06_1e-2_A_S1_L001_R2_001.fastq.gz
	 * SM=sample001 RG=rg0013 SO=queryname O=Consensus_testdata.bam
	 * samtools view -s 0.00001 Consensus_testdata.bam > TagReadWithRabiesBarcodes_testdata_paired.bam
	 */
	
	private static final File TEST_BAM_CONSENSUS = new File("testdata/org/broadinstitute/dropseq/sbarro/TagReadWithRabiesBarcodes_testdata_paired.bam");
	private static final File EXPECTED_REPORT_CONSENSUS = new File("testdata/org/broadinstitute/dropseq/sbarro/TagReadWithRabiesBarcodes_testdata.consensus_report.txt");
	private static final File EXPECTED_BAM_CONSENSUS = new File("testdata/org/broadinstitute/dropseq/sbarro/TagReadWithRabiesBarcodes_output_consensus.bam");
	
	
	@Test
	public void doWorkTest() throws IOException {
		TagReadWithRabiesBarcodes t = new TagReadWithRabiesBarcodes();
		t.INPUT=TEST_BAM;
		t.OUTPUT=File.createTempFile("TagReadWithRabiesBarcodesTest", ".bam");
		t.OUTPUT.deleteOnExit();
		t.REPORT=File.createTempFile("TagReadWithRabiesBarcodesTest", ".report.txt");
		t.REPORT.deleteOnExit();
		t.GFP_ANCHOR_SEQUENCE="CGGCATGGACGAGCTGTACAAGTAAGCTA";  //the default, but important that the data uses the default.
		t.BASE_QUALITY_REPORT=File.createTempFile("TagReadWithRabiesBarcodesTest", ".bq_report.txt");
		
		int result = t.doWork();
		Assert.assertEquals(result, 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_REPORT, t.REPORT));
		
		TestUtils.assertSamFilesSame(EXPECTED_BAM, t.OUTPUT, false);
	}
		
	@Test
	public void doWorkConsensusTest() throws IOException {
		TagReadWithRabiesBarcodes t = new TagReadWithRabiesBarcodes();
		t.INPUT=TEST_BAM_CONSENSUS;
		t.OUTPUT=File.createTempFile("TagReadWithRabiesBarcodesTest", ".consensus.bam");
		t.OUTPUT.deleteOnExit();
		t.REPORT=File.createTempFile("TagReadWithRabiesBarcodesTest", ".consensus.report.txt");
		t.REPORT.deleteOnExit();
		t.GENERATE_CONSENSUS=true;
		t.GFP_ANCHOR_SEQUENCE="CGGCATGGACGAGCTGTACAAGTAAGCTA";  //the default, but important that the data uses the default.
		
		int result = t.doWork();
		Assert.assertEquals(result, 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_REPORT_CONSENSUS, t.REPORT));
		
		TestUtils.assertSamFilesSame(EXPECTED_BAM_CONSENSUS, t.OUTPUT, false);
	}
	
}
