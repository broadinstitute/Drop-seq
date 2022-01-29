/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;



public class DetectDoubletsTest {

	private static final File TESTDATA_DIR = new File("testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/sampleassignment/multisample");

	@Test (enabled=true)
	// I'm not sure why "big" is a much smaller data set than "small", but it is.
	public void testBig () throws IOException {
		DetectDoublets assigner = new DetectDoublets();
		File EXPECTED_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.txt.unguided_donors.txt");
		File EXPECTED_PAIR_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.unguided.perDonor.txt");
		File EXPECTED_PER_SNP_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.unguided.perSNP.txt.gz");

		// OUTPUT_PER_SNP
		assigner.INPUT_BAM=Collections.singletonList(new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG_retagged.bam"));
		assigner.VCF=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG.vcf");
		assigner.SINGLE_DONOR_LIKELIHOOD_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.txt");
		assigner.CELL_BC_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG_retagged.cellBarcodes.txt");
		assigner.SAMPLE_FILE=new File(TESTDATA_DIR, "donors.txt");
		assigner.OUTPUT= File.createTempFile("DetectDoublets", ".txt");
		assigner.OUTPUT_ALL_PAIRS=File.createTempFile("DetectDoublets", ".pairs.txt");
		assigner.OUTPUT_PER_SNP=File.createTempFile("DetectDoublets", ".per_snp.txt.gz");
		// assigner.USE_MISSING_DATA=false;
		assigner.OUTPUT.deleteOnExit();
		assigner.OUTPUT_ALL_PAIRS.deleteOnExit();
		assigner.OUTPUT_PER_SNP.deleteOnExit();

		assigner.FIXED_ERROR_RATE=0.1;
		assigner.GQ_THRESHOLD=30;
		int ret = assigner.doWork();
		Assert.assertTrue(ret==0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_OUTPUT, assigner.OUTPUT));
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_PAIR_OUTPUT, assigner.OUTPUT_ALL_PAIRS));
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_PER_SNP_OUTPUT, assigner.OUTPUT_PER_SNP));
	}
	
	@Test (enabled=true)
	// I'm not sure why "big" is a much smaller data set than "small", but it is.
	public void testBig2 () throws IOException {
		DetectDoublets assigner = new DetectDoublets();
		File EXPECTED_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.txt.unguided_donors.txt");
		File EXPECTED_PAIR_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.unguided.perDonor.txt");
		File EXPECTED_PER_SNP_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.unguided.perSNP.txt.gz");

		// OUTPUT_PER_SNP
		assigner.INPUT_BAM=Collections.singletonList(new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG_retagged.bam"));
		assigner.VCF=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG.vcf");
		assigner.SINGLE_DONOR_LIKELIHOOD_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.txt");
		assigner.CELL_BC_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG_retagged.cellBarcodes.txt");
		assigner.SAMPLE_FILE=new File(TESTDATA_DIR, "donors.txt");
		assigner.OUTPUT= File.createTempFile("DetectDoublets", ".txt");
		assigner.OUTPUT_ALL_PAIRS=File.createTempFile("DetectDoublets", ".pairs.txt");
		assigner.OUTPUT_PER_SNP=File.createTempFile("DetectDoublets", ".per_snp.txt.gz");
		// assigner.USE_MISSING_DATA=false;
		assigner.OUTPUT.deleteOnExit();
		assigner.OUTPUT_ALL_PAIRS.deleteOnExit();
		assigner.OUTPUT_PER_SNP.deleteOnExit();
		assigner.TEST_UNEXPECTED_DONORS=false;
		assigner.FIXED_ERROR_RATE=0.1;
		assigner.GQ_THRESHOLD=30;
		int ret = assigner.doWork();
		Assert.assertTrue(ret==0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_OUTPUT, assigner.OUTPUT));
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_PAIR_OUTPUT, assigner.OUTPUT_ALL_PAIRS));
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_PER_SNP_OUTPUT, assigner.OUTPUT_PER_SNP));
	}
	
	
	@Test (enabled=true)
	public void testBigWithContamination () throws IOException {
		DetectDoublets assigner = new DetectDoublets();
		File EXPECTED_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.contamination.sampleAssignments.txt");
		File EXPECTED_PAIR_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.contamination.sampleAssignments.perDonor.txt");
		File EXPECTED_PER_SNP_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.contamination.sampleAssignments.perSNP.txt.gz");

		// OUTPUT_PER_SNP
		assigner.INPUT_BAM=Collections.singletonList(new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG_retagged.bam"));
		assigner.VCF=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG.vcf");
		assigner.SINGLE_DONOR_LIKELIHOOD_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.txt");
		assigner.CELL_BC_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG_retagged.cellBarcodes.txt");
		assigner.SAMPLE_FILE=new File(TESTDATA_DIR, "donors.txt");
		assigner.ALLELE_FREQUENCY_ESTIMATE_FILE=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.minor_allele_freq.txt");;
		assigner.CELL_CONTAMINATION_ESTIMATE_FILE=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.contamination.txt");;
				
		assigner.OUTPUT= File.createTempFile("DetectDoublets", ".txt");
		assigner.OUTPUT_ALL_PAIRS=File.createTempFile("DetectDoublets", ".pairs.txt");
		assigner.OUTPUT_PER_SNP=File.createTempFile("DetectDoublets", ".per_snp.txt.gz");
		// assigner.USE_MISSING_DATA=false;
		assigner.OUTPUT.deleteOnExit();
		assigner.OUTPUT_ALL_PAIRS.deleteOnExit();
		assigner.OUTPUT_PER_SNP.deleteOnExit();

		assigner.FIXED_ERROR_RATE=0.1;
		assigner.GQ_THRESHOLD=30;
		int ret = assigner.doWork();
		Assert.assertTrue(ret==0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_OUTPUT, assigner.OUTPUT));
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_PAIR_OUTPUT, assigner.OUTPUT_ALL_PAIRS));
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_PER_SNP_OUTPUT, assigner.OUTPUT_PER_SNP));
	}
	
	@Test (enabled=true)
	public void testBigWithContaminationScaledLikelihood () throws IOException {
		DetectDoublets assigner = new DetectDoublets();
		File EXPECTED_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.contamination.sampleAssignments.txt");
		File EXPECTED_PAIR_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.contamination.sampleAssignments.perDonor.txt");
		File EXPECTED_PER_SNP_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.contamination.sampleAssignments.perSNP.txt.gz");

		// OUTPUT_PER_SNP
		assigner.INPUT_BAM=Collections.singletonList(new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG_retagged.bam"));
		assigner.VCF=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG.vcf");
		assigner.SINGLE_DONOR_LIKELIHOOD_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.txt");
		assigner.CELL_BC_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG_retagged.cellBarcodes.txt");
		assigner.SAMPLE_FILE=new File(TESTDATA_DIR, "donors.txt");
		assigner.ALLELE_FREQUENCY_ESTIMATE_FILE=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.minor_allele_freq.txt");;
		assigner.CELL_CONTAMINATION_ESTIMATE_FILE=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG.contamination.txt");;
				
		assigner.OUTPUT= File.createTempFile("DetectDoublets", ".txt");
		assigner.OUTPUT_ALL_PAIRS=File.createTempFile("DetectDoublets", ".pairs.txt");
		assigner.OUTPUT_PER_SNP=File.createTempFile("DetectDoublets", ".per_snp.txt.gz");
		// assigner.USE_MISSING_DATA=false;
		assigner.OUTPUT.deleteOnExit();
		assigner.OUTPUT_ALL_PAIRS.deleteOnExit();
		assigner.OUTPUT_PER_SNP.deleteOnExit();

		assigner.FIXED_ERROR_RATE=0.1;
		assigner.GQ_THRESHOLD=30;
		int ret = assigner.doWork();
		Assert.assertTrue(ret==0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_OUTPUT, assigner.OUTPUT));
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_PAIR_OUTPUT, assigner.OUTPUT_ALL_PAIRS));
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_PER_SNP_OUTPUT, assigner.OUTPUT_PER_SNP));
	}

	@Test (enabled=true)
	public void testSmall () throws IOException {
		File EXPECTED_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG2.sampleAssignments.txt.unguided_donors.txt");
		File EXPECTED_PAIR_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG2.sampleAssignments.unguided.perDonor.txt");

		DetectDoublets assigner = new DetectDoublets();
		assigner.INPUT_BAM=Collections.singletonList(new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2_retagged.bam"));
		assigner.VCF=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2.vcf");
		assigner.SINGLE_DONOR_LIKELIHOOD_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2.sampleAssignments.txt");
		assigner.CELL_BC_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2_retagged.cellBarcodes.txt");
		assigner.SAMPLE_FILE=new File(TESTDATA_DIR, "donors.txt");
		assigner.OUTPUT= File.createTempFile("DetectDoublets", ".txt");
		assigner.OUTPUT_ALL_PAIRS=File.createTempFile("DetectDoublets", ".pairs.txt");
		assigner.OUTPUT.deleteOnExit();
		assigner.SCALE_LIKELIHOODS=false;
		assigner.OUTPUT_ALL_PAIRS.deleteOnExit();
		// assigner.USE_MISSING_DATA=false;

		assigner.FIXED_ERROR_RATE=0.1;
		assigner.GQ_THRESHOLD=30;
		int ret = assigner.doWork();
		Assert.assertTrue(ret==0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_OUTPUT, assigner.OUTPUT));
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_PAIR_OUTPUT, assigner.OUTPUT_ALL_PAIRS));
	}
	
	@Test (enabled=true)
	public void testSmallScaledLikes () throws IOException {
		File EXPECTED_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG2.sampleAssignments.txt.unguided_donors.txt");
		File EXPECTED_PAIR_OUTPUT=new File(TESTDATA_DIR, "TEST_TTTGCGCGGAGC:ATTGTTTAGGAG2.sampleAssignments.unguided.scaledlikes.perDonor.txt");

		DetectDoublets assigner = new DetectDoublets();
		assigner.INPUT_BAM=Collections.singletonList(new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2_retagged.bam"));
		assigner.VCF=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2.vcf");
		assigner.SINGLE_DONOR_LIKELIHOOD_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2.sampleAssignments.txt");
		assigner.CELL_BC_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2_retagged.cellBarcodes.txt");
		assigner.SAMPLE_FILE=new File(TESTDATA_DIR, "donors.txt");
		assigner.OUTPUT= File.createTempFile("DetectDoublets", ".txt");
		assigner.OUTPUT_ALL_PAIRS=File.createTempFile("DetectDoublets", ".pairs.txt");
		assigner.OUTPUT.deleteOnExit();
		assigner.SCALE_LIKELIHOODS=true;
		assigner.OUTPUT_ALL_PAIRS.deleteOnExit();
		// assigner.USE_MISSING_DATA=false;

		assigner.FIXED_ERROR_RATE=0.1;
		assigner.GQ_THRESHOLD=30;
		int ret = assigner.doWork();
		Assert.assertTrue(ret==0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_OUTPUT, assigner.OUTPUT));
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_PAIR_OUTPUT, assigner.OUTPUT_ALL_PAIRS));
	}
	
	/**
	 * With split BAMs and small sets of selected cells, not finding any barcodes in a BAM is no longer a failure.
	 */
	@Test
	public void testNoBarcodesFound () throws IOException {
		File EXPECTED_OUTPUT=new File(TESTDATA_DIR, "empty.sampleAssignments.txt.unguided_donors.txt");
		File EXPECTED_PAIR_OUTPUT=new File(TESTDATA_DIR, "empty.sampleAssignments.unguided.perDonor.txt");

		DetectDoublets assigner = new DetectDoublets();
		assigner.INPUT_BAM=Collections.singletonList(new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2_retagged.bam"));
		assigner.VCF=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2.vcf");
		assigner.SINGLE_DONOR_LIKELIHOOD_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2.sampleAssignments.txt");
		assigner.CELL_BC_FILE=new File(TESTDATA_DIR, "matches_nothing.cellBarcodes.txt");
		assigner.SAMPLE_FILE=new File(TESTDATA_DIR, "donors.txt");
		assigner.OUTPUT= File.createTempFile("DetectDoublets", ".txt");
		assigner.OUTPUT_ALL_PAIRS=File.createTempFile("DetectDoublets", ".pairs.txt");
		assigner.OUTPUT.deleteOnExit();
		assigner.OUTPUT_ALL_PAIRS.deleteOnExit();
		// assigner.USE_MISSING_DATA=false;

		assigner.FIXED_ERROR_RATE=0.1;
		assigner.GQ_THRESHOLD=30;
		int ret = assigner.doWork();
		Assert.assertTrue(ret==0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_OUTPUT, assigner.OUTPUT));
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_PAIR_OUTPUT, assigner.OUTPUT_ALL_PAIRS));
	}

	@Test (enabled=true)
	public void testEmptyVcf () throws IOException {
		File EXPECTED_OUTPUT=new File(TESTDATA_DIR, "empty.sampleAssignments.txt.unguided_donors.txt");
		File EXPECTED_PAIR_OUTPUT=new File(TESTDATA_DIR, "empty.sampleAssignments.unguided.perDonor.txt");

		DetectDoublets assigner = new DetectDoublets();
		assigner.INPUT_BAM=Collections.singletonList(new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2_retagged.bam"));
		assigner.VCF=new File(TESTDATA_DIR, "empty.vcf");
		assigner.SINGLE_DONOR_LIKELIHOOD_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2.sampleAssignments.txt");
		assigner.CELL_BC_FILE=new File(TESTDATA_DIR, "TTTGCGCGGAGC:ATTGTTTAGGAG2_retagged.cellBarcodes.txt");
		assigner.SAMPLE_FILE=new File(TESTDATA_DIR, "donors.txt");
		assigner.OUTPUT= File.createTempFile("DetectDoublets", ".txt");
		assigner.OUTPUT_ALL_PAIRS=File.createTempFile("DetectDoublets", ".pairs.txt");
		assigner.OUTPUT.deleteOnExit();
		assigner.OUTPUT_ALL_PAIRS.deleteOnExit();
		// assigner.USE_MISSING_DATA=false;

		assigner.FIXED_ERROR_RATE=0.1;
		assigner.GQ_THRESHOLD=30;
		int ret = assigner.doWork();
		Assert.assertTrue(ret==1);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_OUTPUT, assigner.OUTPUT));
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_PAIR_OUTPUT, assigner.OUTPUT_ALL_PAIRS));
	}

	@Test (enabled=true)
	public void testFilterDonorMap() {
		Map<String, String> test = ImmutableMap.<String, String>builder()
				.put("k1", "v1")
				.put("k2", "v2")
			    .put("k3", "v3").build();
		
		Map<String, String> expected = ImmutableMap.<String, String>builder()
				.put("k1", "v1")
				.put("k2", "v2").build();
		
		List<String> list = Arrays.asList("v1", "v2");
		DetectDoublets dd = new DetectDoublets();
		Map<String,String> actual = dd.filterDonorMap(test, list);
		Assert.assertEquals(actual, expected);
		
		
	}
	
	

}
