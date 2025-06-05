package org.broadinstitute.dropseqrna.eqtl;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.eqtl.PrepareEqtlGenotypeData;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class PrepareEqtlGenotypeDataTest {

	private final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/eqtl/");
	private final File INPUT_VCF = new File(TEST_DATA_DIR, "integration.vcf.gz");
	private final File DONOR_LIST_FILE = new File(TEST_DATA_DIR, "heat1-6.donor_list.txt");

	private final File THREE_DONOR_FILE=new File(TEST_DATA_DIR, "3_samples.txt");
	private final File THREE_DONOR_EXPECTED_GENOTYPE_MATRIX =
			new File(TEST_DATA_DIR, "3_samples.genotypes.txt");
	private final File THREE_DONOR_EXPECTED_SNPs =
			new File(TEST_DATA_DIR, "3_samples.snp_locations.txt");
	private final File THREE_DONOR_EXPECTED_GENOTYPE_BED =
			new File(TEST_DATA_DIR, "3_samples.genotypes.bed");
	private final File THREE_DONOR_EXPECTED_TPED = new File(TEST_DATA_DIR, "3_samples.tped");
	private final File THREE_DONOR_EXPECTED_TFAM = new File(TEST_DATA_DIR, "3_samples.tfam");
	private final File THREE_DONOR_EXPECTED_REF_ALLELE = new File(TEST_DATA_DIR, "3_samples.ref_allele");

	private final File ATAC_SEQ_INTERVALS = new File(TEST_DATA_DIR, "atac_seq.intervals");
	private final File GTF = new File(TEST_DATA_DIR, "atac_seq_test.gtf");

	private final File EXPECTED_SNPS_GTF_FILTERED=new File(TEST_DATA_DIR, "3_samples_gtf_filtered.snp_locations.txt");
	private final File EXPECTED_SNPS_INTERVAL_FILTERED=new File(TEST_DATA_DIR, "3_samples_atac_seq_filtered.snp_locations.txt");
	private final File EXPECTED_SNPS_BOTH_FILTERED=new File(TEST_DATA_DIR, "3_samples_both_filtered.snp_locations.txt");


	@Test
	public void testFilteringIntegration () {
		List<String> requestedDonorList =ParseBarcodeFile.readCellBarcodeFile(this.DONOR_LIST_FILE);
		// dHefine some thresholds.
		double maf=0.1;
		double hwe=0.1;

		final VCFFileReader vcfReader = new VCFFileReader(this.INPUT_VCF, false);

		PeekableIterator<VariantContext> iter =
				new PrepareEqtlGenotypeData().getVCFIterator(
						vcfReader, null, requestedDonorList, maf, hwe, 30, 0.90, null, false);

		Set<Integer> expectedSites = new HashSet<> (Arrays.asList(16051347,16051497));
		Set <Integer> observedSites = new HashSet<>();
		while (iter.hasNext()) {
			VariantContext vc = iter.next();
			observedSites.add(vc.getStart());
		}

		Assert.assertEquals(observedSites, expectedSites);
		iter.close();
		CloserUtil.close(vcfReader);
	}

	// validate that the dosage of the alt allele is being encoded in the output genotype matrix.
	// For the VCF: GT genotype, encoded as alleles values separated by either of / or |, e.g. The allele values are 0 for the reference allele (what is in the reference sequence), 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on.
	//
	// EX:  bcftools view -r 22:16050822 -s CHB10_P24_140728,UCLA9_P16_150821,CHB11_P25_140722 integration.vcf.gz |cut -f 1-5,10-
	// #CHROM	POS	ID	REF	ALT	CHB10_P24_140728	UCLA9_P16_150821	CHB11_P25_140722
	// 22	16050822	rs190022386	G	A	0/0:32,0:32:96:0,96,933	1/1:1,37:38:99:1168,104,0	0/1:44,20:64:99:450,0,1027
	@Test
	public void fullTestGenotypeMatrix() throws IOException {
		PrepareEqtlGenotypeData p = new PrepareEqtlGenotypeData();
		p.INPUT_VCF = new PicardHtsPath(INPUT_VCF);
		p.SAMPLE_FILE=THREE_DONOR_FILE;
		p.GENOTYPE_MATRIX=File.createTempFile("PrepareEqtlGenotypeData.", ".genotypes.txt");
		p.SNP_LOCATIONS=File.createTempFile("PrepareEqtlGenotypeData.", ".snp_locations.txt");
		p.GENOTYPE_MATRIX.deleteOnExit();
		p.SNP_LOCATIONS.deleteOnExit();
		int result = p.doWork();
		Assert.assertEquals(result, 0);
		Assert.assertTrue(TestUtils.testFilesSame(THREE_DONOR_EXPECTED_GENOTYPE_MATRIX, p.GENOTYPE_MATRIX));
		Assert.assertTrue(TestUtils.testFilesSame(THREE_DONOR_EXPECTED_SNPs, p.SNP_LOCATIONS));
	}

	@Test
	public void fullTestGenotypeBed() throws IOException {
		PrepareEqtlGenotypeData p = new PrepareEqtlGenotypeData();
		p.INPUT_VCF = new PicardHtsPath(INPUT_VCF);
		p.SAMPLE_FILE=THREE_DONOR_FILE;
		p.GENOTYPE_BED=File.createTempFile("PrepareEqtlGenotypeData.", ".genotypes.bed");
		p.GENOTYPE_BED.deleteOnExit();
		int result = p.doWork();
		Assert.assertEquals(result, 0);
		Assert.assertTrue(TestUtils.testFilesSame(THREE_DONOR_EXPECTED_GENOTYPE_BED, p.GENOTYPE_BED));
	}

	/**
	 * 	Validate that the dosage of the alt allele is being encoded in the output plink tped.
	 * 	For the VCF: GT genotype, encoded as alleles values separated by either of / or |, e.g.
	 * 	The allele values are 0 for the reference allele (what is in the reference sequence),
	 * 	1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on.
	 *
	 * 	EX:  bcftools view -r 22:16050822 -s CHB10_P24_140728,UCLA9_P16_150821,CHB11_P25_140722 integration.vcf.gz |cut -f 1-5,10-
	 * 	#CHROM	POS	ID	REF	ALT	CHB10_P24_140728	UCLA9_P16_150821	CHB11_P25_140722
	 * 	22	16050822	rs190022386	G	A	0/0:32,0:32:96:0,96,933	1/1:1,37:38:99:1168,104,0	0/1:44,20:64:99:450,0,1027
	 */
	@Test
	public void fullTestPlinkTped() throws IOException {
		final PrepareEqtlGenotypeData p = new PrepareEqtlGenotypeData();
		final File tped = File.createTempFile(
				"PrepareEqtlGenotypeData.", PrepareEqtlGenotypeData.TPED_FILE_EXTENSION);
		final File tfam = PrepareEqtlGenotypeData.resolvePlinkSibling(
				tped, PrepareEqtlGenotypeData.TFAM_FILE_EXTENSION);
		final File refAllele = PrepareEqtlGenotypeData.resolvePlinkSibling(
				tped, PrepareEqtlGenotypeData.REF_ALLELE_FILE_EXTENSION);
		p.INPUT_VCF = new PicardHtsPath(INPUT_VCF);
		p.SAMPLE_FILE=THREE_DONOR_FILE;
		p.PLINK_TPED=tped;
		tped.deleteOnExit();
		tfam.deleteOnExit();
		refAllele.deleteOnExit();
		int result = p.doWork();
		Assert.assertEquals(result, 0);
		Assert.assertTrue(TestUtils.testFilesSame(THREE_DONOR_EXPECTED_TPED, tped));
		Assert.assertTrue(TestUtils.testFilesSame(THREE_DONOR_EXPECTED_TFAM, tfam));
		Assert.assertTrue(TestUtils.testFilesSame(THREE_DONOR_EXPECTED_REF_ALLELE, refAllele));
	}

	/**
	 * 	Validate that the dosage of the alt allele is being encoded in the output genotype matrix and
	 * 	plink tped.
	 * 	For the VCF: GT genotype, encoded as alleles values separated by either of / or |, e.g.
	 * 	The allele values are 0 for the reference allele (what is in the reference sequence),
	 * 	1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on.
	 *
	 * 	EX:  bcftools view -r 22:16050822 -s CHB10_P24_140728,UCLA9_P16_150821,CHB11_P25_140722 integration.vcf.gz |cut -f 1-5,10-
	 * 	#CHROM	POS	ID	REF	ALT	CHB10_P24_140728	UCLA9_P16_150821	CHB11_P25_140722
	 * 	22	16050822	rs190022386	G	A	0/0:32,0:32:96:0,96,933	1/1:1,37:38:99:1168,104,0	0/1:44,20:64:99:450,0,1027
	 */
	@Test
	public void fullTestGenotypeMatrixBedAndTped() throws IOException {
		PrepareEqtlGenotypeData p = new PrepareEqtlGenotypeData();

		final File tped = File.createTempFile(
				"PrepareEqtlGenotypeData.", PrepareEqtlGenotypeData.TPED_FILE_EXTENSION);
		final File tfam = PrepareEqtlGenotypeData.resolvePlinkSibling(
				tped, PrepareEqtlGenotypeData.TFAM_FILE_EXTENSION);
		final File refAllele = PrepareEqtlGenotypeData.resolvePlinkSibling(
				tped, PrepareEqtlGenotypeData.REF_ALLELE_FILE_EXTENSION);

		p.INPUT_VCF = new PicardHtsPath(INPUT_VCF);
		p.SAMPLE_FILE=THREE_DONOR_FILE;
		p.GENOTYPE_MATRIX = File.createTempFile("PrepareEqtlGenotypeData.", ".genotypes.txt");
		p.SNP_LOCATIONS = File.createTempFile("PrepareEqtlGenotypeData.", ".snp_locations.txt");
		p.GENOTYPE_BED = File.createTempFile("PrepareEqtlGenotypeData.", ".genotypes.bed");
		p.GENOTYPE_MATRIX.deleteOnExit();
		p.SNP_LOCATIONS.deleteOnExit();
		p.GENOTYPE_BED.deleteOnExit();
		p.PLINK_TPED=tped;
		tped.deleteOnExit();
		tfam.deleteOnExit();
		refAllele.deleteOnExit();
		int result = p.doWork();
		Assert.assertEquals(result, 0);
		Assert.assertTrue(TestUtils.testFilesSame(THREE_DONOR_EXPECTED_GENOTYPE_MATRIX, p.GENOTYPE_MATRIX));
		Assert.assertTrue(TestUtils.testFilesSame(THREE_DONOR_EXPECTED_SNPs, p.SNP_LOCATIONS));
		Assert.assertTrue(TestUtils.testFilesSame(THREE_DONOR_EXPECTED_GENOTYPE_BED, p.GENOTYPE_BED));
		Assert.assertTrue(TestUtils.testFilesSame(THREE_DONOR_EXPECTED_TPED, tped));
		Assert.assertTrue(TestUtils.testFilesSame(THREE_DONOR_EXPECTED_TFAM, tfam));
		Assert.assertTrue(TestUtils.testFilesSame(THREE_DONOR_EXPECTED_REF_ALLELE, refAllele));
	}

	@Test
	public void fullTestWithIntervals() throws IOException {
		PrepareEqtlGenotypeData p = new PrepareEqtlGenotypeData();
		p.INPUT_VCF = new PicardHtsPath(INPUT_VCF);
		p.GQ_THRESHOLD = 0;
		p.HWE_PVALUE = 1e-300;
		p.MAF = 0;
		p.FRACTION_SAMPLES_PASSING = 0;
		p.SAMPLE_FILE = THREE_DONOR_FILE;
		p.GENOTYPE_MATRIX = File.createTempFile("PrepareEqtlGenotypeData.", ".genotypes.txt");
		p.SNP_LOCATIONS = File.createTempFile("PrepareEqtlGenotypeData.", ".snp_locations.txt");
		p.GENOTYPE_MATRIX.deleteOnExit();
		p.SNP_LOCATIONS.deleteOnExit();
		p.INTERVAL_FILE = this.ATAC_SEQ_INTERVALS;
		int result = p.doWork();
		Assert.assertEquals(result, 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_SNPS_INTERVAL_FILTERED, p.SNP_LOCATIONS));
	}
}
