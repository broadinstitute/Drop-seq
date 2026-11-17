/*
 * MIT License
 *
 * Copyright 2026 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.broadinstitute.dropseqrna.eqtl;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import picard.nio.PicardHtsPath;

/**
 * End-to-end test of PrepareEqtlData using small existing fixtures: 38 donors x 100 genes of meta cell expression,
 * and a 5-variant chr22 VCF containing those donors.  Covariates are generated to match the meta cell donors.
 */
public class PrepareEqtlDataTest {

	private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/eqtl/");
	private static final File META_CELL = new File(TEST_DATA_DIR, "100_genes.38_donors.meta_cells.txt");
	private static final File GTF = new File(TEST_DATA_DIR, "100_genes.48_donors.gtf");
	private static final File SD_FILE = new File(TEST_DATA_DIR, "hg19.dict");
	private static final File VCF = new File(TEST_DATA_DIR, "38_donors.chr1_10_22.vcf.gz");

	private static List<String> metaCellDonors() throws IOException {
		final String header = Files.readAllLines(META_CELL.toPath()).get(0);
		final String[] f = header.split("\t");
		return new ArrayList<>(Arrays.asList(f).subList(1, f.length));
	}

	private static File writeCovariates(final List<String> donors) throws IOException {
		final File f = File.createTempFile("covariates", ".txt");
		f.deleteOnExit();
		final StringBuilder h = new StringBuilder("id");
		final StringBuilder sex = new StringBuilder("SEX");
		final StringBuilder age = new StringBuilder("AGE");
		for (int i = 0; i < donors.size(); i++) {
			h.append('\t').append(donors.get(i));
			sex.append('\t').append(1 + i % 2);
			age.append('\t').append(20 + i);
		}
		Files.write(f.toPath(), Arrays.asList(h.toString(), sex.toString(), age.toString()));
		return f;
	}

	private static File temp(final String prefix) throws IOException {
		final File f = File.createTempFile(prefix, ".txt");
		f.deleteOnExit();
		return f;
	}

	private PrepareEqtlData makeBase(final List<String> covariateDonors) throws IOException {
		final PrepareEqtlData p = new PrepareEqtlData();
		p.META_CELL_FILE = Arrays.asList(META_CELL);
		p.COVARIATE_FILE = Arrays.asList(writeCovariates(covariateDonors));
		p.INPUT_VCF = new PicardHtsPath(VCF);
		p.ANNOTATIONS_FILE = GTF;
		p.SEQUENCE_DICTIONARY = SD_FILE;
		p.REMOVE_PCT_EXPRESSION = 0d;
		p.MAF = 0.01;
		p.HWE_PVALUE = 0d;
		p.COVARIATE_MATRIX = temp("covariate_matrix");
		p.EXPRESSION_FILE = temp("expression");
		p.GENE_LOCATION_FILE = temp("gene_location");
		p.GENOTYPE_MATRIX = temp("genotypes");
		p.SNP_LOCATIONS = temp("snps");
		p.OUT_DONOR_LIST = temp("donors");
		return p;
	}

	private static final File EXPECTED_DIR = new File(TEST_DATA_DIR, "PrepareEqtlData");

	private static void assertMatchesExpected(final String expectedName, final File actual) throws IOException {
		final File expected = new File(EXPECTED_DIR, expectedName);
		if (System.getenv("EQTL_REGEN") != null) {
			Files.createDirectories(EXPECTED_DIR.toPath());
			Files.copy(actual.toPath(), expected.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
		}
		Assert.assertTrue(TestUtils.testFilesSame(expected, actual), "Output differs from " + expected);
	}

	private static String[] header(final File f) throws IOException {
		return Files.readAllLines(f.toPath()).get(0).split("\t");
	}

	@Test
	public void testEndToEnd() throws IOException {
		final List<String> donors = metaCellDonors();
		final PrepareEqtlData p = makeBase(donors);
		p.ADD_BATCH_COVARIATE = true;
		p.SNP_GENE_MAP = temp("snp_gene_map");
		p.CIS_DIST = 1_000_000;
		p.NS_CIS_DIST = 1_000_000;
		Assert.assertNull(p.customCommandLineValidation());
		Assert.assertEquals(p.doWork(), 0);

		// Expression and covariate matrices must have the same donors in the same order.
		final String[] expressionHeader = header(p.EXPRESSION_FILE);
		final String[] covariateHeader = header(p.COVARIATE_MATRIX);
		Assert.assertEquals(expressionHeader.length - 1, donors.size());
		Assert.assertEquals(
				Arrays.asList(covariateHeader).subList(1, covariateHeader.length),
				Arrays.asList(expressionHeader).subList(1, expressionHeader.length));
		Assert.assertEquals(Files.readAllLines(p.OUT_DONOR_LIST.toPath()), donors);

		// SEX, AGE and the batch covariate.
		Assert.assertTrue(Files.readAllLines(p.COVARIATE_MATRIX.toPath()).size() >= 3);
		// Header plus at least one variant passed filtering, and genotypes match snp locations.
		final List<String> snps = Files.readAllLines(p.SNP_LOCATIONS.toPath());
		final List<String> genotypes = Files.readAllLines(p.GENOTYPE_MATRIX.toPath());
		Assert.assertTrue(snps.size() > 1);
		Assert.assertEquals(genotypes.size(), snps.size());
		Assert.assertEquals(genotypes.get(0).split("\t").length - 1, donors.size());
		Assert.assertTrue(Files.readAllLines(p.SNP_GENE_MAP.toPath()).size() > 1);

		// Content must not change.
		assertMatchesExpected("endToEnd.expression.txt", p.EXPRESSION_FILE);
		assertMatchesExpected("endToEnd.gene_location.txt", p.GENE_LOCATION_FILE);
		assertMatchesExpected("endToEnd.covariates.txt", p.COVARIATE_MATRIX);
		assertMatchesExpected("endToEnd.genotypes.txt", p.GENOTYPE_MATRIX);
		assertMatchesExpected("endToEnd.snp_locations.txt", p.SNP_LOCATIONS);
		assertMatchesExpected("endToEnd.donors.txt", p.OUT_DONOR_LIST);
		assertMatchesExpected("endToEnd.snp_gene_map.txt", p.SNP_GENE_MAP);
	}

	@Test
	public void testRejectedDonors() throws IOException {
		final List<String> donors = metaCellDonors();
		final PrepareEqtlData p = makeBase(donors);
		final File rejected = temp("rejected");
		Files.write(rejected.toPath(), donors.subList(0, 2));
		p.REJECTED_DONOR_LIST = rejected;
		p.GENOTYPE_BED = temp("genotypes_bed");
		p.PLINK_TPED = File.createTempFile("plink", ".tped");
		p.PLINK_TPED.deleteOnExit();
		p.EXPRESSION_BED_FILE = temp("expression_bed");
		Assert.assertEquals(p.doWork(), 0);
		Assert.assertEquals(header(p.EXPRESSION_FILE).length - 1, donors.size() - 2);
		Assert.assertEquals(header(p.COVARIATE_MATRIX).length - 1, donors.size() - 2);

		assertMatchesExpected("rejected.expression.txt", p.EXPRESSION_FILE);
		assertMatchesExpected("rejected.expression.bed", p.EXPRESSION_BED_FILE);
		assertMatchesExpected("rejected.covariates.txt", p.COVARIATE_MATRIX);
		assertMatchesExpected("rejected.genotypes.bed", p.GENOTYPE_BED);
		assertMatchesExpected("rejected.tped", p.PLINK_TPED);
		assertMatchesExpected("rejected.tfam", new File(p.PLINK_TPED.getPath().replaceAll("\\.tped$", ".tfam")));
	}

	@Test
	public void testMismatchedCovariateDonorsFailsValidation() throws IOException {
		final List<String> donors = metaCellDonors();
		final PrepareEqtlData p = makeBase(donors.subList(0, donors.size() - 1));
		Assert.assertEquals(p.doWork(), 1);
	}

	@Test
	public void testSnpGeneMapValidationErrors() throws IOException {
		final PrepareEqtlData p = makeBase(metaCellDonors());
		// SNP_GENE_MAP requested without any distance or interval criteria.
		p.SNP_GENE_MAP = temp("snp_gene_map");
		Assert.assertNotNull(p.customCommandLineValidation());
	}

	@Test
	public void testGetSequenceDictionaryFile() {
		Assert.assertEquals(PrepareEqtlData.getSequenceDictionaryFile(SD_FILE, GTF), SD_FILE);
		Assert.assertThrows(RuntimeException.class,
				() -> PrepareEqtlData.getSequenceDictionaryFile(new File("does_not_exist.dict"), GTF));
		// 100_genes.48_donors.dict does not exist next to the GTF.
		Assert.assertThrows(RuntimeException.class, () -> PrepareEqtlData.getSequenceDictionaryFile(null, GTF));
		Assert.assertEquals(
				PrepareEqtlData.getSequenceDictionaryFile(null, new File(TEST_DATA_DIR, "hg19.gtf.gz")), SD_FILE);
	}
}
