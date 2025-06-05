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

package org.broadinstitute.dropseqrna.eqtl;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.dropseqrna.eqtl.DonorCovariates;
import org.testng.Assert;
import org.testng.annotations.Test;

public class DonorCovariatesTest {

	private final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/eqtl/");
	private final File COVARS1 = new File (TEST_DATA_DIR, "covars_test1.txt");
	private final File COVARS2 = new File (TEST_DATA_DIR, "covars_test2.txt");
	private final File COVARS3 = new File (TEST_DATA_DIR, "covars_test3.txt");
	
	@Test
	public void merge() {
		final DonorCovariates result = DonorCovariates.parseFile(COVARS1, ValidationStringency.SILENT);
		final DonorCovariates other = DonorCovariates.parseFile(COVARS2, ValidationStringency.SILENT);
		result.merge(other);
		// contains 6 donors
		Assert.assertEquals(5, result.getDonorNames().size());
		// 5 phenotypes
		Assert.assertEquals(5, result.getCovariates().size());
		// has a value from both files.
		Assert.assertEquals(result.getValue("CSES12_P14_140616", "AGE"), "8");
		Assert.assertEquals(result.getValue("CSES12_P14_140616", "NUM_CELLS"), "44");
		// has no value because donors are disjoint
		Assert.assertNull(result.getValue("CT4_P15_140622", "SEX"));
		Assert.assertNull(result.getValue("CT2_P11_140622", "AVERAGE_UMIS_PER_CELL"));
	}

	@Test
	public void parseFile() {
		final DonorCovariates result = DonorCovariates.parseFile(COVARS1, ValidationStringency.SILENT);
		Assert.assertNotNull(result);
		// donors are sorted alphanumeric by set.
		Set<String> donors = new HashSet<> (Arrays.asList("CHB5_P25_140801", "CHB9_P24_140728", "CSES12_P14_140616", "CT2_P11_140622"));
		Assert.assertEquals(new HashSet<String>(result.getDonorNames()), donors);
		Assert.assertEquals(result.getCovariates(), new HashSet<String>(Arrays.asList("SEX", "AGE")));
		// random values
		Assert.assertEquals(result.getValue("CSES12_P14_140616", "AGE"), "8");
	}

	@Test
	public void parseFileLenient() {
		final DonorCovariates result = DonorCovariates.parseFile(COVARS3, ValidationStringency.LENIENT);
		Assert.assertNotNull(result);
		final List<String> donors =
				Arrays.asList("CHB5_P25_140801", "CHB9_P24_140728", "CSES12_P14_140616", "CT2_P11_140622");
		Assert.assertEqualsNoOrder(result.getDonorNames(), donors);
		Assert.assertEqualsNoOrder(result.getCovariates(), Arrays.asList("SEX", "AGE"));
		Assert.assertEquals(result.getValues("SEX", donors), Arrays.asList("2", "NA", "2", "1"));
		Assert.assertEquals(result.getValues("AGE", donors), Arrays.asList("5", "5", "8", "NA"));
	}

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void parseFileStrict() {
		DonorCovariates.parseFile(COVARS3, ValidationStringency.STRICT);
	}
}
