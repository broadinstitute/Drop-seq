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
import java.io.IOException;
import java.util.*;
import java.util.stream.Stream;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.eqtl.EqtlCovariate;
import org.testng.Assert;
import org.testng.annotations.Test;

import picard.PicardException;

public class EqtlCovariateTest {
	private final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/eqtl/");
	private final File validCovars = new File(TEST_DATA_DIR, "valid.covariates.txt");
	private final File invalidCovars = new File(TEST_DATA_DIR, "invalid.covariates.txt");
	private final File missingCovars = new File(TEST_DATA_DIR, "missing.covariates.txt");
	private final File nonnumericCovars =
			new File(TEST_DATA_DIR, "nonnumeric.covariates.txt");

	@Test
	public void testMerge () {
		List<String> donors = Arrays.asList("A", "B", "C");
		List<String> donors2 = Arrays.asList("D", "E", "F");

		EqtlCovariate d = new EqtlCovariate(donors);
		String [] values1 = {"1","2","3"};
		d.setValues("attribute1", values1);
		d.toString();

		EqtlCovariate d2 = new EqtlCovariate(donors2);

		String [] values2 = {"4","5","6"};
		d2.setValues("attribute1", values2);

		d.mergeDonors(d2);
		// test donor size.
		Assert.assertEquals(6, d.donorNames().size());
		// test number of attributes
		Assert.assertEquals(1, d.getAttributes().size());

		String[] expectedValues = Stream.of(values1, values2).flatMap(Stream::of).toArray(String[]::new);

		String [] resultValues = d.getValues("attribute1");
		Assert.assertEquals(expectedValues, resultValues);
	}

	@Test(expectedExceptions=TranscriptomeException.class)
	// fail because attributes are different.
	public void testFailedMerge1 () {
		List<String> donors = Arrays.asList("A", "B", "C");
		List<String> donors2 = Arrays.asList("D", "E", "F");

		EqtlCovariate d = new EqtlCovariate(donors);
		String [] values1 = {"1","2","3"};
		d.setValues("attribute1", values1);
		d.toString();

		EqtlCovariate d2 = new EqtlCovariate(donors2);

		String [] values2 = {"4","5","6"};
		d2.setValues("attribute2", values2);

		boolean safeToMerge=d.validateMergeDonors(d2);
		Assert.assertFalse(safeToMerge);
		d.mergeDonors(d2);

	}

	@Test
	// test read/write.
	public void testRoundTrip () {
		File f=null;
		try {
			f = File.createTempFile("tmp_covariate", ".txt");
			f.deleteOnExit();
		} catch (IOException e) {
			e.printStackTrace();
		}

		List<String> donors = Arrays.asList("A", "B", "C");
		EqtlCovariate d = new EqtlCovariate(donors);
		String [] values1 = {"1","2","3"};
		d.setValues("attribute1", values1);
		String [] values2 = {"4","5","6"};
		d.setValues("attribute2", values2);

		d.writeFile(f);
		EqtlCovariate d2=EqtlCovariate.parseFile(f);
		for (String a: d.getAttributes()) {
			String [] v1 = d.getValues(a);
			String [] v2 = d2.getValues(a);
			Assert.assertEquals(v1, v2);
		}

	}

	@Test(expectedExceptions=PicardException.class)
	public void testInvalidFile () {
		// one row has too many values for the attribute.
		EqtlCovariate r = EqtlCovariate.parseFile(this.invalidCovars);
	}

	@Test
	public void testMissingDataFile () {
		EqtlCovariate r = EqtlCovariate.parseFile(this.missingCovars);
	}

	@Test
	public void testValidFile () {
		EqtlCovariate r = EqtlCovariate.parseFile(this.validCovars);
		Assert.assertNotNull(r);
	}

	/**
	 *  id      CHB12_P23_140809        CHB9_P24_140728 Genea42_P19_150107      HS346_P25_140704
	 *	foo     A       B       C       D
	 *	bar     E       F       G       H
	 */
	@Test
	public void testValidFileWithIgnored () {
		Set<String> ignoredDonors = new HashSet<>(Arrays.asList("Genea42_P19_150107", "CHB9_P24_140728"));
		List <String> expectedDonors = Arrays.asList("CHB12_P23_140809", "HS346_P25_140704");

		EqtlCovariate r = EqtlCovariate.parseFile(
				this.validCovars,
				ignoredDonors,
				ValidationStringency.SILENT
		);
		List<String> donors = r.donorNames();
		Assert.assertEquals(expectedDonors, donors);
	}

	@Test(expectedExceptions=TranscriptomeException.class)
	public void testNonNumericStrict() {
		EqtlCovariate.parseFile(this.nonnumericCovars,
				Collections.emptySet(),
				ValidationStringency.STRICT
		);
	}

	@Test
	public void testNonNumericLenient() {
		final EqtlCovariate covariate = EqtlCovariate.parseFile(
				this.nonnumericCovars,
				Collections.emptySet(),
				ValidationStringency.LENIENT
		);
		Assert.assertEquals(covariate.donorNames(), Arrays.asList("CHB12_P23_140809", "CHB9_P24_140728", "Genea42_P19_150107", "HS346_P25_140704"));
		Assert.assertEquals(covariate.getAttributes(), Arrays.asList("foo", "bar"));
		Assert.assertEquals(covariate.getValues("foo"), new String[]{"1", "B", "3", "4"});
		Assert.assertEquals(covariate.getValues("bar"), new String[]{"5.6", "7.8", "8.9", "H"});
	}

	@Test
	public void testNonNumericStrictIgnore() {
		Set<String> ignoredDonors = new HashSet<>(Arrays.asList("HS346_P25_140704", "CHB9_P24_140728"));
		final EqtlCovariate covariate = EqtlCovariate.parseFile(
				this.nonnumericCovars,
				ignoredDonors,
				ValidationStringency.STRICT
		);
		Assert.assertEquals(covariate.donorNames(), Arrays.asList("CHB12_P23_140809", "Genea42_P19_150107"));
		Assert.assertEquals(covariate.getAttributes(), Arrays.asList("foo", "bar"));
		Assert.assertEquals(covariate.getValues("foo"), new String[]{"1", "3"});
		Assert.assertEquals(covariate.getValues("bar"), new String[]{"5.6", "8.9"});
	}
}
