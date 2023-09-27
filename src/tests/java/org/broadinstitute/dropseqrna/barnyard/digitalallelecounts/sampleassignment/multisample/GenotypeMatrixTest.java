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
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.broadinstitute.dropseqrna.vcftools.filters.SimpleDiploidVariantContextFilter;
import org.testng.annotations.Test;

import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;


public class GenotypeMatrixTest {

	private final static String rootDir="testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/sampleassignment/";
	
	private final File VCF = new File(rootDir+"/test.vcf.gz");

	@Test()
	public void testMatrix1 () {

		final VCFFileReader vcfReader = new VCFFileReader(this.VCF, false);
		Set<String> samples = new HashSet<>(vcfReader.getFileHeader().getSampleNamesInOrder());
		Iterator<VariantContext> filteredIter=vcfReader.iterator();
		filteredIter = new SimpleDiploidVariantContextFilter(filteredIter, true, false, 2, true);

		GenotypeMatrix m = new GenotypeMatrix(filteredIter, 0, samples);


		// filtered out.
		Assert.assertNull(m.getGenotype(new Interval("1", 10439, 10439), "FOO"));
		Assert.assertNull(m.getGenotype(new Interval("1", 10439, 10439), "BAR"));

		Assert.assertEquals(GenotypeType.NO_CALL, m.getGenotype(new Interval("1", 12672, 12672), "FOO"));
		Assert.assertEquals(GenotypeType.NO_CALL, m.getGenotype(new Interval("1", 12672, 12672), "BAR"));

		// filtered out.
		Assert.assertNull(m.getGenotype(new Interval("1", 243901, 243901), "FOO"));
		Assert.assertNull(m.getGenotype(new Interval("1", 243901, 243901), "BAR"));

		// data!
		Assert.assertEquals(GenotypeType.HET, m.getGenotype(new Interval("1", 76227022, 76227022), "FOO"));
		Assert.assertEquals(GenotypeType.HET, m.getGenotype(new Interval("1", 76227022, 76227022), "BAR"));

		// data!
		Assert.assertEquals(GenotypeType.HOM_REF, m.getGenotype(new Interval("1", 150199123, 150199123), "FOO"));
		Assert.assertEquals(GenotypeType.HOM_REF, m.getGenotype(new Interval("1", 150199123, 150199123), "BAR"));

		vcfReader.close();
	}

	@Test()
	public void testMatrixGCFiltered () {

		final VCFFileReader vcfReader = new VCFFileReader(this.VCF, false);

		Iterator<VariantContext> filteredIter=vcfReader.iterator();
		filteredIter = new SimpleDiploidVariantContextFilter(filteredIter, true, false, 2, true);
		Set<String> samples = new HashSet<>(vcfReader.getFileHeader().getSampleNamesInOrder());
		GenotypeMatrix m = new GenotypeMatrix(filteredIter, 30, samples);


		// filtered out.
		Assert.assertNull(m.getGenotype(new Interval("1", 10439, 10439), "FOO"));
		Assert.assertNull(m.getGenotype(new Interval("1", 10439, 10439), "BAR"));

		Assert.assertEquals(GenotypeType.NO_CALL, m.getGenotype(new Interval("1", 12672, 12672), "FOO"));
		Assert.assertEquals(GenotypeType.NO_CALL, m.getGenotype(new Interval("1", 12672, 12672), "BAR"));

		// filtered out.
		Assert.assertNull(m.getGenotype(new Interval("1", 243901, 243901), "FOO"));
		Assert.assertNull(m.getGenotype(new Interval("1", 243901, 243901), "BAR"));

		// data!
		Assert.assertEquals(GenotypeType.NO_CALL, m.getGenotype(new Interval("1", 76227022, 76227022), "FOO"));
		Assert.assertEquals(GenotypeType.HET, m.getGenotype(new Interval("1", 76227022, 76227022), "BAR"));

		// data!
		Assert.assertEquals(GenotypeType.HOM_REF, m.getGenotype(new Interval("1", 150199123, 150199123), "FOO"));
		Assert.assertEquals(GenotypeType.HOM_REF, m.getGenotype(new Interval("1", 150199123, 150199123), "BAR"));

		vcfReader.close();
	}

	@Test
	public void testMatrix2 () {
		File vcf = new File (rootDir + "/multisample/TTTGCGCGGAGC:ATTGTTTAGGAG.vcf");
		final VCFFileReader vcfReader = new VCFFileReader(vcf, false);
		Set<String> samples = new HashSet<>(vcfReader.getFileHeader().getSampleNamesInOrder());
		Iterator<VariantContext> filteredIter=vcfReader.iterator();
		filteredIter = new SimpleDiploidVariantContextFilter(filteredIter, true, false, 2, true);

		GenotypeMatrix m = new GenotypeMatrix(filteredIter, 30, samples);
		Interval interval = new Interval ("1", 1337334, 1337334);
		GenotypeType t = m.getGenotype(interval, "HUES53");
		Assert.assertNotNull(t);
		double [] freqs = m.getGenotypeFrequencies(interval);
		Assert.assertNotNull(freqs);
		double [] expected = new double [] {0.6, 0.2, 0.2};
		for (int i=0; i<expected.length; i++)
			Assert.assertEquals(expected[i], freqs[i], 0.0001);

		vcfReader.close();
	}

}
