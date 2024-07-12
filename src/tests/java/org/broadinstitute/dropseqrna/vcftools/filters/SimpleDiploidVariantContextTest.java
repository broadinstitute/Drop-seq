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
package org.broadinstitute.dropseqrna.vcftools.filters;

import java.io.File;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class SimpleDiploidVariantContextTest {

	private final static String rootDir="testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/sampleassignment";
	
	private final File vcfFile = new File(rootDir+ "/test.vcf.gz");
	private final File symbolicAllelesVCFFile = new File(rootDir+ "/symbolic_alleles.vcf.gz");

	@Test
	public void testSymbolicAllelesFiltered() {
		VCFFileReader reader = new VCFFileReader(this.symbolicAllelesVCFFile, false);
		PeekableIterator<VariantContext> vcfIterator = new PeekableIterator<VariantContext>(reader.iterator());
		SimpleDiploidVariantContextFilter filter = new SimpleDiploidVariantContextFilter(vcfIterator);
		while (vcfIterator.hasNext()) {
			VariantContext data = vcfIterator.next();
			Assert.assertTrue(filter.filterOut(data));
		}
		reader.close();
	}

	@Test
	public void testFilter() {
		// in the test data, we expect the variants at these positions to be filtered.
		// add a SNP where the alt allele has 2 alt bases 1:243901
		// 1       243901  .       G       A,C
		// add a SNP where there are 2 alt alleles - 1:10439 (not a snp, an indel)
		// 1       10439   rs112766696     AC      A
		// add a SNP where the site is filtered - 1:12672  (VQSRTrancheSNP99.80to99.90 )

		VCFFileReader reader = new VCFFileReader(this.vcfFile, false);
		PeekableIterator<VariantContext> vcfIterator = new PeekableIterator<VariantContext>(reader.iterator());
		SimpleDiploidVariantContextFilter filter = new SimpleDiploidVariantContextFilter(vcfIterator);
		while (vcfIterator.hasNext())
			testIsFiltered(vcfIterator.next(), filter);
		reader.close();
	}


	private void testIsFiltered (final VariantContext data, final SimpleDiploidVariantContextFilter filter) {
		if (data.getStart()==10439)
			Assert.assertTrue(filter.filterOut(data));
		if (data.getStart()==12672)
			Assert.assertTrue(filter.filterOut(data));
		if (data.getStart()==243901)
			Assert.assertTrue(filter.filterOut(data));
		if (data.getStart()==76227022)
			Assert.assertFalse(filter.filterOut(data));
		if (data.getStart()==150199123)
			Assert.assertFalse(filter.filterOut(data));
	}

}
