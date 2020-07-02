package org.broadinstitute.dropseqrna.vcftools.filters;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class AlleleFrequencyTagFilterTest {
	@Test
	public void testFilterOut() {
		Iterator<VariantContext> underlyingIterator = Collections.emptyIterator();
		String alleleFreqTag = "AF_TAG";

		VariantContext vc = getTestData(alleleFreqTag, 0.05);
		AlleleFrequencyTagFilter f = new AlleleFrequencyTagFilter(underlyingIterator, alleleFreqTag, 0.05);
		Assert.assertFalse(f.filterOut(vc));

		vc = getTestData(alleleFreqTag, 0.10);
		Assert.assertFalse(f.filterOut(vc));

		vc = getTestData(alleleFreqTag, 0.01);
		Assert.assertTrue(f.filterOut(vc));


	}

	private VariantContext getTestData (final String alleleFreqTag, final double alleleFreq) {
		VariantContextBuilder b = new VariantContextBuilder();
		VariantContext vc = b.alleles(Arrays.asList( Allele.create("A", true),  Allele.create("T", false))).start(1).stop(1).chr("1").attribute(alleleFreqTag, alleleFreq).make();
		return vc;
	}
}
