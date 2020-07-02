package org.broadinstitute.dropseqrna.censusseq;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.broadinstitute.dropseqrna.vcftools.filters.MonomorphicVariantContextFilter;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class MonomorphicVariantContextFilterTest {
	@Test
	public void testFilter() {

		VariantContextBuilder b = new VariantContextBuilder();
		Allele a1 = Allele.create("A", true);
		Allele a2 = Allele.create("T", false);

		final List<Allele> allelesHet = new ArrayList<>(Arrays.asList(a1, a2));
		final List<Allele> allelesRef = new ArrayList<>(Arrays.asList(a1, a1));
		final List<Allele> allelesAlt = new ArrayList<>(Arrays.asList(a2, a2));

		// this has a het singleton.
		Collection<Genotype> genotypes = new ArrayList<>();
		genotypes.add(GenotypeBuilder.create("donor1", allelesRef));
		genotypes.add(GenotypeBuilder.create("donor2", allelesHet));
		genotypes.add(GenotypeBuilder.create("donor3", allelesRef));
		genotypes.add(GenotypeBuilder.create("donor4", allelesAlt));

		VariantContext vc = b.alleles(allelesHet).start(1).stop(1).chr("1").genotypes(genotypes).make();
		Assert.assertNotNull(vc);

		Iterator<VariantContext> underlyingIterator = Collections.emptyIterator();

		// depending on which donors are included in the list, this VC is retained or
		// filtered.

		// filter when 2 refs.
		List<String> vcfSamples = Arrays.asList("donor1", "donor3");
		MonomorphicVariantContextFilter f = new MonomorphicVariantContextFilter(underlyingIterator, vcfSamples);
		boolean t1 = f.filterOut(vc);
		Assert.assertTrue(t1);

		// don't  filter when 2 refs+het.
		vcfSamples = Arrays.asList("donor1", "donor2", "donor3");
		f = new MonomorphicVariantContextFilter(underlyingIterator, vcfSamples);
		t1 = f.filterOut(vc);
		Assert.assertFalse(t1);

		// don't filter when 2 refs+alt.
		vcfSamples = Arrays.asList("donor1", "donor3", "donor4");
		f = new MonomorphicVariantContextFilter(underlyingIterator, vcfSamples);
		t1 = f.filterOut(vc);
		Assert.assertFalse(t1);
	}
}
