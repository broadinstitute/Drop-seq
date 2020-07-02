package org.broadinstitute.dropseqrna.vcftools.filters;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class CallRateVariantContextFilterTest {
	@Test
	// test all donors.
	public void testFilter1() {
		VariantContextBuilder b = new VariantContextBuilder();
		Allele a1 = Allele.create("A", true);
		Allele a2 = Allele.create("T", false);

		final List<Allele> allelesHet = new ArrayList<>(Arrays.asList(a1, a2));
		final List<Allele> allelesRef = new ArrayList<>(Arrays.asList(a1, a1));
		final List<Allele> allelesAlt = new ArrayList<>(Arrays.asList(a2, a2));

		// GQ=30 Genotypes
		Collection<Genotype> genotypes = new ArrayList<>();
		genotypes.add(new GenotypeBuilder("donor1", allelesRef).GQ(30).make());
		genotypes.add(new GenotypeBuilder("donor2", allelesHet).GQ(30).make());
		genotypes.add(new GenotypeBuilder("donor3", allelesRef).GQ(30).make());
		genotypes.add(new GenotypeBuilder("donor4", allelesAlt).GQ(30).make());

		VariantContext vc = b.alleles(allelesHet).start(1).stop(1).chr("1").genotypes(genotypes).make();
		Assert.assertNotNull(vc);
		Iterator<VariantContext> underlyingIterator = Collections.emptyIterator();

		// test against 4 well called donors.  Don't filter.
		CallRateVariantContextFilter f = new CallRateVariantContextFilter(underlyingIterator, 30, 0.9);
		Assert.assertFalse(f.filterOut(vc));

		// test against 4 well called donors + 1 donor with low GQ.  Filter.
		genotypes.add(new GenotypeBuilder("donor5", allelesAlt).GQ(10).make());
		vc = b.alleles(allelesHet).start(1).stop(1).chr("1").genotypes(genotypes).make();
		f = new CallRateVariantContextFilter(underlyingIterator, 30, 0.9);
		Assert.assertTrue(f.filterOut(vc));

		// test against 4 well called donors, skip the 5th because it's not on the donor list.
		List<String> donors = Arrays.asList("donor1", "donor2", "donor3", "donor4");
		f = new CallRateVariantContextFilter(underlyingIterator, 30, 0.9, donors);
		Assert.assertFalse(f.filterOut(vc));

		// add the 5th donor back and filter.
		donors = Arrays.asList("donor1", "donor2", "donor3", "donor4", "donor5");
		f = new CallRateVariantContextFilter(underlyingIterator, 30, 0.9, donors, true);
		Assert.assertTrue(f.filterOut(vc));

		//construct a no call genotype.  This also fails and gets filtered.
		genotypes.add(new GenotypeBuilder("donor6", Arrays.asList(Allele.NO_CALL)).make());
		vc = b.alleles(allelesHet).start(1).stop(1).chr("1").genotypes(genotypes).make();
		donors = Arrays.asList("donor1", "donor2", "donor3", "donor4", "donor6");
		f = new CallRateVariantContextFilter(underlyingIterator, 30, 0.9, donors, true);
		Assert.assertTrue(f.filterOut(vc));


	}
}
