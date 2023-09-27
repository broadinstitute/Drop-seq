package org.broadinstitute.dropseqrna.vcftools.filters;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.broadinstitute.dropseqrna.vcftools.filters.GenotypeGQFilter;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class GenotypeGQFilterTest {


	@Test
	public void testFilter() {
		VariantContextBuilder b = new VariantContextBuilder();
		Allele a1 = Allele.create("A", true);
		Allele a2 = Allele.create("T", false);

		final List<Allele> allelesRef = new ArrayList<>(Arrays.asList(a1, a1));

		// GQ=30 Genotypes
		Collection<Genotype> genotypes = new ArrayList<>();
		genotypes.add(new GenotypeBuilder("donor1", allelesRef).GQ(35).make());
		genotypes.add(new GenotypeBuilder("donor2", allelesRef).GQ(20).make());
		genotypes.add(new GenotypeBuilder("donor3", allelesRef).GQ(31).make());
		genotypes.add(new GenotypeBuilder("donor4", allelesRef).GQ(28).make());
		genotypes.add(new GenotypeBuilder("donor5", Arrays.asList(Allele.NO_CALL)).GQ(33).make());

		VariantContext vc = b.alleles(Arrays.asList(a1,a2)).start(1).stop(1).chr("1").genotypes(genotypes).make();
		GenotypesContext gc = vc.getGenotypes();

		GenotypeGQFilter f = new GenotypeGQFilter(gc.iterator(), 30);
		while (f.hasNext()) {
			Genotype g = f.next();
			String sampleName = g.getSampleName();
			Assert.assertTrue(sampleName.equals("donor1") || sampleName.equals("donor3"));
		}

	}
}
