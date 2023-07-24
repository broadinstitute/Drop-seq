package org.broadinstitute.dropseqrna.utils;

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

public class VariantContextSingletonFilterTest {

  @Test
  public void filterOutHetSinglton() {
	  VariantContextBuilder b = new VariantContextBuilder();
	  Allele a1 = Allele.create("A", true);
	  Allele a2 = Allele.create("T", false);

	  final List<Allele> allelesHet = new ArrayList<>(Arrays.asList(a1,a2));
	  final List<Allele> allelesRef = new ArrayList<>(Arrays.asList(a1,a1));
	  
	  // this has a het singleton.
	  Collection<Genotype> genotypes = new ArrayList<>();
	  genotypes.add(GenotypeBuilder.create("donor1", allelesRef));
	  genotypes.add(GenotypeBuilder.create("donor2", allelesHet));
	  genotypes.add(GenotypeBuilder.create("donor3", allelesRef));

	  VariantContext vc = b.alleles(allelesHet).start(1).stop(1).chr("1").genotypes(genotypes).make();
	  Assert.assertNotNull(vc);

	  Iterator<VariantContext> underlyingIterator = Collections.emptyIterator();

	  VariantContextSingletonFilter f = new VariantContextSingletonFilter(underlyingIterator, true);
	  boolean t1 = f.filterOut(vc);
	  Assert.assertFalse(t1);
	  f.close();

  }

  @Test
  public void filterOutHetSinglton2() {
	  VariantContextBuilder b = new VariantContextBuilder();
	  Allele a1 = Allele.create("A", true);
	  Allele a2 = Allele.create("T", false);

	  final List<Allele> allelesHet = new ArrayList<>(Arrays.asList(a1,a2));
	  final List<Allele> allelesRef = new ArrayList<>(Arrays.asList(a1,a1));
	  

	  // this does not have a het singleton.
	  Collection<Genotype> genotypes = new ArrayList<>();
	  genotypes.add(GenotypeBuilder.create("donor1", allelesRef));
	  genotypes.add(GenotypeBuilder.create("donor2", allelesHet));
	  genotypes.add(GenotypeBuilder.create("donor3", allelesHet));

	  VariantContext vc = b.alleles(allelesHet).start(1).stop(1).chr("1").genotypes(genotypes).make();
	  Assert.assertNotNull(vc);

	  Iterator<VariantContext> underlyingIterator = Collections.emptyIterator();

	  VariantContextSingletonFilter f = new VariantContextSingletonFilter(underlyingIterator, true);
	  boolean t1 = f.filterOut(vc);
	  Assert.assertTrue(t1);
	  f.close();
	  
  }

  @Test
  public void filterOutHetSinglto32() {
	  VariantContextBuilder b = new VariantContextBuilder();
	  Allele a1 = Allele.create("A", true);
	  Allele a2 = Allele.create("T", false);

	  final List<Allele> allelesHet = new ArrayList<>(Arrays.asList(a1,a2));
	  final List<Allele> allelesRef = new ArrayList<>(Arrays.asList(a1,a1));
	  final List<Allele> allelesAlt = new ArrayList<>(Arrays.asList(a2,a2));

	  // this does not have a het singleton because it has an alt.
	  Collection<Genotype> genotypes = new ArrayList<>();
	  genotypes.add(GenotypeBuilder.create("donor1", allelesRef));
	  genotypes.add(GenotypeBuilder.create("donor2", allelesRef));
	  genotypes.add(GenotypeBuilder.create("donor3", allelesAlt));

	  VariantContext vc = b.alleles(allelesHet).start(1).stop(1).chr("1").genotypes(genotypes).make();
	  Assert.assertNotNull(vc);

	  Iterator<VariantContext> underlyingIterator = Collections.emptyIterator();

	  VariantContextSingletonFilter f = new VariantContextSingletonFilter(underlyingIterator, true);
	  boolean t1 = f.filterOut(vc);
	  Assert.assertTrue(t1);
	  f.close();

  }

}
