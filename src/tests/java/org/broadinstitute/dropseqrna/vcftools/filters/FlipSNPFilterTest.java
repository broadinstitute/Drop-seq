package org.broadinstitute.dropseqrna.vcftools.filters;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class FlipSNPFilterTest {
	@Test
	public void testFilter() {

		List<List<String>> allBaseCombinations = getAllBaseCombinations();
		List<List<String>> filteredCombinations = getFilteredBaseCombinations();
		allBaseCombinations.removeAll(filteredCombinations);

		List<VariantContext> filtered = filteredCombinations.stream().map(x-> getVC(x)).collect(Collectors.toList());
		List<VariantContext> unfiltered = allBaseCombinations.stream().map(x-> getVC(x)).collect(Collectors.toList());

		Iterator<VariantContext> underlyingIterator = Collections.emptyIterator();
		FlipSNPFilter f = new FlipSNPFilter(underlyingIterator);
		for (VariantContext vc: filtered)
			Assert.assertTrue(f.filterOut(vc));
		for (VariantContext vc: unfiltered)
			Assert.assertFalse(f.filterOut(vc));

	}

	private VariantContext getVC (final List<String> bases) {
		VariantContextBuilder b = new VariantContextBuilder();
		return (b.alleles(Arrays.asList(Allele.create(bases.get(0), true), Allele.create(bases.get(1), false))).start(1).stop(1).chr("1").make());
	}

	private List<List<String>> getAllBaseCombinations () {
		List<String> bases = Arrays.asList("A", "C", "G", "T");
		List<List<String>> result = new ArrayList<>();
		for (String b1: bases)
			for (String b2: bases)
				if (!b1.equals(b2))
					result.add(Arrays.asList(b1, b2));

		return (result);
	}

	private List<List<String>> getFilteredBaseCombinations () {
		List<List<String>> result = new ArrayList<>();
		result.add(Arrays.asList("A", "T"));
		result.add(Arrays.asList("T", "A"));
		result.add(Arrays.asList("C", "G"));
		result.add(Arrays.asList("G", "C"));
		return (result);
	}
}
