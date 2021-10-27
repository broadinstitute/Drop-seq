package org.broadinstitute.dropseqrna.vcftools.filters;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class MinorAlleleFreqVariantContextFilterTest {

	private final static String rootDir="testdata/org/broadinstitute/dropseq/vcftools/filters";
	
	private final File MAF_VCF = new File (rootDir+ "/test_maf.vcf.gz");

	@Test
	public void testMAF () {
		final VCFFileReader vcfReader = new VCFFileReader(this.MAF_VCF, false);
		for (VariantContext vc: vcfReader) {
			double maf = MinorAlleleFreqVariantContextFilter.calculateMinorAlleleFrequency(vc, 30);
			checkMAFResult(vc, maf);
		}
		vcfReader.close();
	}

	@Test
	public void testMAFFiltering () {
		final VCFFileReader vcfReader = new VCFFileReader(this.MAF_VCF, false);
		MinorAlleleFreqVariantContextFilter filter = new MinorAlleleFreqVariantContextFilter(vcfReader.iterator(), 0.1, true);
		// sites that have an MAF > threshold of 0.1
		Set<Integer> expectedSites = new HashSet<> (Arrays.asList(16050822,16051347));
		Set <Integer> observedSites = new HashSet<>();
		for (VariantContext vc: filter)
			observedSites.add(vc.getStart());
		Assert.assertTrue(expectedSites.containsAll(observedSites));
		Assert.assertTrue(observedSites.containsAll(expectedSites));
		filter.close();
		vcfReader.close();
	}

	private void checkMAFResult (final VariantContext vc, final double maf) {
		int start = vc.getStart();
		switch (start) {
		case 16050822: Assert.assertEquals(maf, 0.235849,0.0001); break;
		case 16050933: Assert.assertEquals(maf, 0.0849057,0.0001); break;
		case 16051249: Assert.assertEquals(maf, 0.0943396,0.0001); break;
		case 16051347: Assert.assertEquals(maf, 0.301887,0.0001); break;
		case 16055207: Assert.assertEquals(maf, 0.00943396,0.0001); break;
		case 16055484: Assert.assertEquals(maf, 0.0,0.0001); break;
		default: break;
		}
	}
}
