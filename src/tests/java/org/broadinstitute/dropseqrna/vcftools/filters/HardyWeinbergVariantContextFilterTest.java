package org.broadinstitute.dropseqrna.vcftools.filters;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class HardyWeinbergVariantContextFilterTest {

	private final File HWE_VCF = new File ("testdata/org/broadinstitute/dropseq/vcftools/filters/test_hwe.vcf.gz");

	@Test
	public void testHWE () {
		final VCFFileReader vcfReader = new VCFFileReader(this.HWE_VCF, false);
		for (VariantContext vc: vcfReader) {
			double hwe = HardyWeinbergVariantContextFilter.getHardyWeinbergPvalue(vc);
			checkHWEResult(vc, hwe);
		}
		vcfReader.close();
	}

	@Test
	public void testHWEFiltering () {
		final VCFFileReader vcfReader = new VCFFileReader(this.HWE_VCF, false);
		HardyWeinbergVariantContextFilter filter = new HardyWeinbergVariantContextFilter(vcfReader.iterator(), 0.001);
		// sites that have an MAF > threshold of 0.1
		Set<Integer> expectedSites = new HashSet<> (Arrays.asList(16050822,16050933,16051249,16051347,16051497,16051556));
		Set <Integer> observedSites = new HashSet<>();
		for (VariantContext vc: filter)
			observedSites.add(vc.getStart());
		Assert.assertEquals(observedSites, expectedSites);
		filter.close();
		vcfReader.close();
	}

	private void checkHWEResult (final VariantContext vc, final double hwe) {
		int start = vc.getStart();
		switch (start) {
		case 16050822: Assert.assertEquals(hwe, 0.0475866507,0.000001); break;
		case 16050933: Assert.assertEquals(hwe, 1.0,0.000001); break;
		case 16051249: Assert.assertEquals(hwe, 1.0,0.000001); break;
		case 16051347: Assert.assertEquals(hwe, 1.0,0.000001); break;
		case 16051497: Assert.assertEquals(hwe, 1.0,0.000001); break;
		case 16051556: Assert.assertEquals(hwe, 1.0,0.000001); break;
		case 16092075: Assert.assertEquals(hwe, 0.000298946,0.000001); break;
		case 16201235: Assert.assertEquals(hwe, 0.000554785,0.000001); break;
		case 16210537: Assert.assertEquals(hwe, 0.000554785,0.000001); break;
		case 16292829: Assert.assertEquals(hwe, 0.0006397894,0.000001); break;
		case 16963545: Assert.assertEquals(hwe, 0.00007108298,0.000001); break;
		case 17027476: Assert.assertEquals(hwe, 0.00001378143,0.000001); break;
		//
		default: break;
		}
	}
}
