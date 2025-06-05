package org.broadinstitute.dropseqrna.eqtl;

import java.io.File;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.dropseqrna.eqtl.CalculateXReactivationCovariate;
import org.broadinstitute.dropseqrna.eqtl.EqtlCovariate;
import org.testng.Assert;
import org.testng.annotations.Test;

public class CalculateXReactivationCovariateTest {
	private final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/eqtl/");
	private final File contigGroupFile=new File(TEST_DATA_DIR, "GRCh38_GRCm38.contig_groups.yaml");
	private final File gtfFile=new File(TEST_DATA_DIR, "GRCh38_GRCm38.smaller.gtf");
	private final File sequenceDictionaryFile = new File(TEST_DATA_DIR, "GRCh38_GRCm38.dict");
	private final File metaCellFile=new File(TEST_DATA_DIR, "d14-42_NGN2.meta_cell.expression.smaller.txt");
	private final File escapeGenesFile= new File(TEST_DATA_DIR, "x_escape_balaton_2015.txt.gz");
	private final File expected = new File(TEST_DATA_DIR, "d14-42_NGN2.maf_0.20_cisDist_10kb.expected_covars.txt");
	private final File expectedWithEscape = new File(TEST_DATA_DIR, "d14-42_NGN2.maf_0.20_cisDist_10kb.expected_covars_with_escape.txt");
	
	
	@Test
	// test without escape genes.
	public void getXReactivationCovariate() {
		CalculateXReactivationCovariate c = new CalculateXReactivationCovariate();
		EqtlCovariate result= c.getXReactivationCovariate(this.metaCellFile, this.gtfFile, this.contigGroupFile,
				this.sequenceDictionaryFile, null, ValidationStringency.STRICT);
		EqtlCovariate expected = EqtlCovariate.parseFile(this.expected);
		
		String [] actual = result.getValues(CalculateXReactivationCovariate.FRACTION_X_COVAR_NAME);		
		String [] exp = expected.getValues(CalculateXReactivationCovariate.FRACTION_X_COVAR_NAME);
		
		Assert.assertEquals(expected.donorNames(), result.donorNames());
		
		for (int i=0; i<actual.length; i++) {
			double a=Double.parseDouble(actual[i]);
			double e = Double.parseDouble(exp[i]);
			Assert.assertEquals(a, e, 0.0001);
		}						
	}
	
	@Test
	public void getXReactivationCovariateWithEscape() {
		CalculateXReactivationCovariate c = new CalculateXReactivationCovariate();
		EqtlCovariate result= c.getXReactivationCovariate(this.metaCellFile, this.gtfFile, this.contigGroupFile,
				this.sequenceDictionaryFile, escapeGenesFile, ValidationStringency.STRICT);
		EqtlCovariate expected = EqtlCovariate.parseFile(this.expectedWithEscape);
		
		String [] actual = result.getValues(CalculateXReactivationCovariate.FRACTION_X_COVAR_NAME);		
		String [] exp = expected.getValues(CalculateXReactivationCovariate.FRACTION_X_COVAR_NAME);
		
		Assert.assertEquals(expected.donorNames(), result.donorNames());
		
		for (int i=0; i<actual.length; i++) {
			double a=Double.parseDouble(actual[i]);
			double e = Double.parseDouble(exp[i]);
			Assert.assertEquals(a, e, 0.0001);
		}						
	}
}
