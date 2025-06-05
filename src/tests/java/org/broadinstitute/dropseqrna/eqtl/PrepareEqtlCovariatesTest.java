package org.broadinstitute.dropseqrna.eqtl;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;

import org.broadinstitute.dropseqrna.eqtl.PrepareEqtlCovariates;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class PrepareEqtlCovariatesTest {

	private final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/eqtl/");
	private final File REF_COVARS1 = new File (TEST_DATA_DIR, "covars_test1.txt");
	private final File REF_COVARS2 = new File (TEST_DATA_DIR, "covars_test2.txt");
	private final File META_CELL=new File(TEST_DATA_DIR, "test1_2.meta_cell.expression.txt"); 
	private final File EQTL_COVARS_RESULT1=new File(TEST_DATA_DIR, "test1_2.expected_covars.txt");
	
	
	@Test
	public void testMultipleReferenceFiles() throws IOException {
		PrepareEqtlCovariates clp = new PrepareEqtlCovariates();
		clp.META_CELL_FILE=META_CELL;
		clp.COVARIATE_REFERENCE_FILE=Arrays.asList(REF_COVARS1, REF_COVARS2);
		clp.COVARIATE_FILES=Collections.emptyList();
		clp.OUTPUT=File.createTempFile("covars", ".txt");
		clp.OUTPUT.deleteOnExit();
		clp.REFERENCE_COVARIATE=Arrays.asList("SEX", "AGE", "NUM_CELLS", "AVERAGE_UMIS_PER_CELL");
		int ret = clp.doWork();
		Assert.assertEquals(ret, 0);
		Assert.assertTrue(TestUtils.testFilesSame(clp.OUTPUT, this.EQTL_COVARS_RESULT1));
		
	}
}
