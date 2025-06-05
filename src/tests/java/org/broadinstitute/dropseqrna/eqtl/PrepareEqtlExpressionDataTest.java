/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.eqtl;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import java.util.Collections;

import org.broadinstitute.dropseqrna.eqtl.PrepareEqtlExpressionData;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class PrepareEqtlExpressionDataTest {

	private final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/eqtl/");
	private final File META_CELL=new File(TEST_DATA_DIR, "100_genes.48_donors.meta_cells.txt");
	private final File META_CELL_SMALLER=new File(TEST_DATA_DIR, "100_genes.38_donors.meta_cells.txt");
	private final File SD_FILE=new File(TEST_DATA_DIR, "hg19.dict");
	private final File ANNOTATIONS_FILE=new File(TEST_DATA_DIR, "100_genes.48_donors.gtf");
	private final File REJECTED_DONOR_LIST=new File(TEST_DATA_DIR, "100_genes.48_donors.rejected_donors.txt");
	
	private final File EXPRESSION_NO_NORM=new File(TEST_DATA_DIR, "100_genes.38_donors.expression.txt");
	private final File EXPRESSION_BED_FILE = new File(TEST_DATA_DIR, "100_genes.38_donors.expression.bed"
	);
	
	@Test 
	public void testExcludedDonors () throws IOException {
		PrepareEqtlExpressionData p = new PrepareEqtlExpressionData();
		p.ANNOTATIONS_FILE=this.ANNOTATIONS_FILE;		
		p.SD_FILE=SD_FILE;
		p.META_CELL_FILE=Arrays.asList(META_CELL_SMALLER);		
		p.REMOVE_PCT_EXPRESSION=0d;
		p.EXPRESSION_FILE=File.createTempFile("expression", ".txt");
		p.GENE_LOCATION_FILE=File.createTempFile("gene_location", ".txt");
		p.OUT_DONOR_LIST=File.createTempFile("donor_list", ".txt");
		p.EXPRESSION_FILE.deleteOnExit();
		p.GENE_LOCATION_FILE.deleteOnExit();
		p.OUT_DONOR_LIST.deleteOnExit();
		int ret = p.doWork();
		Assert.assertEquals(ret, 0);
		
		// The original meta cells file, with an argument to get rid of a subset of donors.
		// This should be equivalent to the first run.
		PrepareEqtlExpressionData p2 = new PrepareEqtlExpressionData();
		p2.ANNOTATIONS_FILE=this.ANNOTATIONS_FILE;		
		p2.SD_FILE=SD_FILE;
		p2.META_CELL_FILE=Arrays.asList(META_CELL);
		p2.REMOVE_PCT_EXPRESSION=0d;
		p2.EXPRESSION_FILE=File.createTempFile("expression_smaller", ".txt");
		p2.GENE_LOCATION_FILE=File.createTempFile("gene_location_smaller", ".txt");
		p2.OUT_DONOR_LIST=File.createTempFile("donor_list_smaller", ".txt");
		p2.REJECTED_DONOR_LIST=REJECTED_DONOR_LIST;
		p2.EXPRESSION_FILE.deleteOnExit();
		p2.GENE_LOCATION_FILE.deleteOnExit();
		p2.OUT_DONOR_LIST.deleteOnExit();
		int ret2 = p2.doWork();
		Assert.assertEquals(ret2, 0);
		
		Assert.assertTrue(TestUtils.testFilesSame(p.EXPRESSION_FILE, p2.EXPRESSION_FILE));
	}
	
	/*
	 * This test validates that when TRANSCRIPTS_PER_CELL = -1, the meta cell counts are passed through without normalization.
	 */
	@Test 
	public void testDisableNormalization () throws IOException {
		PrepareEqtlExpressionData p = new PrepareEqtlExpressionData();
		p.ANNOTATIONS_FILE=this.ANNOTATIONS_FILE;		
		p.SD_FILE=SD_FILE;
		p.META_CELL_FILE=Arrays.asList(META_CELL_SMALLER);		
		p.REMOVE_PCT_EXPRESSION=0d;
		p.EXPRESSION_FILE=File.createTempFile("expression", ".txt");
		p.GENE_LOCATION_FILE=File.createTempFile("gene_location", ".txt");
		p.OUT_DONOR_LIST=File.createTempFile("donor_list", ".txt");
		p.EXPRESSION_FILE.deleteOnExit();
		p.GENE_LOCATION_FILE.deleteOnExit();
		p.OUT_DONOR_LIST.deleteOnExit();
		p.TRANSCRIPTS_PER_CELL=-1;
		int ret = p.doWork();
		Assert.assertEquals(ret, 0);
				
		Assert.assertTrue(TestUtils.testFilesSame(p.EXPRESSION_FILE, EXPRESSION_NO_NORM));
	}

	@Test
	public void testGeneCounts() throws IOException {
		PrepareEqtlExpressionData p = new PrepareEqtlExpressionData();
		p.ANNOTATIONS_FILE = this.ANNOTATIONS_FILE;
		p.SD_FILE = SD_FILE;
		p.META_CELL_FILE = Collections.singletonList(META_CELL_SMALLER);
		p.REMOVE_PCT_EXPRESSION = 0d;
		p.TRANSCRIPTS_PER_CELL = -1;
		p.EXPRESSION_BED_FILE = File.createTempFile("expression", ".bed");
		p.OUT_DONOR_LIST = File.createTempFile("donor_list", ".txt");
		p.EXPRESSION_BED_FILE.deleteOnExit();
		p.OUT_DONOR_LIST.deleteOnExit();
		int ret = p.doWork();
		Assert.assertEquals(ret, 0);

		Assert.assertTrue(TestUtils.testFilesSame(p.EXPRESSION_BED_FILE, EXPRESSION_BED_FILE));
	}
}
