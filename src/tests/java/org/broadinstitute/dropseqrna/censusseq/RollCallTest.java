package org.broadinstitute.dropseqrna.censusseq;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import org.testng.Assert;

public class RollCallTest {

	private static final List<File> IN_BAM = new ArrayList<File> (Collections.singletonList(new File("testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.bam")));
	private static final File IN_VCF = new File("testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.vcf.gz");
	
	private static final File OUT_ROLL_CALL = new File("testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.roll_call.txt");
	private static final File OUT_ROLL_CALL_VERBOSE = new File("testdata/org/broadinstitute/dropseq/censusseq/10_donors_chr22.selected_sites.roll_call.verbose.txt.gz");

	@BeforeClass
	public void beforeSuite() {
		TestUtils.setInflaterDeflaterIfMacOs();
	}
	
	
	@Test
	// Tests full path and result files.  Barebones but useful.
	// the math is checked more stringently in other unit tests.
	public void testRollCall() throws IOException {
		RollCall f = new RollCall();
		f.INPUT_BAM=IN_BAM;
		f.INPUT_VCF=IN_VCF;
		f.KNOWN_DONOR_TAG="ZS";
		f.MIN_BASE_QUALITY=null;
		f.OUTPUT=File.createTempFile("testRollCall.", ".roll_call.txt");
		f.OUTPUT.deleteOnExit();
		// f.USE_JDK_DEFLATER=true;		
		f.OUTPUT_VERBOSE=File.createTempFile("testRollCall.", ".roll_call.verbose.txt.gz");
		f.OUTPUT_VERBOSE.deleteOnExit();
		String TMP_DIR=f.OUTPUT.getParent();
		f.TMP_DIR=Arrays.asList(new File (TMP_DIR));
		int ret = f.doWork();
		Assert.assertTrue(ret==0);
		Assert.assertTrue(TestUtils.testFilesSame(OUT_ROLL_CALL, f.OUTPUT)); 
		Assert.assertTrue(TestUtils.testFilesSame(OUT_ROLL_CALL_VERBOSE, f.OUTPUT_VERBOSE));

	}

}
