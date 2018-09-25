package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.junit.Assert;
import org.testng.annotations.Test;

public class BamTagOfTagCountsTest {

	private static final File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.bam");
	private static final File EXPECTED_OUT1 = new File("testdata/org/broadinstitute/transcriptome/barnyard/tag_of_tag_XC_XM.txt");
	private static final File EXPECTED_OUT2 = new File("testdata/org/broadinstitute/transcriptome/barnyard/tag_of_tag_XC_NM.txt");

	@Test
	public void testDoWorkXCXM() {
		File outFile=null;
		try {
			outFile = File.createTempFile("BamTagOfTagCounts.", ".out.txt");
		} catch (IOException e) {
			e.printStackTrace();
		}

		BamTagOfTagCounts b = new BamTagOfTagCounts();
		b.INPUT=IN_FILE;
		b.PRIMARY_TAG="XC";
		b.READ_QUALITY=0;
		b.SECONDARY_TAG="XM";
		b.OUTPUT=outFile;
		int r = b.doWork();
		Assert.assertTrue(r==0);

		try {
			Assert.assertTrue(FileUtils.contentEquals(outFile, EXPECTED_OUT1));
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	@Test
	public void testDoWorkXCNM() {
		File outFile=null;
		try {
			outFile = File.createTempFile("BamTagOfTagCounts.", ".out.txt");
		} catch (IOException e) {
			e.printStackTrace();
		}

		BamTagOfTagCounts b = new BamTagOfTagCounts();
		b.INPUT=IN_FILE;
		b.PRIMARY_TAG="XC";
		b.READ_QUALITY=0;
		b.SECONDARY_TAG="NM";
		b.OUTPUT=outFile;
		int r = b.doWork();
		Assert.assertTrue(r==0);

		try {
			Assert.assertTrue(FileUtils.contentEquals(outFile, EXPECTED_OUT2));
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
