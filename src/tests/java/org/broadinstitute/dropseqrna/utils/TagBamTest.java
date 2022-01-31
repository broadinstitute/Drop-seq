package org.broadinstitute.dropseqrna.utils;

import java.io.File;
import java.io.IOException;

import org.testng.annotations.Test;

import junit.framework.Assert;

public class TagBamTest {

	
	File INPUT = new File ("testdata/org/broadinstitute/dropseq/utils/unpaired_reads_tagged.bam");
	File EXPECTED = new File  ("testdata/org/broadinstitute/dropseq/utils/unpaired_reads_tagged_tag_result.bam");
	
	@Test
	public void doWorkTest() throws IOException {
		TagBam t = new TagBam();
		t.INPUT=INPUT;
		t.OUTPUT=File.createTempFile("TagBam", "test.bam");
		t.TAG_NAME="ZT";
		t.TAG_VALUE="TEST";
		t.OUTPUT.deleteOnExit();
		int result = t.doWork();
		Assert.assertEquals(0, result);
		
		TestUtils.assertSamRecordsSame(t.OUTPUT, EXPECTED);
		
		
	}
}
