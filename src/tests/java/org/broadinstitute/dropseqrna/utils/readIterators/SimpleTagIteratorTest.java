package org.broadinstitute.dropseqrna.utils.readIterators;

import htsjdk.samtools.SAMRecord;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.testng.annotations.Test;

public class SimpleTagIteratorTest {

	File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.bam");

	@Test
	public void testGetBatch() {
		List<String>sortingTags = new ArrayList<String>();
		sortingTags.add("GE");
		sortingTags.add("XC");
		
		SimpleTagIterator iter = new SimpleTagIterator(IN_FILE, sortingTags, 10, true);
		Collection<SAMRecord> recs=null;
		
		recs.size();
		
	}
}
