package org.broadinstitute.dropseqrna.utils.readIterators;

import htsjdk.samtools.SAMRecord;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.junit.Assert;
import org.testng.annotations.Test;

public class SimpleTagIteratorTest {

	File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.bam");

	@Test
	public void testGetByBatch() {
		List<String>sortingTags = new ArrayList<String>();
		sortingTags.add("GE");
		sortingTags.add("XC");
		
		SimpleTagIterator iter = new SimpleTagIterator(IN_FILE, sortingTags, 10, true);
		
		while (iter.hasNext()) {
			Collection<SAMRecord> recs=iter.next();
			Assert.assertTrue(testAllRecordsSamTags(recs, sortingTags));
		}
		
	}
	
	private boolean testAllRecordsSamTags(Collection<SAMRecord> recs, List<String> sortingTags) {
		
		SAMRecord firstRec = recs.iterator().next();
		List<String> values = DEIteratorUtils.getValuesForTags(sortingTags, firstRec);
		System.out.println(buildLogString(values, sortingTags));
		for (SAMRecord r: recs) {
			List<String> newVals = DEIteratorUtils.getValuesForTags(sortingTags, firstRec);
			boolean notEqualFlag = DEIteratorUtils.testTagsNotEqual(values, newVals);
			if (notEqualFlag) return (false);
		}
		return true;		
	}
	
	private String buildLogString (List<String> values, List<String> sortingTags) {
		Assert.assertEquals(values.size(), sortingTags.size());
		
		StringBuilder b = new StringBuilder();
		for (int i=0; i<values.size(); i++) {
			b.append(sortingTags.get(i));
			b.append(" [" + values.get(i) +"] ");
		}
		return (b.toString());
	}
}
