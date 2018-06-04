package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class TagValueFilteringIteratorTest {

	private static final File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.bam");

	@Test
	public void testCellBarcodeFiltering () {
		String cellBarcodeTag = "XC";

		SamReader inputSam = SamReaderFactory.makeDefault().open(IN_FILE);
		String [] desiredBarcodes = {"ATCAGGGACAGA", "TTGCCTTACGCG", "TACAATTAAGGC"};

		TagValueFilteringIterator<String> iter = new TagValueFilteringIterator<String>(inputSam.iterator(), cellBarcodeTag, Arrays.asList(desiredBarcodes));
		Set<String> barcodesFound = new HashSet<String>();
		while (iter.hasNext()) {
			SAMRecord r = iter.next();
			String v = r.getStringAttribute(cellBarcodeTag);
			barcodesFound.add(v);
		}

		// should be the same size.
		Assert.assertEquals(barcodesFound.size(), desiredBarcodes.length);
		// should contain the same answers.
		for (String s: desiredBarcodes)
			Assert.assertTrue(barcodesFound.contains(s));

	}
}
