/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.broadinstitute.dropseqrna.spermseq.metrics.duplicates;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import junit.framework.Assert;
import org.broadinstitute.dropseqrna.spermseq.metrics.duplicates.SpermSeqMarkDuplicates.PCRDuplicateMetrics;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;
import org.broadinstitute.dropseqrna.utils.MultiComparator;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
import org.broadinstitute.dropseqrna.utils.readiterators.MissingTagFilteringIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class ReadClusterIteratorTest {

	File input = new File ("testdata/org/broadinstitute/spermseq/metrics/duplicates/TGATTAGGG_GAGGGGGGAGGGATAG_chr1.bam");
	private final String MOLECULAR_BARCODE_TAG="XM";
	private final String CELL_BARCODE_TAG="XC";
	private final String POS_TAG="ZZ";

	@Test
	public void test() {
		SpermSeqMarkDuplicates d = new SpermSeqMarkDuplicates();
		String cellBarcode="TGATTAGGG";

		Map<String, PCRDuplicateMetrics> metricsMap = new HashMap<String, PCRDuplicateMetrics>();
		PCRDuplicateMetrics m = d.new PCRDuplicateMetrics();
		metricsMap.put(cellBarcode, m);

		int windowSize=5000;
		GroupingIterator<SAMRecord> iter = getIterator();
		ReadClusterIterator rci = new ReadClusterIterator(iter, windowSize, CELL_BARCODE_TAG, MOLECULAR_BARCODE_TAG);
		while (rci.hasNext()) {
			ReadCluster rc = rci.next();
			Assert.assertNotNull(rc);
			testCluster(rc);
			Collection<SAMRecord> recs = d.markDuplicates(rc.getReads(), metricsMap);
			testClusterDuplicates(rc, recs);
		}
		rci.close();
	}


	private void testClusterDuplicates (final ReadCluster rc, final Collection<SAMRecord> recs) {
		// no recs of size 1 can have a duplicate.
		if (recs.size()==1)
			Assert.assertFalse(recs.iterator().next().getDuplicateReadFlag());

		if (rc.getClusterInterval().getStart()==608915)
			testDuplicates(recs, "HK5NTBGXX:2:11201:4839:16014");

		if (rc.getClusterInterval().getStart()==21158132)
			testDuplicates(recs, "HK5NTBGXX:1:12307:15437:10947");

		if (rc.getClusterInterval().getStart()==69503037)
			testDuplicates(recs, "HK5NTBGXX:3:12611:10046:8846");

		if (rc.getClusterInterval().getStart()==77465030)
			testDuplicates(recs, "HK5NTBGXX:3:12405:26846:8987");

		if (rc.getClusterInterval().getStart()==79479439)
			testDuplicates(recs, "HK5NTBGXX:3:21608:19950:4590");


	}

	/**
	 * Test a set of reads to see if the reads are dupe marked correctly.
	 * Only the notDuplicateReadName read will not be marked as a duplicate.
	 * @param recs
	 * @param notDuplicateReadName
	 */
	private void testDuplicates (final Collection<SAMRecord> recs, final String notDuplicateReadName) {
		for (SAMRecord r: recs) {
			String currentReadName = r.getReadName();
			if (currentReadName.equals(notDuplicateReadName))
				Assert.assertFalse(r.getDuplicateReadFlag());
			else
				Assert.assertTrue(r.getDuplicateReadFlag());
		}
	}

	private void testCluster (final ReadCluster rc) {
		if (rc.getClusterInterval().getStart()==608915)
			Assert.assertEquals(2, rc.getReads().size());
		else if (rc.getClusterInterval().getStart()==21158132)
			Assert.assertEquals(3, rc.getReads().size());
		else if (rc.getClusterInterval().getStart()==63628199)
			Assert.assertEquals(18, rc.getReads().size());
		else if (rc.getClusterInterval().getStart()==69503037)
			Assert.assertEquals(22, rc.getReads().size());
		else if (rc.getClusterInterval().getStart()==77465030)
			Assert.assertEquals(12, rc.getReads().size());
		else if (rc.getClusterInterval().getStart()==79479439)
			Assert.assertEquals(5, rc.getReads().size());
		else
			Assert.assertEquals(1, rc.getReads().size());


	}

	private GroupingIterator<SAMRecord> getIterator() {
		SamReader inputSam = SamReaderFactory.makeDefault().open(input);
		final Iterator<SAMRecord> filteringIterator = new MissingTagFilteringIterator(inputSam.iterator(), this.MOLECULAR_BARCODE_TAG, this.CELL_BARCODE_TAG);
	    // sort by Cell and molecular barcode and position.
	    @SuppressWarnings("unchecked")
		final MultiComparator<SAMRecord> comparator = new MultiComparator<SAMRecord>(
	            new StringTagComparator(this.CELL_BARCODE_TAG), new StringTagComparator(this.MOLECULAR_BARCODE_TAG), new IntervalTagComparator(this.POS_TAG, inputSam.getFileHeader().getSequenceDictionary()));
	    // add the position tag.
	    final ReadDuplicateWrapper sortingIteratorWrapper = new ReadDuplicateWrapper(filteringIterator, POS_TAG);
	    final CloseableIterator<SAMRecord> sortingIterator =
	            SamRecordSortingIteratorFactory.create(inputSam.getFileHeader(), sortingIteratorWrapper, comparator, null);
	    @SuppressWarnings("unchecked")
	    final MultiComparator<SAMRecord> groupingComparator = new MultiComparator<SAMRecord>(
	            new StringTagComparator(this.CELL_BARCODE_TAG), new StringTagComparator(this.MOLECULAR_BARCODE_TAG));

	    final GroupingIterator<SAMRecord> groupingIterator = new GroupingIterator<>(sortingIterator, groupingComparator);
	    return groupingIterator;
	}

}
