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
import htsjdk.samtools.metrics.MetricsFile;
import junit.framework.Assert;
import org.broadinstitute.dropseqrna.spermseq.metrics.duplicates.SpermSeqMarkDuplicates.DuplicateStrategy;
import org.broadinstitute.dropseqrna.spermseq.metrics.duplicates.SpermSeqMarkDuplicates.PCRDuplicateMetrics;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class SpermSeqMarkDuplicatesTest {


	File INPUT = new File ("testdata/org/broadinstitute/spermseq/metrics/duplicates/test_sorted.bam");

	@Test
	// tests which reads are marked as duplicates by the read position strategy.
	public void testDetectDuplicatesByReadPositionStrategy() throws IOException {
		String [] duplicateReadNames={"READ1:2", "READ2:3"};
		Set<String> dupes = new HashSet<String>(Arrays.asList(duplicateReadNames));

		SpermSeqMarkDuplicates d = new SpermSeqMarkDuplicates();
		d.INPUT=Arrays.asList(INPUT);
		d.OUTPUT=File.createTempFile("testDetectDuplicatesByReadPositionStrategy.", ".bam");
		d.OUTPUT.deleteOnExit();
		d.OUTPUT_STATS=File.createTempFile("testDetectDuplicatesByReadPositionStrategy.", ".pcr_duplicate_metrics");
		d.OUTPUT_STATS.deleteOnExit();
		Assert.assertEquals(0, d.doWork());

		SamReader inputSam = SamReaderFactory.makeDefault().open(d.OUTPUT);
		for (SAMRecord r: inputSam) {
			boolean duplicateReadFlag = r.getDuplicateReadFlag();
			String readName = r.getReadName();
			if (dupes.contains(readName))
				Assert.assertTrue(duplicateReadFlag);
			else
				Assert.assertFalse(duplicateReadFlag);
		}
		final List<SpermSeqMarkDuplicates.PCRDuplicateMetrics> beans = MetricsFile.readBeans(d.OUTPUT_STATS);
		Assert.assertEquals(1, beans.size());
		Assert.assertEquals(dupes.size(), beans.get(0).NUM_DUPLICATES);
	}

	@Test
	public void testMetrics () {
		SpermSeqMarkDuplicates d = new SpermSeqMarkDuplicates();
		d.INPUT=Arrays.asList(INPUT);
		d.CELL_BARCODE_TAG="XC";
		d.MOLECULAR_BARCODE_TAG="XM";
		d.NUM_BARCODES=10;
		Map<String, PCRDuplicateMetrics> metricsMap = d.getPerCellMetricsMap(d.getCellBarcodes ());

		// generate the data to update the metricsMap.
		SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(d.INPUT, false, SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE));
		GroupingIterator<SAMRecord> iter  = d.getCellGeneIterator(headerAndIterator);

		Assert.assertTrue(iter.hasNext());
		while (iter.hasNext()) {
			Collection<SAMRecord> batch = iter.next();
			Collection<SAMRecord> result = d.markDuplicates(batch, metricsMap);

		}
		Assert.assertNotNull(metricsMap);

		PCRDuplicateMetrics f1 = metricsMap.get("FAKE1");
		Assert.assertEquals(3, f1.NUM_READS);
		Assert.assertEquals(1, f1.NUM_DUPLICATES);

		PCRDuplicateMetrics f2 = metricsMap.get("FAKE2");
		Assert.assertEquals(3, f2.NUM_READS);
		Assert.assertEquals(1, f2.NUM_DUPLICATES);


	}



}
