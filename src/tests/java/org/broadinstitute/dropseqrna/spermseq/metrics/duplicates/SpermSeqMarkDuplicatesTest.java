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
import junit.framework.Assert;
import org.broadinstitute.dropseqrna.spermseq.metrics.duplicates.SpermSeqMarkDuplicates.DuplicateStrategy;
import org.broadinstitute.dropseqrna.spermseq.metrics.duplicates.SpermSeqMarkDuplicates.PCRDuplicateMetrics;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class SpermSeqMarkDuplicatesTest {


	File input = new File ("testdata/org/broadinstitute/spermseq/metrics/duplicates/test_sorted.bam");
	File output = new File ("testdata/org/broadinstitute/spermseq/metrics/duplicates/test_output.bam");
	File outputStats = new File ("testdata/org/broadinstitute/spermseq/metrics/duplicates/test_output.stats");

	private SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE);

	@Test(enabled=true)
	// tests which reads are marked as duplicates by the read position strategy.
	public void testDetectDuplicatesByReadPositionStrategy() {
		String [] duplicateReadNames={"READ1:2", "READ2:3"};
		Set<String> dupes = new HashSet<String>(Arrays.asList(duplicateReadNames));

		SpermSeqMarkDuplicates d = new SpermSeqMarkDuplicates();
		List<File> INPUT = new ArrayList<File> ();
		INPUT.add(this.input);
		d.INPUT=INPUT;
		d.STRATEGY=DuplicateStrategy.READ_POSITION;
		d.OUTPUT=this.output;
		d.OUTPUT_STATS=this.outputStats;
		d.doWork();

		SamReader inputSam = SamReaderFactory.makeDefault().open(output);
		for (SAMRecord r: inputSam) {
			boolean duplicateReadFlag = r.getDuplicateReadFlag();
			String readName = r.getReadName();
			if (dupes.contains(readName))
				Assert.assertTrue(duplicateReadFlag);
			else
				Assert.assertFalse(duplicateReadFlag);
		}

		// cleanup.
		output.delete();
		outputStats.delete();
	}

	@Test(enabled=false)
	// tests which reads are marked as duplicates by the read position strategy.
	public void testDetectDuplicatesByClusterStrategy() {
		String [] duplicateReadNames={"READ1:2", "READ2:3"};
		Set<String> dupes = new HashSet<String>(Arrays.asList(duplicateReadNames));

		SpermSeqMarkDuplicates d = new SpermSeqMarkDuplicates();
		List<File> INPUT = new ArrayList<File> ();
		INPUT.add(this.input);
		d.INPUT=INPUT;
		d.STRATEGY=DuplicateStrategy.CLUSTER;
		d.OUTPUT=this.output;
		d.OUTPUT_STATS=this.outputStats;
		d.doWork();

		SamReader inputSam = SamReaderFactory.makeDefault().open(output);
		for (SAMRecord r: inputSam) {
			boolean duplicateReadFlag = r.getDuplicateReadFlag();
			String readName = r.getReadName();
			if (dupes.contains(readName))
				Assert.assertTrue(duplicateReadFlag);
			else
				Assert.assertFalse(duplicateReadFlag);
		}

		// cleanup.
		output.delete();
		outputStats.delete();
	}


	@Test
	public void testMetrics () {
		List<File> INPUT = new ArrayList<File> ();
		INPUT.add(this.input);

		SpermSeqMarkDuplicates d = new SpermSeqMarkDuplicates();
		d.INPUT=INPUT;
		d.CELL_BARCODE_TAG="XC";
		d.MOLECULAR_BARCODE_TAG="XM";
		d.NUM_BARCODES=10;
		Map<String, PCRDuplicateMetrics> metricsMap = d.getPerCellMetricsMap(d.getCellBarcodes ());

		// generate the data to update the metricsMap.
		SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(INPUT, false, samReaderFactory);
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
