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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import picard.annotation.LocusFunction;

public class MultiCellDigitalAlleleCountsTest {

	private final File smallBAMFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/smallTest_snpUMIPileUp_retagged.sam");

	private final File bigBAMFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_5_cell_2_snp_testdata_retagged.bam");

	private final File cellBCFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_cells_cell_barcodes.txt");
	private final File snpIntervalsFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_cells_2snps.intervals");

	private final File cellClusterFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/clusters.txt");

	private final String cellBarcodeTag = "XC";
	private final String molBCTag = "XM";
	private final String snpTag = "YS";
	private final String functionTag="XF";
	private final int readMQ = 10;
	private final boolean assignReadsToAllGenes = true;

	private String GENE_NAME_TAG="gn";
	private String GENE_STRAND_TAG="gs";
	private String GENE_FUNCTION_TAG="gf";
	private StrandStrategy STRAND_STRATEGY = StrandStrategy.SENSE;
	private List<LocusFunction> LOCUS_FUNCTION_LIST=new ArrayList<LocusFunction>(Arrays.asList(LocusFunction.CODING, LocusFunction.UTR));

	@Test
	public void getMetaAnalysis() {
		List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
		IntervalList snpIntervals = IntervalList.fromFile(snpIntervalsFile);
		int baseQualityThreshold=10;

		SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
				new SamHeaderAndIterator(smallBAMFile), snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
				LOCUS_FUNCTION_LIST, STRAND_STRATEGY, GeneFunctionCommandLineBase.DEFAULT_FUNCTIONAL_STRATEGY, cellBarcodeTag,
				molBCTag, snpTag, functionTag, readMQ, assignReadsToAllGenes, cellBarcodes, null, SortOrder.SNP_GENE);

		MultiCellDigitalAlleleCountsIterator multiIter = new MultiCellDigitalAlleleCountsIterator(new DigitalAlleleCountsIterator(sbpi, baseQualityThreshold));
		MultiCellDigitalAlleleCounts r = multiIter.next();

		DigitalAlleleCounts c = r.getMetaAnalysis();
		c.toString();
		Assert.assertNotNull(c);
	}


	@Test
	public void getMetaAnalysisSpecificCells() {

		Map<String, Set<String>> clusterMap = ParseBarcodeFile.readCellClusterFile(this.cellClusterFile);

		List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
		IntervalList snpIntervals = IntervalList.fromFile(snpIntervalsFile);
		int baseQualityThreshold=10;

		SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
				new SamHeaderAndIterator(bigBAMFile), snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
				LOCUS_FUNCTION_LIST, STRAND_STRATEGY, GeneFunctionCommandLineBase.DEFAULT_FUNCTIONAL_STRATEGY, cellBarcodeTag,
				molBCTag, snpTag, functionTag, readMQ, assignReadsToAllGenes, cellBarcodes, null, SortOrder.SNP_GENE);

		MultiCellDigitalAlleleCountsIterator multiIter = new MultiCellDigitalAlleleCountsIterator(new DigitalAlleleCountsIterator(sbpi, baseQualityThreshold));
		MultiCellDigitalAlleleCounts r = multiIter.next();
		for (String cell: r.getCells()) {
			DigitalAlleleCounts one = r.getDigitalAlleleCounts(cell);
			// System.out.println(one);
			one.getCell();
		}

		List<String> clusterIDs=new ArrayList<> (clusterMap.keySet());
		Collections.sort(clusterIDs);

		for (String clusterID: clusterIDs) {
			Set<String> barcodes = clusterMap.get(clusterID);
			DigitalAlleleCounts dac = r.getMetaAnalysis(clusterID, barcodes);
			if (clusterID.equals("1")) {
				ObjectCounter<Character> o = dac.getUMIAlleleCount();
				Assert.assertEquals(o.getCountForKey('A'), 37);
				Assert.assertEquals(o.getCountForKey('G'), 31);
				ObjectCounter<Character> z = dac.getReadCounts();
				Assert.assertEquals(z.getCountForKey('A'), 189);
				Assert.assertEquals(z.getCountForKey('G'), 175);

			}
			if (clusterID.equals("2")) {
				ObjectCounter<Character> o = dac.getUMIAlleleCount();
				Assert.assertEquals(o.getCountForKey('A'), 10);
				Assert.assertEquals(o.getCountForKey('G'), 5);
				ObjectCounter<Character> z = dac.getReadCounts();
				Assert.assertEquals(z.getCountForKey('A'), 47);
				Assert.assertEquals(z.getCountForKey('G'), 7);
			}
		}
		DigitalAlleleCounts c = r.getMetaAnalysis();
		c.toString();
		Assert.assertNotNull(c);
		multiIter.close();
	}

	@Test
	public void testUMIPurityFilter () {
		Interval snpInterval = new Interval("1", 1, 1);

		DigitalAlleleCounts dac = new DigitalAlleleCounts(snpInterval, "FOO", "CELL", 10, 'N', 'N');
		ObjectCounter<Character> umi1 = new ObjectCounter<>();
		umi1.incrementByCount('T', 4);
		umi1.incrementByCount('A', 1);
		dac.addReadsForUMI("FILTERED_UMI", umi1);
		Assert.assertEquals(dac.getUMIPurity("FILTERED_UMI"), 0.8, 0.001);


		ObjectCounter<Character> umi2 = new ObjectCounter<>();
		umi2.incrementByCount('T', 4);
		dac.addReadsForUMI("RETAINED_UMI", umi2);
		Assert.assertEquals(dac.getUMIPurity("RETAINED_UMI"), 1, 0.001);

		MultiCellDigitalAlleleCounts d = new MultiCellDigitalAlleleCounts(snpInterval, "FOO");
		d.add(dac);
		DigitalAlleleCounts dac1 = d.getDigitalAlleleCounts("CELL");

		double purity = dac1.getMeanUMIPurity();
		Assert.assertEquals(purity, 0.9, 0.01);

		d.filterDataByUMIPurity(1);

		double purity2 = d.getDigitalAlleleCounts("CELL").getMeanUMIPurity();
		Assert.assertEquals(purity2, 1, 0.01);

	}


}
