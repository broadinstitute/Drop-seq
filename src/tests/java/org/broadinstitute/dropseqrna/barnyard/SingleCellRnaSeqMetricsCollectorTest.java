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
package org.broadinstitute.dropseqrna.barnyard;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.metrics.MetricsFile;
import picard.analysis.RnaSeqMetrics;
import picard.analysis.directed.RnaSeqMetricsCollector;

public class SingleCellRnaSeqMetricsCollectorTest {

	/**
	 * Data for this test is generated from the following procedure (data is from the Barnyard_Runs2014/9-27-14_NextSeq/bams directory):
	 * /broad/mccarroll/software/dropseq/prod/BamTagHistogram I=/broad/mccarroll/evan/Barnyard_Runs2015/4-19-15_NextSeq/mm10/bams/P3Bipolars.bam O=P3Bipolars_reads_histogram.txt TAG=XC READ_QUALITY=10
	 * head -n 10 P3Bipolars_reads_histogram.txt |cut -f 2 |tail -n 2 > 2mouse_cells.txt
	 * /broad/mccarroll/software/dropseq/jn_branch/FilterBamByTag I=/broad/mccarroll/evan/Barnyard_Runs2015/4-19-15_NextSeq/mm10/bams/P3Bipolars.bam O=2MouseCells.bam TAG=XC TAG_VALUES_FILE=2mouse_cells.txt
	 * java -jar /seq/software/picard/current/bin/picard.jar DownsampleSam I=2MouseCells.bam O=2MouseCells.bam_downsampled.bam P=0.05
	 * /broad/mccarroll/software/dropseq/jn_branch/FilterBamByTag I=2MouseCells.bam_downsampled.bam O=2MouseCells.bam_downsampled_CGTCACTTGCAC.bam TAG=XC TAG_VALUE=CGTCACTTGCAC
	 * /broad/mccarroll/software/dropseq/jn_branch/FilterBamByTag I=2MouseCells.bam_downsampled.bam O=2MouseCells.bam_downsampled_TTCGCCCGGCTT.bam TAG=XC TAG_VALUE=TTCGCCCGGCTT
	 * java -jar /seq/software/picard/current/bin/picard.jar CollectRnaSeqMetrics I=2MouseCells.bam_downsampled_TTCGCCCGGCTT.bam O=TTCGCCCGGCTT.rna_seq_metrics.txt REF_FLAT=/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.refFlat RIBOSOMAL_INTERVALS=/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.rRNA.intervals STRAND_SPECIFICITY=NONE
	 * java -jar /seq/software/picard/current/bin/picard.jar CollectRnaSeqMetrics I=2MouseCells.bam_downsampled_CGTCACTTGCAC.bam O=CGTCACTTGCAC.rna_seq_metrics.txt REF_FLAT=/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.refFlat RIBOSOMAL_INTERVALS=/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.rRNA.intervals STRAND_SPECIFICITY=NONE
	 * /broad/mccarroll/software/dropseq/prod/SingleCellRnaSeqMetricsCollector I=2MouseCells.bam_downsampled.bam O=2MouseCells.bam_downsampled.rna_seq_metrics.txt CELL_BARCODE_TAG=XC ANNOTATIONS_FILE=/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.refFlat RIBOSOMAL_INTERVALS=/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.rRNA.intervals NUM_CORE_BARCODES=2 READ_MQ=0
	 *
	 */

	private File IN_BAM=  new File("testdata/org/broadinstitute/transcriptome/barnyard/2MouseCells.bam_downsampled.bam");

	@Test(enabled=true)
	public void test1() {
		SingleCellRnaSeqMetricsCollector c = new SingleCellRnaSeqMetricsCollector();
		String cellBarcodeTag = "XC";
		List<String> cellBarcodes = new ArrayList<>();
		cellBarcodes.add("TTCGCCCGGCTT");
		cellBarcodes.add("CGTCACTTGCAC");

		File inBAM=  new File("testdata/org/broadinstitute/transcriptome/barnyard/2MouseCells.bam_downsampled.bam");
		File annotationsFile = new File ("testdata/org/broadinstitute/transcriptome/barnyard/mm10.refFlat.gz");
		File rRNAIntervalsFile = new File ("testdata/org/broadinstitute/transcriptome/barnyard/mm10.rRNA.intervals");
		int readMQ=0;

		// c.MT_SEQUENCE=Arrays.asList("MT");
		RnaSeqMetricsCollector collector = c.getRNASeqMetricsCollector(cellBarcodeTag, cellBarcodes, inBAM,
	    		RnaSeqMetricsCollector.StrandSpecificity.NONE,0.8, readMQ, annotationsFile, rRNAIntervalsFile);

		final MetricsFile<RnaSeqMetrics, Integer> file = new MetricsFile<>();
		collector.addAllLevelsToFile(file);
		List<RnaSeqMetrics> metrics =  file.getMetrics();

		for (RnaSeqMetrics m: metrics) {
			if (m.SAMPLE.equals("TTCGCCCGGCTT")) {
				Assert.assertEquals(149785, m.PF_BASES);
				Assert.assertEquals(107319, m.PF_ALIGNED_BASES);
				Assert.assertEquals(new Long (761), m.RIBOSOMAL_BASES);
				Assert.assertEquals(53790L, m.CODING_BASES);
				Assert.assertEquals(26818L, m.UTR_BASES);
				Assert.assertEquals(9804L, m.INTRONIC_BASES);
				Assert.assertEquals(16146L, m.INTERGENIC_BASES);
			}
			if (m.SAMPLE.equals("CGTCACTTGCAC")) {
				Assert.assertEquals(163896, m.PF_BASES);
				Assert.assertEquals(120171, m.PF_ALIGNED_BASES);
				Assert.assertEquals(new Long (338), m.RIBOSOMAL_BASES);
				Assert.assertEquals(50169L, m.CODING_BASES);
				Assert.assertEquals(39045L, m.UTR_BASES);
				Assert.assertEquals(10622L, m.INTRONIC_BASES);
				Assert.assertEquals(19997L, m.INTERGENIC_BASES);
			}
		}
	}

	@Test(enabled=true)
	public void testDoWork() {
		SingleCellRnaSeqMetricsCollector c = new SingleCellRnaSeqMetricsCollector();

		File annotationsFile = new File ("testdata/org/broadinstitute/transcriptome/barnyard/mm10.refFlat.gz");
		File rRNAIntervalsFile = new File ("testdata/org/broadinstitute/transcriptome/barnyard/mm10.rRNA.intervals");
		File cellBarcodeFile=new File ("testdata/org/broadinstitute/transcriptome/barnyard/SingleCellRnaSeqMetricsCollector.cellBarcodes.txt");
		File expectedOutFile=new File ("testdata/org/broadinstitute/transcriptome/barnyard/SingleCellRnaSeqMetricsCollector.expected_output.txt");

		File outFile=null;
		try {
			outFile = File.createTempFile("SingleCellRnaSeqMetricsCollector.", ".output.txt");
		} catch (IOException e) {
			e.printStackTrace();
		}

		outFile.deleteOnExit();

		c.INPUT=this.IN_BAM;
		c.MT_SEQUENCE=Arrays.asList("MT");
		c.OUTPUT=outFile;
		c.ANNOTATIONS_FILE=annotationsFile;
		c.RIBOSOMAL_INTERVALS=rRNAIntervalsFile;
		c.STRAND_SPECIFICITY = RnaSeqMetricsCollector.StrandSpecificity.NONE;
		c.CELL_BC_FILE=cellBarcodeFile;
		c.READ_MQ=0;
		int r = c.doWork();
		Assert.assertTrue(r==0);

		try {
			Assert.assertTrue (FileUtils.contentEquals(outFile, expectedOutFile));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}



	}

	@Test
	public void testGetCellBarcodes () {
		// test the alternate path for getting cell barcodes.
		SingleCellRnaSeqMetricsCollector c = new SingleCellRnaSeqMetricsCollector();
		List<String> cb = c.getCellBarcodes(null, this.IN_BAM, "XC", 0, 10);
		Assert.assertNotNull(cb);
		List<String> expected = Arrays.asList("CGTCACTTGCAC", "TTCGCCCGGCTT");
		Assert.assertEquals(expected, cb);

	}


}
