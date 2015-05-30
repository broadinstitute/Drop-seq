package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.metrics.MetricsFile;

import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Assert;
import org.testng.annotations.Test;

import picard.analysis.RnaSeqMetrics;
import picard.analysis.directed.RnaSeqMetricsCollector;

public class SingleCellRnaSeqMetricsCollectorTest {
	
	/**
	 * Data for this test is generated from the following procedure (data is from the Barnyard_Runs2014/9-27-14_NextSeq/bams directory):
	 * /broad/mccarroll/software/dropseq/prod/BAMTagHistogram I=/broad/mccarroll/evan/Barnyard_Runs2015/4-19-15_NextSeq/mm10/bams/P3Bipolars.bam O=P3Bipolars_reads_histogram.txt TAG=XC READ_QUALITY=10
	 * head -n 10 P3Bipolars_reads_histogram.txt |cut -f 2 |tail -n 2 > 2mouse_cells.txt
	 * /broad/mccarroll/software/dropseq/jn_branch/FilterBAMByTag I=/broad/mccarroll/evan/Barnyard_Runs2015/4-19-15_NextSeq/mm10/bams/P3Bipolars.bam O=2MouseCells.bam TAG=XC TAG_VALUES_FILE=2mouse_cells.txt
	 * java -jar /seq/software/picard/current/bin/picard.jar DownsampleSam I=2MouseCells.bam O=2MouseCells.bam_downsampled.bam P=0.05
	 * /broad/mccarroll/software/dropseq/jn_branch/FilterBAMByTag I=2MouseCells.bam_downsampled.bam O=2MouseCells.bam_downsampled_CGTCACTTGCAC.bam TAG=XC TAG_VALUE=CGTCACTTGCAC
	 * /broad/mccarroll/software/dropseq/jn_branch/FilterBAMByTag I=2MouseCells.bam_downsampled.bam O=2MouseCells.bam_downsampled_TTCGCCCGGCTT.bam TAG=XC TAG_VALUE=TTCGCCCGGCTT
	 * java -jar /seq/software/picard/current/bin/picard.jar CollectRnaSeqMetrics I=2MouseCells.bam_downsampled_TTCGCCCGGCTT.bam O=TTCGCCCGGCTT.rna_seq_metrics.txt REF_FLAT=/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.refFlat RIBOSOMAL_INTERVALS=/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.rRNA.intervals STRAND_SPECIFICITY=NONE
	 * java -jar /seq/software/picard/current/bin/picard.jar CollectRnaSeqMetrics I=2MouseCells.bam_downsampled_CGTCACTTGCAC.bam O=CGTCACTTGCAC.rna_seq_metrics.txt REF_FLAT=/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.refFlat RIBOSOMAL_INTERVALS=/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.rRNA.intervals STRAND_SPECIFICITY=NONE
	 * /broad/mccarroll/software/dropseq/prod/SingleCellRnaSeqMetricsCollector I=2MouseCells.bam_downsampled.bam O=2MouseCells.bam_downsampled.rna_seq_metrics.txt CELL_BARCODE_TAG=XC ANNOTATIONS_FILE=/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.refFlat RIBOSOMAL_INTERVALS=/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.rRNA.intervals NUM_CORE_BARCODES=2 READ_MQ=0
	 * 
	 * Unfortunately this test relies on external data as well, so it can't be run as an automatic test so it's left disabled for normal testing.
	 */
  
	@Test(enabled=false)
	public void test1() {
		SingleCellRnaSeqMetricsCollector c = new SingleCellRnaSeqMetricsCollector();
		String cellBarcodeTag = "XC";
		Set<String> cellBarcodes = new HashSet<String>();
		cellBarcodes.add("TTCGCCCGGCTT");
		cellBarcodes.add("CGTCACTTGCAC");
		
		File inBAM=  new File("testdata/org/broadinstitute/transcriptome/barnyard/2MouseCells.bam_downsampled.bam");
		File annotationsFile = new File ("/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.refFlat");
		File rRNAIntervalsFile = new File ("/broad/mccarroll/software/metadata/individual_reference/mm10/mm10.rRNA.intervals");
		int readMQ=0;
		
		RnaSeqMetricsCollector collector = c.getRNASeqMetricsCollector(cellBarcodeTag, cellBarcodes, inBAM, 
	    		RnaSeqMetricsCollector.StrandSpecificity.NONE,0.8, readMQ, annotationsFile, rRNAIntervalsFile);
		
		final MetricsFile<RnaSeqMetrics, Integer> file = new MetricsFile<RnaSeqMetrics, Integer>();
		collector.addAllLevelsToFile(file);
		List<RnaSeqMetrics> metrics =  file.getMetrics();
		
		for (RnaSeqMetrics m: metrics) {
			if (m.SAMPLE.equals("TTCGCCCGGCTT")) {
				Assert.assertEquals(149785, m.PF_BASES);
				Assert.assertEquals(107319, m.PF_ALIGNED_BASES);
				Assert.assertEquals(new Long (762), m.RIBOSOMAL_BASES);
				Assert.assertEquals(53790L, m.CODING_BASES);
				Assert.assertEquals(26818L, m.UTR_BASES);
				Assert.assertEquals(9804L, m.INTRONIC_BASES);
				Assert.assertEquals(16146L, m.INTERGENIC_BASES);
			}
			if (m.SAMPLE.equals("CGTCACTTGCAC")) {
				Assert.assertEquals(163896, m.PF_BASES);
				Assert.assertEquals(120171, m.PF_ALIGNED_BASES);
				Assert.assertEquals(new Long (339), m.RIBOSOMAL_BASES);
				Assert.assertEquals(50169L, m.CODING_BASES);
				Assert.assertEquals(39045L, m.UTR_BASES);
				Assert.assertEquals(10622L, m.INTRONIC_BASES);
				Assert.assertEquals(19997L, m.INTERGENIC_BASES);
			}
		}
		
		
	}
}
