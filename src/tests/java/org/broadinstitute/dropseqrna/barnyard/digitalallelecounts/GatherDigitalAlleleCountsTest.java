package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.GatherDigitalAlleleCounts;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.annotations.Test;

import htsjdk.samtools.util.TestUtil;
import org.testng.Assert;

public class GatherDigitalAlleleCountsTest {

	private final static String rootDir="testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/";
	
	private final File largeBAMFile = new File(rootDir, "hek_5_cell_2_snp_testdata_retagged.bam");
	private final File cellBCFile = new File(rootDir, "hek_cells_cell_barcodes.txt");	
	private final File clusterFile = new File(rootDir, "hek_cluster_cell_barcodes.txt");	
			
	private final File vcfFile = new File(rootDir,"hek_cells_2snps.vcf");
	private final File sampleFile = new File(rootDir,"hek_cells_2snps.sample_list");
	
	private final File EXPECTED_OUTFILE = new File (rootDir, "hek_5_cell_2_snp_testdata_retagged.dac.txt");
	private final File EXPECTED_SUMMARY = new File (rootDir, "hek_5_cell_2_snp_testdata_retagged.dac.summary.txt");
	private final File EXPECTED_OUTFILE_CLUSTER = new File (rootDir, "hek_5_cell_2_snp_testdata_retagged.cluster.dac.txt");
	private final File EXPECTED_ALLELE_FREQ_OUTPUT= new File (rootDir, "hek_5_cell_2_snp_testdata_retagged.dac.allele_freq.txt");

	// new data for including untagged reads
	private final File untaggedReadsBAM = new File (rootDir, "includeUntaggedReads.bam");
	private final File untaggedReadsVCF = new File (rootDir, "includeUntaggedReads.vcf.gz");
	private final File untaggedReadsCellBarcodes = new File (rootDir, "includeUntaggedReads.cell_barcodes.txt");
	 
	private final File EXPECTED_EMPTY_SUMMARY=new File (rootDir, "includeUntaggedReadsFalse.summary.txt");
	private final File EXPECTED_UNSPLIT_SUMMARY=new File (rootDir, "includeUntaggedReadsTrue.summary.txt");
	private final File EXPECTED_SPLIT_SUMMARY=new File (rootDir, "includeUntaggedReadsTrueSplitByStrand.summary.txt");
	
	
	// See MultiCellDigitalAlleleCountsIteratorTest:bigDataTestED1 for more fine grained validation of the large BAM file data.
	@Test 
	public void testDoWork() throws IOException {
		GatherDigitalAlleleCounts g = new GatherDigitalAlleleCounts();
		g.CELL_BC_FILE=cellBCFile;
		g.INPUT=Collections.singletonList(largeBAMFile);
		g.VCF=vcfFile;
		g.HET_SNPS_ONLY=true;
		g.SAMPLE_FILE=sampleFile;
		g.OUTPUT = File.createTempFile("GatherDigitalAlleleCountsTest.", ".dac.txt");
		g.OUTPUT.deleteOnExit();
		g.SUMMARY = File.createTempFile("GatherDigitalAlleleCountsTest.", ".dac.summary.txt");
		g.SUMMARY.deleteOnExit();
		g.ALLELE_FREQUENCY_OUTPUT=File.createTempFile("GatherDigitalAlleleCountsTest.", ".dac.allele_freq.txt");
		g.ALLELE_FREQUENCY_OUTPUT.deleteOnExit();
		
		int ret = g.doWork();
		Assert.assertTrue(ret==0);

		Assert.assertTrue(FileUtils.contentEquals(g.OUTPUT, EXPECTED_OUTFILE));
		Assert.assertTrue(FileUtils.contentEquals(g.SUMMARY, EXPECTED_SUMMARY));
		Assert.assertTrue(FileUtils.contentEquals(g.ALLELE_FREQUENCY_OUTPUT, EXPECTED_ALLELE_FREQ_OUTPUT));
	}

	@Test
	public void testDoWorkCluster() throws IOException {
		GatherDigitalAlleleCounts g = new GatherDigitalAlleleCounts();
		g.CELL_BC_FILE=cellBCFile;
		g.INPUT=Collections.singletonList(largeBAMFile);
		g.VCF=vcfFile;
		g.HET_SNPS_ONLY=true;
		g.SAMPLE_FILE=sampleFile;		
		g.CLUSTER_FILE=clusterFile;
		File out = File.createTempFile("GatherDigitalAlleleCountsTest.", ".dac.txt");
		File summary = File.createTempFile("GatherDigitalAlleleCountsTest.", ".dac.summary.txt");
		File cluster = File.createTempFile("GatherDigitalAlleleCountsTest.", ".dac.txt");
		summary.deleteOnExit();
		cluster.deleteOnExit();
		out.deleteOnExit();
		g.OUTPUT=out;
		g.SUMMARY=summary;
		g.CLUSTER_OUTPUT=cluster;
		int ret = g.doWork();
		Assert.assertTrue(ret==0);

		Assert.assertTrue(FileUtils.contentEquals(cluster, EXPECTED_OUTFILE_CLUSTER));
		Assert.assertTrue(FileUtils.contentEquals(out, EXPECTED_OUTFILE));
		Assert.assertTrue(FileUtils.contentEquals(summary, EXPECTED_SUMMARY));

	}		
	
	@Test
	public void testIncludeUntaggedReadsNull() throws IOException {
		// without the new argument, the output is empty because all reads are intergenic.
		GatherDigitalAlleleCounts g = new GatherDigitalAlleleCounts();
		g.CELL_BC_FILE=untaggedReadsCellBarcodes;
		g.INPUT=Collections.singletonList(untaggedReadsBAM);
		g.VCF=untaggedReadsVCF;
		g.HET_SNPS_ONLY=false;
		g.SAMPLEs=Arrays.asList("AN10162");
		g.SUMMARY = File.createTempFile("GatherDigitalAlleleCountsTest.", ".dac.summary.txt");
		g.SUMMARY.deleteOnExit();
		int ret = g.doWork();
		Assert.assertTrue(ret==0);		
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_EMPTY_SUMMARY, g.SUMMARY));
	}
	
	@Test
	public void testIncludeUntaggedReads() throws IOException {
		// test with intergenic reads, contigs NOT split.
		GatherDigitalAlleleCounts g = new GatherDigitalAlleleCounts();
		g.CELL_BC_FILE=untaggedReadsCellBarcodes;
		g.INPUT=Collections.singletonList(untaggedReadsBAM);
		g.VCF=untaggedReadsVCF;
		g.HET_SNPS_ONLY=false;
		g.SAMPLEs=Arrays.asList("AN10162");
		g.SUMMARY = File.createTempFile("GatherDigitalAlleleCountsTest.", ".dac.summary.txt");
		g.SUMMARY.deleteOnExit();
		g.INCLUDE_UNTAGGED_READS=true;
		g.SPLIT_CONTIG_BY_STRAND=false;
		int ret = g.doWork();
		Assert.assertTrue(ret==0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_UNSPLIT_SUMMARY, g.SUMMARY));
		
	}
	
	@Test
	public void testIncludeUntaggedReadsSplit() throws IOException {
		// test with intergenic reads, contigs split.
		GatherDigitalAlleleCounts g = new GatherDigitalAlleleCounts();	
		g = new GatherDigitalAlleleCounts();
		g.CELL_BC_FILE=untaggedReadsCellBarcodes;
		g.INPUT=Collections.singletonList(untaggedReadsBAM);
		g.VCF=untaggedReadsVCF;
		g.HET_SNPS_ONLY=false;
		g.SAMPLEs=Arrays.asList("AN10162");
		g.SUMMARY = File.createTempFile("GatherDigitalAlleleCountsTest.", ".dac.summary.txt");
		g.SUMMARY.deleteOnExit();
		g.INCLUDE_UNTAGGED_READS=true;
		g.SPLIT_CONTIG_BY_STRAND=true;
		int ret = g.doWork();
		Assert.assertTrue(ret==0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_SPLIT_SUMMARY, g.SUMMARY));		
	}
	
	
	
	

}
