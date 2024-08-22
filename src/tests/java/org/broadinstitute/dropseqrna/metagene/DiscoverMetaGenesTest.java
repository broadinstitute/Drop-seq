package org.broadinstitute.dropseqrna.metagene;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;

import org.broadinstitute.dropseqrna.barnyard.DigitalExpression;
import org.broadinstitute.dropseqrna.metagene.DiscoverMetaGenes;
import org.broadinstitute.dropseqrna.metagene.MetaGene;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class DiscoverMetaGenesTest {

	
	
	/*
	 *  SMN EXPECTED DATA:
     *  GENE AACAACCAGTTACGGG AAGACAATCCGCAACG
     *  SMN1                1                2
     *  SMN2                1                1
     *  
     *  Additionally there's one read that is mapped at MQ=3, and resolves to an additional metaGene UMI on AAGACAATCCGCAACG.
	 */
	private static File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/metagene");
	private static final File IN_SMN_FILE = new File(TEST_DATA_DIR, "SMN.bam");
	private static final File IN_SMN_CELL_BARCODE_FILE = new File(TEST_DATA_DIR, "SMN_cellbarcodes.txt");
	private static final File EXPECTED_SMN_REPORT = new File(TEST_DATA_DIR, "SMN_expected_report.txt");
	private static final File EXPECTED_SMN_METAGENE_DGE = new File(TEST_DATA_DIR, "SMN_metagene.dge.txt");
	private static final File EXPECTED_SMN_DGE = new File(TEST_DATA_DIR, "SMN.dge.txt");
	
	private static final File SMN_EXTENDED_MODEL_NAME_FILE=new File(TEST_DATA_DIR, "SMN_extended_name.txt");
	private static final File SMN_EXTENDED_MODEL_NAME_METAGENE_DGE = new File(TEST_DATA_DIR, "SMN_metagene_extended.dge.txt");
	
	private static final File IN_SMN_FILE_BIG = new File(TEST_DATA_DIR, "BA46_downsampled_SMN.bam");
	private static final File IN_SMN_CELL_BARCODE_FILE_BIG = new File(TEST_DATA_DIR, "BA46_downsampled_SMN.selectedCellBarcodes.txt");
	private static final File EXPECTED_SMN_REPORT_BIG = new File(TEST_DATA_DIR, "BA46_downsampled_SMN.metagene_report.txt");
		
	@Test(enabled=true)
	public void testDetectSMN() throws IOException {		
		DiscoverMetaGenes dmg = new DiscoverMetaGenes();
		dmg.CELL_BC_FILE=IN_SMN_CELL_BARCODE_FILE;
		dmg.WRITE_ALL_READS=true;
		dmg.INPUT=IN_SMN_FILE;
		dmg.MIN_READ_MQ=2;		
		dmg.REPORT=File.createTempFile("DiscoverMetaGenes", ".report.txt");
		dmg.REPORT.deleteOnExit();
		dmg.OUTPUT=File.createTempFile("DiscoverMetaGenes", ".SMN.bam");
		dmg.OUTPUT.deleteOnExit();
		// take all meta genes since it's only the one.
		dmg.META_GENE_RATIO=0d;
		dmg.MAXIMUM_READ_EDIT_DISTANCE=null;
		TestUtils.setInflaterDeflaterIfMacOs();
		int result = dmg.doWork();
		Assert.assertEquals(0, result);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_SMN_REPORT, dmg.REPORT));
		
		// Extract DGE for the meta gene.  Should be AAGACAATCCGCAACG=1, AACAACCAGTTACGGG=0		
		File dgeOutput =File.createTempFile("SMN", ".digital_expression.txt");		
		dgeOutput.deleteOnExit();
		String [] dgeArgs = {"INPUT="+dmg.OUTPUT, "CELL_BC_FILE="+this.IN_SMN_CELL_BARCODE_FILE, "READ_MQ=2", "GENE_NAME_TAG="+dmg.METAGENE_NAME, 
				"GENE_STRAND_TAG="+dmg.METAGENE_STRAND, "GENE_FUNCTION_TAG="+dmg.METAGENE_FUNCTION, "OUTPUT="+dgeOutput};
											
		int dgeReturn = new DigitalExpression().instanceMain(dgeArgs);
		
		Assert.assertEquals(0, dgeReturn);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_SMN_METAGENE_DGE, dgeOutput));
		
		// Extract DGE for the unique genes.  
		File dgeOutputUnique =File.createTempFile("SMN", ".digital_expression.txt");		
		dgeOutputUnique.deleteOnExit();
		String [] dgeArgs2 = {"INPUT="+dmg.OUTPUT, "CELL_BC_FILE="+this.IN_SMN_CELL_BARCODE_FILE, "READ_MQ=10", "OUTPUT="+dgeOutputUnique};					
		int dgeReturn2 = new DigitalExpression().instanceMain(dgeArgs2);
		Assert.assertEquals(0, dgeReturn2);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_SMN_DGE, dgeOutputUnique));
	}
	
	
	/**
	 * Don't discover metagenes, use a known gene file instead.
	 * @throws IOException
	 */
	@Test(enabled=true)
	public void testSMNKnownGenes() throws IOException {		
		DiscoverMetaGenes dmg = new DiscoverMetaGenes();
		dmg.CELL_BC_FILE=IN_SMN_CELL_BARCODE_FILE;
		dmg.WRITE_ALL_READS=true;
		dmg.INPUT=IN_SMN_FILE;
		dmg.MIN_READ_MQ=2;		
		dmg.KNOWN_META_GENE_FILE=EXPECTED_SMN_REPORT;
		dmg.OUTPUT=File.createTempFile("DiscoverMetaGenes", ".SMN.bam");
		dmg.OUTPUT.deleteOnExit();
		// take all meta genes since it's only the one.
		dmg.META_GENE_RATIO=0d;
		dmg.MAXIMUM_READ_EDIT_DISTANCE=null;
		TestUtils.setInflaterDeflaterIfMacOs();
		int result = dmg.doWork();
		Assert.assertEquals(0, result);		
		
		// Extract DGE for the meta gene.  Should be AAGACAATCCGCAACG=1, AACAACCAGTTACGGG=0		
		File dgeOutput =File.createTempFile("SMN", ".digital_expression.txt");		
		dgeOutput.deleteOnExit();
		String [] dgeArgs = {"INPUT="+dmg.OUTPUT, "CELL_BC_FILE="+this.IN_SMN_CELL_BARCODE_FILE, "READ_MQ=2", "GENE_NAME_TAG="+dmg.METAGENE_NAME, 
				"GENE_STRAND_TAG="+dmg.METAGENE_STRAND, "GENE_FUNCTION_TAG="+dmg.METAGENE_FUNCTION, "OUTPUT="+dgeOutput};
											
		int dgeReturn = new DigitalExpression().instanceMain(dgeArgs);
		
		Assert.assertEquals(0, dgeReturn);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_SMN_METAGENE_DGE, dgeOutput));
		
		// Extract DGE for the unique genes.  
		File dgeOutputUnique =File.createTempFile("SMN", ".digital_expression.txt");		
		dgeOutputUnique.deleteOnExit();
		String [] dgeArgs2 = {"INPUT="+dmg.OUTPUT, "CELL_BC_FILE="+this.IN_SMN_CELL_BARCODE_FILE, "READ_MQ=10", "OUTPUT="+dgeOutputUnique};					
		int dgeReturn2 = new DigitalExpression().instanceMain(dgeArgs2);
		Assert.assertEquals(0, dgeReturn2);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_SMN_DGE, dgeOutputUnique));
	}
	
	/**
	 * Don't discover metagenes, use a known gene file instead.
	 * Use a metagene model that has more genes in the grouping than discovered here, output should report the metagene model not the discovered one.
	 * @throws IOException
	 */
	@Test(enabled=true)
	public void testSMNExtendedModel() throws IOException {		
		DiscoverMetaGenes dmg = new DiscoverMetaGenes();
		dmg.CELL_BC_FILE=IN_SMN_CELL_BARCODE_FILE;
		dmg.WRITE_ALL_READS=true;
		dmg.INPUT=IN_SMN_FILE;
		dmg.MIN_READ_MQ=2;		
		dmg.KNOWN_META_GENE_FILE=SMN_EXTENDED_MODEL_NAME_FILE;
		dmg.OUTPUT=File.createTempFile("DiscoverMetaGenes", ".SMN.bam");
		dmg.OUTPUT.deleteOnExit();
		// take all meta genes since it's only the one.
		dmg.META_GENE_RATIO=0d;
		dmg.MAXIMUM_READ_EDIT_DISTANCE=null;
		TestUtils.setInflaterDeflaterIfMacOs();
		int result = dmg.doWork();
		Assert.assertEquals(0, result);		
		
		// Extract DGE for the meta gene.  Should be AAGACAATCCGCAACG=1, AACAACCAGTTACGGG=0		
		File dgeOutput =File.createTempFile("SMN", ".digital_expression.txt");		
		dgeOutput.deleteOnExit();
		String [] dgeArgs = {"INPUT="+dmg.OUTPUT, "CELL_BC_FILE="+this.IN_SMN_CELL_BARCODE_FILE, "READ_MQ=2", "GENE_NAME_TAG="+dmg.METAGENE_NAME, 
				"GENE_STRAND_TAG="+dmg.METAGENE_STRAND, "GENE_FUNCTION_TAG="+dmg.METAGENE_FUNCTION, "OUTPUT="+dgeOutput};
											
		int dgeReturn = new DigitalExpression().instanceMain(dgeArgs);
		
		Assert.assertEquals(0, dgeReturn);
		Assert.assertTrue(TestUtils.testFilesSame(SMN_EXTENDED_MODEL_NAME_METAGENE_DGE, dgeOutput));
		
	}
	
	@Test
	public void testGetMetaGeneInKnownModel () {
		DiscoverMetaGenes dmg = new DiscoverMetaGenes();
		MetaGene m = new MetaGene(Arrays.asList("SMN1", "SMN2"));
		
		// exact match
		Collection <MetaGene> approvedMetaGenes = Arrays.asList(m);				
		MetaGene result = dmg.getMetaGeneInKnownModel(m, approvedMetaGenes);
		Assert.assertNotNull(result);
		
		// test contains match
		MetaGene k = new MetaGene(Arrays.asList("SMN1", "SMN2", "FOO"));
		approvedMetaGenes = Arrays.asList(k);		
		result = dmg.getMetaGeneInKnownModel(m, approvedMetaGenes);
		Assert.assertNotNull(result);
		
		// no match, but is a substring.
		MetaGene m2 = new MetaGene(Arrays.asList("MN1", "SMN2"));
		result = dmg.getMetaGeneInKnownModel(m2, approvedMetaGenes);
		Assert.assertNull(result);		
	}
	
	@Test
	public void testSMNManySelectedCells() throws IOException {
		DiscoverMetaGenes dmg = new DiscoverMetaGenes();
		dmg.CELL_BC_FILE=IN_SMN_CELL_BARCODE_FILE_BIG;
		dmg.INPUT=IN_SMN_FILE_BIG;
		dmg.MIN_READ_MQ=2;				
		dmg.REPORT=File.createTempFile("DiscoverMetaGenes", ".report");
		dmg.REPORT.deleteOnExit();
		dmg.OUTPUT=File.createTempFile("DiscoverMetaGenes", ".SMN.bam");
		dmg.OUTPUT.deleteOnExit();
		// take all meta genes since it's only the one.
		dmg.META_GENE_RATIO=0d;
		dmg.MAXIMUM_READ_EDIT_DISTANCE=null;
		TestUtils.setInflaterDeflaterIfMacOs();
		int result = dmg.doWork();
		Assert.assertEquals(0, result);		
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_SMN_REPORT_BIG, dmg.REPORT));
		
	}
	
	@Test
	public void testSMNAllCells() throws IOException {
		DiscoverMetaGenes dmg = new DiscoverMetaGenes();		
		dmg.INPUT=IN_SMN_FILE_BIG;
		dmg.MIN_READ_MQ=2;				
		dmg.OUTPUT=File.createTempFile("DiscoverMetaGenes", ".SMN.bam");
		dmg.OUTPUT.deleteOnExit();
		dmg.REPORT=File.createTempFile("DiscoverMetaGenes", ".report");
		dmg.REPORT.deleteOnExit();
		// take all meta genes since it's only the one.
		dmg.META_GENE_RATIO=0d;
		dmg.MAXIMUM_READ_EDIT_DISTANCE=null;
		TestUtils.setInflaterDeflaterIfMacOs();
		int result = dmg.doWork();
		Assert.assertEquals(0, result);		
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_SMN_REPORT_BIG, dmg.REPORT));
		
		
	}
	
	
	
	
}
