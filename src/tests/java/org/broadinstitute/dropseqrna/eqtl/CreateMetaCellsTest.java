package org.broadinstitute.dropseqrna.eqtl;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class CreateMetaCellsTest {
	private final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/eqtl");
	private final File DONOR_MAP = new File(TEST_DATA_DIR, "donor_map.txt");
	private final File EXPECTED_METACELLS_SUM = new File(TEST_DATA_DIR, "metacells.txt.gz");
	
	private final File EXPECTED_METACELLS_MEAN = new File(TEST_DATA_DIR, "metacells.mean.txt.gz");
	private final File EXPECTED_METACELLS_MEDIAN = new File(TEST_DATA_DIR, "metacells.median.txt.gz");
	
	private final File EXPECTED_METACELLS_METRICS = new File (TEST_DATA_DIR, "meta_cell_metrics");
	private final File EXPECTED_SINGLE_METACELL = new File(TEST_DATA_DIR, "single_metacell.txt.gz");
	private final File DGE = new File("testdata/org/broadinstitute/transcriptome/barnyard/digitalexpression/test_with_header3.unpaired.dge.txt.gz");
	private final File SORTED_DGE = new File("testdata/org/broadinstitute/transcriptome/barnyard/digitalexpression/test_with_header3.unpaired.sorted.dge.txt.gz");
	private final File EXPECTED_METACELLS_BY_CLUSTER= new File (TEST_DATA_DIR, "cluster_1_2.meta_cells.txt");
	private final File EXPECTED_METACELLS_BY_CLUSTER_SORTED= new File (TEST_DATA_DIR, "cluster_1_2.meta_cells.sorted.txt");
	private final File MERGED_DGE_HEADER_FILE=new File(TEST_DATA_DIR, "test_with_header3.unpaired.dge.dge_header.txt");
	private final File CLUSTER_ASSIGNMENT_FILE=new File(TEST_DATA_DIR, "test_with_header3.unpaired.dge.assign.txt");
	
	
	
	@Test ()
	public void testMetacellsWithICAClusters() throws IOException {
		final CreateMetaCells clp = makeBasicTestClp();
		Assert.assertEquals(clp.doWork(), 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_METACELLS_BY_CLUSTER, clp.OUTPUT));
		
	}

	private CreateMetaCells makeBasicTestClp() throws IOException {
		final CreateMetaCells clp = new CreateMetaCells();
		clp.INPUT = DGE;
		clp.DONOR_MAP = DONOR_MAP;
		clp.CLUSTER_ASSIGNMENTS_FILE=this.CLUSTER_ASSIGNMENT_FILE;
		clp.CLUSTER_ASSIGNMENT=Arrays.asList("Cluster_1", "Cluster_2");
		clp.MERGED_DGE_HEADER_FILE=MERGED_DGE_HEADER_FILE;
		clp.OUTPUT = File.createTempFile("testMetacellsWithICAClusters.",".metacells.txt.gz");
		clp.OUTPUT.deleteOnExit();
		return clp;
	}

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void testInputSortRequirementNegative() throws IOException {
		final CreateMetaCells clp = makeBasicTestClp();
		clp.GENE_SORT = CreateMetaCells.GeneSort.REQUIRE_SORTED;
		clp.doWork();

	}

	@Test
	public void testInputSortRequirementPositive() throws IOException {
		final CreateMetaCells clp = makeBasicTestClp();
		clp.GENE_SORT = CreateMetaCells.GeneSort.REQUIRE_SORTED;
		clp.INPUT = SORTED_DGE;
		Assert.assertEquals(clp.doWork(), 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_METACELLS_BY_CLUSTER_SORTED, clp.OUTPUT));
	}


	@Test(dataProvider = "testMetacellsWithGeneSortDataProvider")
	public void testMetacellsWithGeneSort(Integer maxRecordsInRam) throws IOException {
		final CreateMetaCells clp = makeBasicTestClp();
		clp.GENE_SORT = CreateMetaCells.GeneSort.SORT;
		if (maxRecordsInRam != null) {
			clp.MAX_RECORDS_IN_RAM = maxRecordsInRam;
		}
		Assert.assertEquals(clp.doWork(), 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_METACELLS_BY_CLUSTER_SORTED, clp.OUTPUT));
	}

	// One invocation that is small enough that SortingCollection spills to disk.
	@DataProvider(name = "testMetacellsWithGeneSortDataProvider")
	public Object[][] testMetacellsWithGeneSortDataProvider() {
		return new Object[][]{
				{null},
				{3}
		};
	}

	@Test
	public void testGetCellBarcodesInClusters() {
		final CreateMetaCells clp = new CreateMetaCells();
		Set<String> clusterLabels = new HashSet<>(Arrays.asList("Cluster_1", "Cluster_2"));
		
		List<String> cellBC = clp.getCellBarcodesInClusters("TEST", clusterLabels, CLUSTER_ASSIGNMENT_FILE);
		Assert.assertEquals(cellBC.size(), 69);
		
		Set<String> clusterLabels1 = new HashSet<>(Collections.singletonList("Cluster_1"));
		
		List<String> cellBC2 = clp.getCellBarcodesInClusters("TEST", clusterLabels1, CLUSTER_ASSIGNMENT_FILE);
		Assert.assertNotNull(cellBC2);
		Assert.assertEquals(cellBC2.size(), 32);						
	}
	
	@Test 
	public void testGetPrefixForExperiment() {
		final CreateMetaCells clp = new CreateMetaCells();
		String uei = clp.getUEI(this.DGE);
		String expectedPrefix="TEST";
		String prefix = clp.getPrefixForExperiment (MERGED_DGE_HEADER_FILE, null, uei);
		Assert.assertEquals(prefix, expectedPrefix);
		
		String expectedPrefix2="FOO";
		String prefix2 = clp.getPrefixForExperiment (MERGED_DGE_HEADER_FILE, expectedPrefix2, uei);
		Assert.assertEquals(prefix2, expectedPrefix2);
		
	}
	
	@Test 
	public void testGetUEI() {
		final CreateMetaCells clp = new CreateMetaCells();
		String uei = clp.getUEI(this.DGE);
		String expectedUEI= "d35Ngn2plusGlia_Ngn2gliaD35B_E8_merged";
		Assert.assertEquals(uei, expectedUEI);
		
	}
	
	
	
	
	@Test
	public void testMetacellsFromDonorMapSum() throws IOException {
		final CreateMetaCells clp = new CreateMetaCells();
		clp.INPUT = DGE;
		clp.DONOR_MAP = DONOR_MAP;
		clp.OUTPUT = File.createTempFile("testMetacellsFromDonorMap.",".metacells.txt.gz");
		clp.METRICS =File.createTempFile("testMetacellsFromDonorMap.",".metrics.txt");
		clp.OUTPUT.deleteOnExit();
		clp.METRICS.deleteOnExit();
		clp.STRATEGY= DonorMergeStrategy.Sum;
		Assert.assertEquals(clp.doWork(), 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_METACELLS_SUM, clp.OUTPUT));
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_METACELLS_METRICS, clp.METRICS));
	}
	
//  test merging by various summaries in R.	
//	dgeFile="/Users/nemesh/dropseqrna/transcriptome_java/testdata/org/broadinstitute/transcriptome/barnyard/digitalexpression/test_with_header3.unpaired.dge.txt.gz"
//	donorMap="/Users/nemesh/dropseqrna/transcriptome_java/testdata/org/broadinstitute/dropseq/private/eqtl/donor_map.txt"
//	dge=read.table(dgeFile, header=T, stringsAsFactors = F)
//	dm=read.table(donorMap, header=T, stringsAsFactors = F)
//	cells=dm[dm$bestSample==1,]$cell
//	idx=match(cells, colnames(dge))
//	oneDonor=dge[,idx]
//	df=data.frame(gene=dge$GENE, sum=apply(oneDonor, 1, sum), mean=apply(oneDonor, 1, mean), median=apply(oneDonor, 1, median), stringsAsFactors = F)

	@Test
	public void testMetacellsFromDonorMapMean() throws IOException {
		final CreateMetaCells clp = new CreateMetaCells();
		clp.INPUT = DGE;
		clp.DONOR_MAP = DONOR_MAP;
		clp.OUTPUT = File.createTempFile("testMetacellsFromDonorMap.",".metacells.txt.gz");
		clp.METRICS =File.createTempFile("testMetacellsFromDonorMap.",".metrics.txt");
		clp.STRATEGY=DonorMergeStrategy.Mean;
		clp.INTEGER_FORMAT=false;
		clp.OUTPUT.deleteOnExit();
		clp.METRICS.deleteOnExit();
		Assert.assertEquals(clp.doWork(), 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_METACELLS_MEAN, clp.OUTPUT));
	}
	
	@Test
	public void testMetacellsFromDonorMapMedian() throws IOException {
		final CreateMetaCells clp = new CreateMetaCells();
		clp.INPUT = DGE;
		clp.DONOR_MAP = DONOR_MAP;
		clp.OUTPUT = File.createTempFile("testMetacellsFromDonorMap.",".metacells.txt.gz");
		clp.METRICS =File.createTempFile("testMetacellsFromDonorMap.",".metrics.txt");
		clp.STRATEGY=DonorMergeStrategy.Median;
		clp.INTEGER_FORMAT=false;
		clp.OUTPUT.deleteOnExit();
		clp.METRICS.deleteOnExit();
		Assert.assertEquals(clp.doWork(), 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_METACELLS_MEDIAN, clp.OUTPUT));
	}
	


	@Test
	public void testMetacellsWithSingleLabel() throws IOException {
		final CreateMetaCells clp = new CreateMetaCells();
		clp.INPUT = DGE;
		clp.SINGLE_METACELL_LABEL = "metacell_label";
		clp.OUTPUT = File.createTempFile("testMetacellsWithSingleLabel.",".metacells.txt.gz");
		clp.OUTPUT.deleteOnExit();
		Assert.assertEquals(clp.doWork(), 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_SINGLE_METACELL, clp.OUTPUT));
	}
	
	
	
}
