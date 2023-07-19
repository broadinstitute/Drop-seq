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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPInfoCollection;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.AssignCellsToSamples;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.CellCollectionSampleLikelihoodCollection;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.CellSampleLikelihoodCollection;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.SampleGenotypeProbabilities;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;

public class AssignCellsToSamplesTest {

	private final static String rootDir="testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/sampleassignment/";
	
	private final File INPUT_BAM = new File(rootDir+"/small_data_retagged.sam");

	private final File VCF = new File (rootDir+"/test_missing_data.vcf");
	private final File VCFC = new File (rootDir+"/test_missing_data.vcf.gz");
	private final File VCF_NON_CANONICAL=new File(rootDir+"/test_non_canonical_base.vcf.gz");
	
	private final File CONTAMINATION_FILE=new File (rootDir+"/small_data.contamination.txt");
	private final File MAF_ESIMATE_FILE=new File (rootDir+"/small_data.minor_allele_freq.txt");
	
	private final File EXPECTED_OUTPUT=new File (rootDir+"/small_data.donor_assignment.txt");
	private final File EXPECTED_VERBOSE_OUTPUT= new File (rootDir+"/small_data.donor_assignment.verbose.txt.gz");
	
	private final File ANSWER_KEY=new File (rootDir+"/small_data.donor_assignment.answer_key.txt");
	private final File EXPECTED_OUTPUT_ANSWER_KEY=new File (rootDir+"/small_data.donor_assignment_with_key.txt");
	
	private final File EXPECTED_CONTAM_OUTPUT= new File (rootDir+"/small_data_contam.donor_assignment.txt");
	private final File EXPECTED_CONTAM_VERBOSE_OUTPUT= new File (rootDir+"/small_data_contam.donor_assignment.verbose.txt.gz");
	
	private final File EXPECTED_MAXLIKE_OUTPUT=new File (rootDir+"/small_data.donor_assignment.maxLike_10.txt");	
	private final File EXPECTED_MAXLIKE_VERBOSE_OUTPUT= new File (rootDir+"/small_data.donor_assignment.maxLike_10.verbose.gz");
	
	
	@Test 
	// Two SNPs (rs3,rs4) have non-canonical bases (*,N).  
	// These SNPs should be filtered out, such that the data appears to be the original VCF.
	public void testNonCannonicalBase() throws IOException {
		AssignCellsToSamples assigner = new AssignCellsToSamples();
		assigner.GQ_THRESHOLD=30;
		assigner.INPUT_BAM=Collections.singletonList(this.INPUT_BAM);
		assigner.VCF=this.VCF_NON_CANONICAL;
		assigner.NUM_BARCODES=10; // more than enough
		assigner.OUTPUT=File.createTempFile("AssignCellsToSamples", ".output");
		assigner.VCF_OUTPUT=File.createTempFile("AssignCellsToSamples", ".vcf");
		assigner.BAM_OUTPUT=File.createTempFile("AssignCellsToSamples", ".informative.bam");
		assigner.VERBOSE_OUTPUT=File.createTempFile("AssignCellsToSamples", ".verbose.gz");
		assigner.VERBOSE_BEST_DONOR_OUTPUT=File.createTempFile("AssignCellsToSamples", ".best_verbose.gz");
		int result = assigner.doWork();
		Assert.assertEquals(result, 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_OUTPUT, assigner.OUTPUT));	
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_VERBOSE_OUTPUT, assigner.VERBOSE_OUTPUT));
		
	}
	
	
	@Test
	public void testEndToEnd() throws IOException {
		AssignCellsToSamples assigner = new AssignCellsToSamples();
		assigner.GQ_THRESHOLD=30;
		assigner.INPUT_BAM=Collections.singletonList(this.INPUT_BAM);
		assigner.VCF=this.VCFC;
		assigner.NUM_BARCODES=10; // more than enough
		assigner.OUTPUT=File.createTempFile("AssignCellsToSamples", ".output");
		assigner.VCF_OUTPUT=File.createTempFile("AssignCellsToSamples", ".vcf");
		assigner.BAM_OUTPUT=File.createTempFile("AssignCellsToSamples", ".informative.bam");
		assigner.VERBOSE_OUTPUT=File.createTempFile("AssignCellsToSamples", ".verbose.gz");
		assigner.VERBOSE_BEST_DONOR_OUTPUT=File.createTempFile("AssignCellsToSamples", ".best_verbose.gz");
		int result = assigner.doWork();
		Assert.assertEquals(result, 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_OUTPUT, assigner.OUTPUT));	
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_VERBOSE_OUTPUT, assigner.VERBOSE_OUTPUT));
		
	}
	
	
	
	// 
	@Test
	public void testEndToEndAnswerKey() throws IOException {
		AssignCellsToSamples assigner = new AssignCellsToSamples();
		assigner.GQ_THRESHOLD=30;
		assigner.INPUT_BAM=Collections.singletonList(this.INPUT_BAM);
		assigner.VCF=this.VCFC;
		assigner.NUM_BARCODES=10; // more than enough
		assigner.OUTPUT=File.createTempFile("AssignCellsToSamples", ".output");
		assigner.BAM_OUTPUT=File.createTempFile("AssignCellsToSamples", ".informative.bam");
		assigner.VERBOSE_OUTPUT=File.createTempFile("AssignCellsToSamples", ".verbose.gz");
		assigner.VERBOSE_BEST_DONOR_OUTPUT=File.createTempFile("AssignCellsToSamples", ".best_verbose.gz");
		assigner.ANSWER_KEY_FILE=ANSWER_KEY;
		
		int result = assigner.doWork();
		Assert.assertEquals(result, 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_OUTPUT_ANSWER_KEY, assigner.OUTPUT));	
		
		
	}

	@Test
	public void testEndToEndMaxLikelihood() throws IOException {
		AssignCellsToSamples assigner = new AssignCellsToSamples();
		assigner.GQ_THRESHOLD=30;
		assigner.INPUT_BAM=Collections.singletonList(this.INPUT_BAM);
		assigner.VCF=this.VCFC;
		assigner.NUM_BARCODES=10; // more than enough
		assigner.MAX_ERROR_RATE=0.1; // -1 per UMI when the donor is homozygous with the wrong allele is the "max" score.
		assigner.OUTPUT=File.createTempFile("AssignCellsToSamples", ".output");
		assigner.BAM_OUTPUT=File.createTempFile("AssignCellsToSamples", ".informative.bam");
		assigner.VERBOSE_OUTPUT=File.createTempFile("AssignCellsToSamples", ".verbose.gz");
		assigner.VERBOSE_BEST_DONOR_OUTPUT=File.createTempFile("AssignCellsToSamples", ".best_verbose.gz");
		int result = assigner.doWork();
		Assert.assertEquals(result, 0);
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_MAXLIKE_OUTPUT, assigner.OUTPUT)); 
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_MAXLIKE_VERBOSE_OUTPUT, assigner.VERBOSE_OUTPUT)); 
	}

	@Test (enabled=true)
	public void testEndToEndWithContamination() throws IOException {

		AssignCellsToSamples assigner = new AssignCellsToSamples();
		assigner.GQ_THRESHOLD=30;
		assigner.INPUT_BAM=Collections.singletonList(this.INPUT_BAM);
		assigner.VCF=this.VCFC;
		assigner.NUM_BARCODES=10; // more than enough
		assigner.OUTPUT=File.createTempFile("AssignCellsToSamples", ".output");
		assigner.BAM_OUTPUT=File.createTempFile("AssignCellsToSamples", ".informative.bam");
		assigner.VERBOSE_OUTPUT=File.createTempFile("AssignCellsToSamples", ".verbose.gz");
		assigner.VERBOSE_BEST_DONOR_OUTPUT=File.createTempFile("AssignCellsToSamples", ".best_verbose.gz");
		assigner.CELL_CONTAMINATION_ESTIMATE_FILE=CONTAMINATION_FILE;
		assigner.ALLELE_FREQUENCY_ESTIMATE_FILE=MAF_ESIMATE_FILE;
		
		int result = assigner.doWork();
		Assert.assertEquals(result, 0);		
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_CONTAM_OUTPUT, assigner.OUTPUT));	
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_CONTAM_VERBOSE_OUTPUT, assigner.VERBOSE_OUTPUT));
		
	}
	
	// set up for real world data to test a more rare bug.
	@Test(enabled=false)
	public void testfilterVerboseToBest() {
		File donorAssignments = new File ("/downloads/dropulation_global_likelihood/d21Ngn2plusGlia_H10.new.donor_assignments.txt");
		File verboseDonorAssignments = new File ("/downloads/dropulation_global_likelihood/d21Ngn2plusGlia_H10.new.verbose.donor_assignments.txt.gz");
		File out = new File ("/downloads/dropulation_global_likelihood/d21Ngn2plusGlia_H10.new.best_donor_verbose.donor_assignments.txt.gz");
		
		CellCollectionSampleLikelihoodCollection likelihoodCollection = CellCollectionSampleLikelihoodCollection.parseFile(donorAssignments);
		
		AssignCellsToSamples.filterVerboseToBest(verboseDonorAssignments, out, likelihoodCollection);
		
	}
	
	@Test
	public void testFilterVerboseToBest () throws IOException {
		File singleDonorLikelihoodFile = new File (rootDir + "/singleDonorAssignments_10_cells.txt");
		File verboseOutput = new File (rootDir + "/singleDonorAssignments_10_cells.verbose.txt.gz");		
		File tempOutput=File.createTempFile("testFilterVerboseToBest", ".txt.gz");
		File EXPECTED_OUTPUT = new File (rootDir + "/singleDonorAssignments_10_cells.verbose_best.txt.gz");
		
		tempOutput.deleteOnExit();
		CellCollectionSampleLikelihoodCollection cslc = CellCollectionSampleLikelihoodCollection.parseFile(singleDonorLikelihoodFile);
		AssignCellsToSamples.filterVerboseToBest(verboseOutput, tempOutput, cslc);
		
		Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_OUTPUT, tempOutput));
		
	}
	

	
	@Test
	// don't add missing values.
	public void testprocessSNP () {

		AssignCellsToSamples assigner = new AssignCellsToSamples();  
		assigner.GQ_THRESHOLD=30;
		assigner.INPUT_BAM=Collections.singletonList(this.INPUT_BAM);
		assigner.VCF=this.VCF;
		assigner.NUM_BARCODES=10;

		final VCFFileReader vcfReader = new VCFFileReader(this.VCF, false);
		List<String> vcfSamples =assigner.getVCFSamples(vcfReader);
		final SNPInfoCollection snpIntervals = SampleAssignmentVCFUtils.getSNPInfoCollection(this.VCF, vcfSamples, false, assigner.GQ_THRESHOLD, assigner.FRACTION_SAMPLES_PASSING, null, null, false, null);

		final PeekableIterator<VariantContext> vcfIterator = SampleAssignmentVCFUtils.getVCFIterator(vcfReader, vcfSamples, assigner.RETAIN_MONOMORPIC_SNPS, assigner.GQ_THRESHOLD, assigner.FRACTION_SAMPLES_PASSING, null, null);		
		List<String> cellBarcodes= assigner.getCellBarcodes();
		PeekableIterator<List<SampleGenotypeProbabilities>> sampleGenotypeIterator = assigner.prepareIterator(snpIntervals, cellBarcodes);

		CellCollectionSampleLikelihoodCollection likelihoodCollection = new CellCollectionSampleLikelihoodCollection(0.1, null, vcfSamples);
		// Set<String> sampleNames = new HashSet<String>();
		// sampleNames.add("FOO");

		while (vcfIterator.hasNext() && sampleGenotypeIterator.hasNext()) {
			VariantContext vc = vcfIterator.next();
			Set<String> sampleNames = vc.getSampleNames();
			List<SampleGenotypeProbabilities> pileUpsSingleSNP = sampleGenotypeIterator.next();
			likelihoodCollection=assigner.processSNP(vc, pileUpsSingleSNP, sampleNames, likelihoodCollection, false);
		}

		// [snp [HUMAN_1:76227022-76227022	-	rs1] cell [ATCAGGGACAGA] [A, A] [37, 8], snp [HUMAN_1:76227022-76227022	-	rs1] cell [TTGCCTTACGCG] [A, G] [27, 27]]
		// FOO: A/A, BAR: A/A, FOO2: A/G, BAR2:G/G, MISSING: NA
		CellSampleLikelihoodCollection cell1 = likelihoodCollection.getLikelihoodCollection("ATCAGGGACAGA");
		Double like=cell1.getLoglikelihood("FOO");
		Assert.assertEquals(like, -0.6935748, 0.000001);
		like = cell1.getLoglikelihood("BAR");
		Assert.assertEquals(like, -0.6935748, 0.000001);
		like = cell1.getLoglikelihood("FOO2");
		Assert.assertEquals(like, -2.60206, 0.000001);
		like = cell1.getLoglikelihood("BAR2");
		Assert.assertEquals(like, -4, 0.000001);
		like = cell1.getLoglikelihood("MISSING");
		Assert.assertEquals(like, 0d);

		CellSampleLikelihoodCollection cell2 = likelihoodCollection.getLikelihoodCollection("TTGCCTTACGCG");
		like=cell2.getLoglikelihood("FOO");
		Assert.assertEquals(like, -1.647817, 0.000001);
		like = cell2.getLoglikelihood("BAR");
		Assert.assertEquals(like, -1.647817, 0.000001);
		like = cell2.getLoglikelihood("FOO2");
		Assert.assertEquals(like, -0.693575, 0.000001);
		like = cell2.getLoglikelihood("BAR2");
		Assert.assertEquals(like, -1.137272, 0.000001);
		like = cell2.getLoglikelihood("MISSING");
		Assert.assertEquals(like, 0d);

	}

	@Test
	// don't add missing values.
	public void testprocessSNPWithMissingValues () {

		AssignCellsToSamples assigner = new AssignCellsToSamples();
		assigner.GQ_THRESHOLD=30;
		assigner.INPUT_BAM=Collections.singletonList(this.INPUT_BAM);
		assigner.VCF=this.VCF;
		assigner.NUM_BARCODES=10;

		final VCFFileReader vcfReader = new VCFFileReader(this.VCF, false);
		List<String> vcfSamples =assigner.getVCFSamples(vcfReader);
		final SNPInfoCollection snpIntervals = SampleAssignmentVCFUtils.getSNPInfoCollection(this.VCF, vcfSamples, false, assigner.GQ_THRESHOLD, assigner.FRACTION_SAMPLES_PASSING, null, null, false, null);
		final PeekableIterator<VariantContext> vcfIterator = SampleAssignmentVCFUtils.getVCFIterator(vcfReader, vcfSamples, assigner.RETAIN_MONOMORPIC_SNPS, assigner.GQ_THRESHOLD, assigner.FRACTION_SAMPLES_PASSING, null, null);
		List<String> cellBarcodes= assigner.getCellBarcodes();
		
		PeekableIterator<List<SampleGenotypeProbabilities>> sampleGenotypeIterator = assigner.prepareIterator(snpIntervals, cellBarcodes);

		CellCollectionSampleLikelihoodCollection likelihoodCollection = new CellCollectionSampleLikelihoodCollection(0.1, null, vcfSamples);
		// Set<String> sampleNames = new HashSet<String>();
		// sampleNames.add("FOO");

		while (vcfIterator.hasNext() && sampleGenotypeIterator.hasNext()) {
			VariantContext vc = vcfIterator.next();
			Set<String> sampleNames = vc.getSampleNames();
			List<SampleGenotypeProbabilities> pileUpsSingleSNP = sampleGenotypeIterator.next();
			likelihoodCollection=assigner.processSNP(vc, pileUpsSingleSNP, sampleNames, likelihoodCollection, true);
		}

		// [snp [HUMAN_1:76227022-76227022	-	rs1] cell [ATCAGGGACAGA] [A, A] [37, 8], snp [HUMAN_1:76227022-76227022	-	rs1] cell [TTGCCTTACGCG] [A, G] [27, 27]]
		// FOO: A/A, BAR: A/A, FOO2: A/G, BAR2:G/G, MISSING: NA
		CellSampleLikelihoodCollection cell1 = likelihoodCollection.getLikelihoodCollection("ATCAGGGACAGA");
		Double like=cell1.getLoglikelihood("FOO");
		Assert.assertEquals(like, -0.6935748, 0.000001);
		like = cell1.getLoglikelihood("BAR");
		Assert.assertEquals(like, -0.6935748, 0.000001);
		like = cell1.getLoglikelihood("FOO2");
		Assert.assertEquals(like, -2.60206, 0.000001);
		like = cell1.getLoglikelihood("BAR2");
		Assert.assertEquals(like, -4, 0.000001);
		like = cell1.getLoglikelihood("MISSING");
		// missing value is set?
		Assert.assertNotSame(like, null);
		// value is correct?
		Assert.assertEquals(like, -1.489455, 0.000001);

		CellSampleLikelihoodCollection cell2 = likelihoodCollection.getLikelihoodCollection("TTGCCTTACGCG");
		like=cell2.getLoglikelihood("FOO");
		Assert.assertEquals(like, -1.647817, 0.000001);
		like = cell2.getLoglikelihood("BAR");
		Assert.assertEquals(like, -1.647817, 0.000001);
		like = cell2.getLoglikelihood("FOO2");
		Assert.assertEquals(like, -0.693575, 0.000001);
		like = cell2.getLoglikelihood("BAR2");
		Assert.assertEquals(like, -1.137272, 0.000001);
		like = cell2.getLoglikelihood("MISSING");
		Assert.assertNotSame(like, null);
		// value is correct?
		Assert.assertEquals(like, -0.9295926, 0.000001);
	}







	@Test
	public void testcompareRecords1() {
		// sgp before vc.
		SampleGenotypeProbabilities sgp = new SampleGenotypeProbabilities(new Interval("chr1", 5, 5), "cell");
		VariantContext vc = getVariantContext("chr1", 10);
		AssignCellsToSamples assign = new AssignCellsToSamples();
		int cmp = assign.compareRecords(sgp, vc, null);
		Assert.assertTrue(cmp<0);
	}

	@Test
	public void testcompareRecords2() {
		// vc before sgp
		SampleGenotypeProbabilities sgp = new SampleGenotypeProbabilities(new Interval("chr1", 10, 10), "cell");
		VariantContext vc = getVariantContext("chr1", 5);
		// sgp before vc.
		AssignCellsToSamples assign = new AssignCellsToSamples();
		int cmp = assign.compareRecords(sgp, vc, null);
		Assert.assertTrue(cmp>0);
	}


	@Test
	public void testcompareRecords3() {
		SampleGenotypeProbabilities sgp = new SampleGenotypeProbabilities(new Interval("chr1", 10, 10), "cell");
		VariantContext vc = getVariantContext("chr1", 10);
		// sgp equal to vc.
		AssignCellsToSamples assign = new AssignCellsToSamples();
		int cmp = assign.compareRecords(sgp, vc, null);
		Assert.assertTrue(cmp==0);
	}

	private VariantContext getVariantContext (final String chr, final int start) {
		VariantContextBuilder b = new VariantContextBuilder();
		String [] all = {"A", "G"};
		b.alleles(all);
		b.start(start);
		b.stop(start);
		b.chr(chr);
		return b.make();
	}
}
