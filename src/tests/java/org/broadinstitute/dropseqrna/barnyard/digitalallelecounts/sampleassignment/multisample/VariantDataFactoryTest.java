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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample;

import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.AssignCellsToSamples;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.SampleGenotypeProbabilities;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.AllPairedSampleAssignmentsForCell;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.FindOptimalDonorMixture;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.GenotypeMatrix;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.SamplePairAssignmentForCell;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.VariantDataFactory;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

public class VariantDataFactoryTest {

	private final static String rootDir="testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/sampleassignment/multisample";
	
	private final File INPUT_BAM = new File(rootDir+"/TTTGCGCGGAGC:ATTGTTTAGGAG_retagged.bam");

	private final File VCF = new File (rootDir+"/TTTGCGCGGAGC:ATTGTTTAGGAG.vcf");

	private static final Log log = Log.getInstance(VariantDataFactoryTest.class);


	@Test(enabled=true)
	public void bigTest() {
		AssignCellsToSamples assigner = new AssignCellsToSamples();
		assigner.GQ_THRESHOLD=30;
		assigner.INPUT_BAM=Collections.singletonList(this.INPUT_BAM);
		assigner.VCF=this.VCF;
		assigner.NUM_BARCODES=10;
		assigner.RETAIN_MONOMORPIC_SNPS=true;

		final VCFFileReader vcfReader = new VCFFileReader(this.VCF, false);
		List<String> vcfSamples =SampleAssignmentVCFUtils.getVCFSamples(vcfReader, null);

		final IntervalList snpIntervals = SampleAssignmentVCFUtils.getSNPIntervals(this.VCF, vcfSamples, false, assigner.GQ_THRESHOLD, 0, null, null);
		final PeekableIterator<VariantContext> vcfIterator = SampleAssignmentVCFUtils.getVCFIterator(vcfReader, vcfSamples, false, assigner.GQ_THRESHOLD, 0.5, null, null);
		List<String> cellBarcodes = assigner.getCellBarcodes();
		PeekableIterator<List<SampleGenotypeProbabilities>> sampleGenotypeIterator = assigner.prepareIterator(snpIntervals, cellBarcodes);

		List<SampleGenotypeProbabilities> allProbs = new ArrayList<>();
		while (sampleGenotypeIterator.hasNext())
			allProbs.addAll(sampleGenotypeIterator.next());


		GenotypeMatrix genotypeMatrix = new GenotypeMatrix(vcfIterator, 30, new HashSet<>(vcfSamples));
		VariantDataFactory f = new VariantDataFactory("TTTGCGCGGAGC:ATTGTTTAGGAG", allProbs, genotypeMatrix, 0.1);
		FindOptimalDonorMixture fodm = new FindOptimalDonorMixture(f);

		AllPairedSampleAssignmentsForCell result = fodm.findBestDonorPair("HUES53", vcfSamples, false);
		SamplePairAssignmentForCell best = result.getBestAssignment();
		Assert.assertEquals(best.getSampleOneSingleLikelihood(),-9.80926731658451, 0.01);
		Assert.assertEquals(best.getSampleTwoSingleLikelihood(),-12.0398723515404, 0.01);
		Assert.assertEquals(best.getDoubletLikelihood(),-8.72621729107288, 0.01);
		Assert.assertEquals(best.getMixture(), 0.643, 0.01);
	}

	@Test(enabled=true)
	public void testMissingData() {
		AssignCellsToSamples assigner = new AssignCellsToSamples();
		assigner.GQ_THRESHOLD=30;
		assigner.INPUT_BAM=Collections.singletonList(this.INPUT_BAM);
		assigner.VCF=this.VCF;
		assigner.NUM_BARCODES=10;
		assigner.RETAIN_MONOMORPIC_SNPS=true;

		final VCFFileReader vcfReader = new VCFFileReader(this.VCF, false);
		List<String> vcfSamples =SampleAssignmentVCFUtils.getVCFSamples(vcfReader, null);

		final IntervalList snpIntervals = SampleAssignmentVCFUtils.getSNPIntervals(this.VCF, vcfSamples, false, assigner.GQ_THRESHOLD, 0, null, null);
		final PeekableIterator<VariantContext> vcfIterator = SampleAssignmentVCFUtils.getVCFIterator(vcfReader, vcfSamples, false, assigner.GQ_THRESHOLD, 0.5, null, null);
		List<String> cellBarcodes = assigner.getCellBarcodes();
		PeekableIterator<List<SampleGenotypeProbabilities>> sampleGenotypeIterator = assigner.prepareIterator(snpIntervals, cellBarcodes);

		List<SampleGenotypeProbabilities> allProbs = new ArrayList<>();
		while (sampleGenotypeIterator.hasNext())
			allProbs.addAll(sampleGenotypeIterator.next());


		GenotypeMatrix genotypeMatrix = new GenotypeMatrix(vcfIterator, 30, new HashSet<>(vcfSamples));
		VariantDataFactory f = new VariantDataFactory("TTTGCGCGGAGC:ATTGTTTAGGAG", allProbs, genotypeMatrix, 0.1, true, null, null, null);
		FindOptimalDonorMixture fodm = new FindOptimalDonorMixture(f);

		AllPairedSampleAssignmentsForCell result = fodm.findBestDonorPair("HUES53", vcfSamples, false);
		SamplePairAssignmentForCell best = result.getBestAssignment();
		Assert.assertEquals(best.getSampleOneSingleLikelihood(),-9.80926731658451, 0.01);
		Assert.assertEquals(best.getSampleTwoSingleLikelihood(),-12.0398723515404, 0.01);
		Assert.assertEquals(best.getDoubletLikelihood(),-8.72621729107288, 0.01);
		Assert.assertEquals(best.getMixture(), 0.643, 0.01);
	}



}
