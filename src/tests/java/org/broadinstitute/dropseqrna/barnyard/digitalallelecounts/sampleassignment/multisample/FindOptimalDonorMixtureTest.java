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

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import org.broadinstitute.dropseqrna.barnyard.BarcodeListRetrieval;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPInfoCollection;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.AssignCellsToSamples;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.SampleGenotypeProbabilities;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.AllPairedSampleAssignmentsForCell;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.DetectDoublets;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.FindOptimalDonorMixture;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.GenotypeMatrix;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.SamplePairAssignmentForCell;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.VariantData;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.VariantDataCollection;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.VariantDataFactory;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

public class FindOptimalDonorMixtureTest {

	private final static String rootDir="testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/sampleassignment/multisample";
	
	private final File INPUT_BAM = new File(rootDir+"/TTTGCGCGGAGC:ATTGTTTAGGAG_retagged.bam");

	private final File VCF = new File (rootDir+"/TTTGCGCGGAGC:ATTGTTTAGGAG.vcf");
	
	//private final File VCF = new File (
	//		"/downloads/testVCF.vcf");

	private static final Log log = Log.getInstance(FindOptimalDonorMixtureTest.class);

	private String [] allDonors = {"NKX", "HUES53", "HUES62", "HUES63", "HUES64"}; 
	private String [] donors = {"HUES53", "NKX"};


	@Test(enabled=true)
	/**
	 * The equivalent in R:
	 * d=data.table(data.frame(pos=c("1:1", "1:2", "1:3", "1:4", "1:5"), refCount=c(2, 1,1,0,1), altCount=c(0,1,1,2,1), sampleOneGenotype=c("ref", "ref", "ref", "het","ref"), sampleTwoGenotype=c("alt", "het", "ref", "alt", "het"), stringsAsFactors=F))
	 * calculateLikelihoodsManyIterations(d, qualityScore=0.9, "one", "two", interval=c(0,1))
	 * sampleOneRatio jointSampleLikelihood       one       two 		num_paired_snps num_paired_informative_snps
     *	 0.5540364             -3.306935 		-3.830847 -4.341392               5                           4
	 */
	public void smallTest1() {
		// build up some reads
		List<VariantData> list = new ArrayList<>();

		list.add(new VariantData(new Interval ("1", 1, 1), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HOM_VAR, new char [] {'A', 'A'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 2, 2), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HET, new char [] {'A', 'T'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 3, 3), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HOM_REF, new char [] {'A', 'T'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 4, 4), 'A', 'T', GenotypeType.HET, GenotypeType.HOM_VAR, new char [] {'T', 'T'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 5, 5), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HET, new char [] {'A', 'T'}, new int [] {10,10}));

		VariantDataCollection vcd = new VariantDataCollection(list, "one", "two", "testCell");
		double likeS1 = vcd.value(1);
		double likeS2 = vcd.value(0);
		SamplePairAssignmentForCell s = vcd.optimizeMixture();
		double likeMixed= s.getDoubletLikelihood();
		double mixture = s.getMixture();
		// the sample mixture is slightly different, but close enough.  It seems like as I increase the number of observations (beyond just this handful) they converge.
		// this may be due to the fact that I stop convergence at 0.001 in the java version, so it can be off by a little bit.

		Assert.assertEquals(likeS1, -3.830847, 0.001);
		Assert.assertEquals(likeS2, -4.341392, 0.001);
		Assert.assertEquals(likeMixed, -3.306935, 0.001);
		Assert.assertEquals(mixture, 0.5540364, 0.015);
		Assert.assertEquals(s.getImpossibleAllelesSampleOne(),3);
		Assert.assertEquals(s.getImpossibleAllelesSampleTwo(),3);
	}

	/*The equivalent in R:
	 * d=data.table(data.frame(pos=c("1:1", "1:2", "1:3", "1:4", "1:5"), refCount=c(2, 1,1,0,1), altCount=c(0,1,1,2,1), sampleOneGenotype=c("ref", "ref", "ref", "het","ref"), sampleTwoGenotype=c("alt", "het", "ref", "alt", "het"), stringsAsFactors=F))
	 * calculateLikelihoodsOneIteration(0.5, d, qualityScore)
	 * [1] -3.313183
	 */
	@Test
	public void smallTestForcedRatio() {
		// build up some reads
		List<VariantData> list = new ArrayList<>();

		list.add(new VariantData(new Interval ("1", 1, 1), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HOM_VAR, new char [] {'A', 'A'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 2, 2), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HET, new char [] {'A', 'T'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 3, 3), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HOM_REF, new char [] {'A', 'T'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 4, 4), 'A', 'T', GenotypeType.HET, GenotypeType.HOM_VAR, new char [] {'T', 'T'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 5, 5), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HET, new char [] {'A', 'T'}, new int [] {10,10}));

		VariantDataCollection vcd = new VariantDataCollection(list, "one", "two", "testCell");
		double likeS1 = vcd.value(1);
		double likeS2 = vcd.value(0);
		double likePerfectMix=vcd.value(0.5);
		SamplePairAssignmentForCell s = vcd.optimizeMixture();
		double likeMixed= s.getDoubletLikelihood();
		double mixture = s.getMixture();
		// the sample mixture is slightly different, but close enough.  It seems like as I increase the number of observations (beyond just this handful) they converge.
		// this may be due to the fact that I stop convergence at 0.001 in the java version, so it can be off by a little bit.

		Assert.assertEquals(likeS1, -3.830847, 0.001);
		Assert.assertEquals(likeS2, -4.341392, 0.001);
		Assert.assertEquals(likeMixed, -3.306935, 0.001);
		Assert.assertEquals(mixture, 0.5540364, 0.015);
		Assert.assertEquals(s.getImpossibleAllelesSampleOne(),3);
		Assert.assertEquals(s.getImpossibleAllelesSampleTwo(),3);

		Assert.assertEquals(likePerfectMix, -3.313183, 0.001);
	}



	@Test(enabled=true)
	/**
	 * calculate the first 4 SNPs.  Then add in the 5th.
	 * d=data.table(data.frame(pos=c("1:1", "1:2", "1:3", "1:4"), refCount=c(2, 1,1,0), altCount=c(0,1,1,2), sampleOneGenotype=c("ref", "ref", "ref", "het"), sampleTwoGenotype=c("alt", "het", "ref", "alt"), stringsAsFactors=F))
	 * calculateLikelihoodsManyIterations(d, qualityScore=0.9, "one", "two", interval=c(0,1))
	 * sampleOneRatio jointSampleLikelihood      one       two num_paired_snps num_paired_informative_snps
	 *       0.676329             -2.586495 -2.78509 -3.739332               4                           3
	 */
	public void missingDataTest() {
		// build up some reads
		List<VariantData> list = new ArrayList<>();

		list.add(new VariantData(new Interval ("1", 1, 1), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HOM_VAR, new char [] {'A', 'A'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 2, 2), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HET, new char [] {'A', 'T'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 3, 3), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HOM_REF, new char [] {'A', 'T'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 4, 4), 'A', 'T', GenotypeType.HET, GenotypeType.HOM_VAR, new char [] {'T', 'T'}, new int [] {10,10}));
		// setting the last SNP to be missing data.
		double missingData = 0.5;
		list.add(new VariantData(new Interval ("1", 5, 5), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.NO_CALL, new char [] {'A', 'A'}, new int [] {10,10}, missingData, null, null, null));

		VariantDataCollection vcd = new VariantDataCollection(list, "one", "two", "testCell");
		double likeS1 = vcd.value(1);
		// mix=1 should be sum(log10(c(0.9, 0.9, 0.9, 0.1, 0.9, 0.1, 0.5, 0.5, 0.9, 0.9)))= -2.876605
		double likeS2 = vcd.value(0);
		// mix=0 should be sum(log10(c(0.1, 0.1, 0.5, 0.5, 0.9, 0.1, 0.9, 0.9, 0.5, 0.5)))= -4.341392
		Assert.assertEquals(likeS1, -2.876605, 0.001);
		Assert.assertEquals(likeS2, -4.341392, 0.001);

		SamplePairAssignmentForCell s = vcd.optimizeMixture();
		double likeMixed= s.getDoubletLikelihood();
		double mixture = s.getMixture();
		Assert.assertEquals(likeMixed, -2.7837, 0.001);
		Assert.assertEquals(mixture, 0.803, 0.01);

	}

	@Test(enabled=true)
	// test a true single donor.!
	public void testPerfectFirstSingleDonor() {
		// build up some reads
		List<VariantData> list = new ArrayList<>();

		list.add(new VariantData(new Interval ("1", 1, 1), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HOM_VAR, new char [] {'A', 'A'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 2, 2), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HET, new char [] {'A', 'A'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 3, 3), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HOM_REF, new char [] {'A', 'A'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 4, 4), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HOM_VAR, new char [] {'A', 'A'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 5, 5), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HET, new char [] {'A', 'A'}, new int [] {10,10}));

		VariantDataCollection vcd = new VariantDataCollection(list, "one", "two", "testCell");
		SamplePairAssignmentForCell s = vcd.optimizeMixture();
		double mixture = s.getMixture();
		Assert.assertEquals(mixture,1, 0.01);
	}

	@Test(enabled=true)
	// test a true single donor.!
	public void testPerfectSecondSingleDonor() {
		// build up some reads
		List<VariantData> list = new ArrayList<>();

		list.add(new VariantData(new Interval ("1", 1, 1), 'A', 'T', GenotypeType.HOM_VAR, GenotypeType.HOM_REF, new char [] {'A', 'A'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 2, 2), 'A', 'T', GenotypeType.HET, GenotypeType.HOM_REF, new char [] {'A', 'A'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 3, 3), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HOM_REF, new char [] {'A', 'A'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 4, 4), 'A', 'T', GenotypeType.HOM_VAR, GenotypeType.HOM_REF, new char [] {'A', 'A'}, new int [] {10,10}));
		list.add(new VariantData(new Interval ("1", 5, 5), 'A', 'T', GenotypeType.HET, GenotypeType.HOM_REF, new char [] {'A', 'A'}, new int [] {10,10}));

		VariantDataCollection vcd = new VariantDataCollection(list, "one", "two", "testCell");
		SamplePairAssignmentForCell s = vcd.optimizeMixture();
		double mixture = s.getMixture();
		Assert.assertEquals(mixture,0, 0.01);
	}


	/**
	 * These results can be replicated by running the R code.
	 * RScript -e 'source("/Users/nemesh/dropseqrna/transcriptome/R/DropSeqGenotyping/MixedSampleUnguidedLikelihoodAnalysisPrototype.R")' -e 'source("/Users/nemesh/dropseqrna/transcriptome/R/DropSeqGenotyping/simpleLikelihoodExample.R")' -e 'system.time(protoTypeFullRun(summaryLikelihoodFile="/downloads/embryonicStemCells/HUES53_NKX/TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.txt", verboseLikelihoodFile="/downloads/embryonicStemCells/HUES53_NKX/TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.txt.verbose.txt.gz", outDoubletCalls="/downloads/embryonicStemCells/HUES53_NKX/TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.txt.unguided_donors.txt", outPerSNPFile="/downloads/embryonicStemCells/HUES53_NKX/TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.unguided.perSNP.txt", outPerDonorPair="/downloads/embryonicStemCells/HUES53_NKX/TTTGCGCGGAGC:ATTGTTTAGGAG.sampleAssignments.unguided.perDonor.txt", nCores=4, topFractionToSearch=0.25, donorFile="/downloads/embryonicStemCells/HUES53_NKX/donors.txt"))'
	 * A copy of the inputs for the R code is in the test directory.
	 * cell    sampleOneMixtureRatio   sampleOne       sampleOneLikelihood     sampleTwo       sampleTwoLikelihood     mixedSample     mixedSampleLikelihood   num_paired_snps num_inform_snps lr_test_stat    sampleOneWrongAlleleCount       sampleTwoWrongAlleleCount       ratio   bestLikelihood  bestSample
	 * TTTGCGCGGAGC:ATTGTTTAGGAG       0.643   HUES53  -9.80926731658451       NKX     -12.0398723515404       HUES53:NKX      -8.72621729107288       28      18      1.0831  5       6       1.08305002551164        -8.72621729107288       HUES53:NKX
	 *
	 */
	@Test(enabled=true)
	public void bigTest() { 
		AssignCellsToSamples assigner = new AssignCellsToSamples();
		assigner.GQ_THRESHOLD=30;
		assigner.INPUT_BAM=Collections.singletonList(this.INPUT_BAM);
		assigner.VCF=this.VCF;
		assigner.NUM_BARCODES=10;
		assigner.RETAIN_MONOMORPIC_SNPS=true;

		final VCFFileReader vcfReader = new VCFFileReader(this.VCF, false);
		List<String> vcfSamples =Arrays.asList(donors);

		final SNPInfoCollection snpIntervals = SampleAssignmentVCFUtils.getSNPInfoCollection(this.VCF, vcfSamples, assigner.RETAIN_MONOMORPIC_SNPS, 30, 0, null, null, false, null);
		final PeekableIterator<VariantContext> vcfIterator = SampleAssignmentVCFUtils.getVCFIterator(vcfReader, vcfSamples, assigner.RETAIN_MONOMORPIC_SNPS, assigner.GQ_THRESHOLD, assigner.FRACTION_SAMPLES_PASSING, null, null);
		List<String> cellBarcodes = assigner.getCellBarcodes();
		PeekableIterator<List<SampleGenotypeProbabilities>> sampleGenotypeIterator = assigner.prepareIterator(snpIntervals, cellBarcodes);

		List<SampleGenotypeProbabilities> allProbs = new ArrayList<>();
		while (sampleGenotypeIterator.hasNext()) {
			List<SampleGenotypeProbabilities> r = sampleGenotypeIterator.next();
			allProbs.addAll(r);
		}


		GenotypeMatrix genotypeMatrix = new GenotypeMatrix(vcfIterator, 30, new HashSet<>(vcfSamples));
		VariantDataFactory f = new VariantDataFactory("TTTGCGCGGAGC:ATTGTTTAGGAG", allProbs, genotypeMatrix, 0.1);
		FindOptimalDonorMixture fodm = new FindOptimalDonorMixture(f);

		AllPairedSampleAssignmentsForCell result = fodm.findBestDonorPair("HUES53", vcfSamples, false);
		SamplePairAssignmentForCell best = result.getBestAssignment();
		Assert.assertEquals(best.getSampleOneSingleLikelihood(),-9.50823732092053, 0.01);
		Assert.assertEquals(best.getSampleTwoSingleLikelihood(),-11.738842355876386, 0.01);
		Assert.assertEquals(best.getDoubletLikelihood(),-8.425187288261675, 0.01);
		Assert.assertEquals(best.getMixture(), 0.643, 0.01);

	}

	@Test(enabled=true)
	public void bigTest3() {
		DetectDoublets assigner = new DetectDoublets();
		assigner.GQ_THRESHOLD=30;
		assigner.INPUT_BAM=Collections.singletonList(this.INPUT_BAM);
		assigner.VCF=this.VCF;		

		final VCFFileReader vcfReader = new VCFFileReader(this.VCF, false);
		List<String> vcfSamples =Arrays.asList(donors);
		List<String> cellBarcodes=new BarcodeListRetrieval().getListCellBarcodesByReadCount(Collections.singletonList(this.INPUT_BAM), "XC", 10, null, 100000000);
		
		final SNPInfoCollection snpIntervals = SampleAssignmentVCFUtils.getSNPInfoCollection(this.VCF, vcfSamples, true, 30, 0, null,null, false, null);
		final PeekableIterator<VariantContext> vcfIterator = SampleAssignmentVCFUtils.getVCFIterator(vcfReader, vcfSamples, true, assigner.GQ_THRESHOLD, assigner.FRACTION_SAMPLES_PASSING, null, log);
		PeekableIterator<List<SampleGenotypeProbabilities>> sampleGenotypeIterator = assigner.prepareIterator(snpIntervals.getIntervalList(), cellBarcodes);

		List<SampleGenotypeProbabilities> allProbs = sampleGenotypeIterator.next();

		GenotypeMatrix genotypeMatrix = new GenotypeMatrix(vcfIterator, assigner.GQ_THRESHOLD, vcfSamples);
		VariantDataFactory f = new VariantDataFactory("TTTGCGCGGAGC:ATTGTTTAGGAG", allProbs, genotypeMatrix, 0.1);
		FindOptimalDonorMixture fodm = new FindOptimalDonorMixture(f);

		AllPairedSampleAssignmentsForCell result = fodm.findBestDonorPair("HUES53", vcfSamples, false);
		SamplePairAssignmentForCell best = result.getBestAssignment();
		Assert.assertEquals(best.getSampleOneSingleLikelihood(),-9.80926731658451, 0.01);
		Assert.assertEquals(best.getSampleTwoSingleLikelihood(),-12.0398723515404, 0.01);
		Assert.assertEquals(best.getDoubletLikelihood(),-8.72621729107288, 0.01);
		Assert.assertEquals(best.getMixture(), 0.643, 0.01);


	}


	/**
	 * These results can be replicated by running the R code.
	 * RScript -e 'source("/Users/nemesh/dropseqrna/transcriptome/R/DropSeqGenotyping/MixedSampleUnguidedLikelihoodAnalysisPrototype.R")' -e 'source("/Users/nemesh/dropseqrna/transcriptome/R/DropSeqGenotyping/simpleLikelihoodExample.R")' -e 'system.time(protoTypeFullRun(summaryLikelihoodFile="/downloads/embryonicStemCells/HUES53_NKX/TTTGCGCGGAGC:ATTGTTTAGGAG2.sampleAssignments.txt", verboseLikelihoodFile="/downloads/embryonicStemCells/HUES53_NKX/TTTGCGCGGAGC:ATTGTTTAGGAG2.sampleAssignments.txt.verbose.txt.gz", outDoubletCalls="/downloads/embryonicStemCells/HUES53_NKX/TTTGCGCGGAGC:ATTGTTTAGGAG2.sampleAssignments.txt.unguided_donors.txt", outPerSNPFile="/downloads/embryonicStemCells/HUES53_NKX/TTTGCGCGGAGC:ATTGTTTAGGAG2.sampleAssignments.unguided.perSNP.txt", outPerDonorPair="/downloads/embryonicStemCells/HUES53_NKX/TTTGCGCGGAGC:ATTGTTTAGGAG2.sampleAssignments.unguided.perDonor.txt", nCores=4, topFractionToSearch=0.25, donorFile="/downloads/embryonicStemCells/HUES53_NKX/donors.txt"))'
	 * A copy of the inputs for the R code is in the test directory.
	 * cell    sampleOneMixtureRatio   sampleOne       sampleOneLikelihood     sampleTwo       sampleTwoLikelihood     mixedSample     mixedSampleLikelihood   num_paired_snps num_inform_snps lr_test_stat    sampleOneWrongAlleleCount       sampleTwoWrongAlleleCount       ratio   bestLikelihood  bestSample
	 * TTTGCGCGGAGC:ATTGTTTAGGAG       0.51    HUES53  -153.499300936702       	NKX     -		156.374146003812       HUES53:NKX      -125.811088971442       448     			285     		27.6882 75      74	27.6882119652593        -125.811088971442       HUES53:NKX
	 * This is a full cell's worth of data.
	 */
	@Test(enabled=true)
	public void bigTest2() {
		File vcf = new File (rootDir+"/TTTGCGCGGAGC:ATTGTTTAGGAG2.vcf");
		List<File> bam = Collections.singletonList(new File (rootDir+"/TTTGCGCGGAGC:ATTGTTTAGGAG_retagged.bam"));
		AssignCellsToSamples assigner = new AssignCellsToSamples();
		assigner.GQ_THRESHOLD=30;
		assigner.INPUT_BAM=bam;
		assigner.VCF=vcf;
		assigner.NUM_BARCODES=10;
		assigner.RETAIN_MONOMORPIC_SNPS=true;

		final VCFFileReader vcfReader = new VCFFileReader(vcf, false);
		List<String> vcfSamples =Arrays.asList(donors);

		final SNPInfoCollection snpIntervals = SampleAssignmentVCFUtils.getSNPInfoCollection(vcf, vcfSamples, assigner.RETAIN_MONOMORPIC_SNPS, 30, 0, null, null, false, null);
		final PeekableIterator<VariantContext> vcfIterator = org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils.getVCFIterator(vcfReader, vcfSamples, assigner.RETAIN_MONOMORPIC_SNPS, assigner.GQ_THRESHOLD, assigner.FRACTION_SAMPLES_PASSING, null, null);
		List<String> cellBarcodes = assigner.getCellBarcodes();
		PeekableIterator<List<SampleGenotypeProbabilities>> sampleGenotypeIterator = assigner.prepareIterator(snpIntervals, cellBarcodes);

		List<SampleGenotypeProbabilities> allProbs = new ArrayList<>();
		while (sampleGenotypeIterator.hasNext())
			allProbs.addAll(sampleGenotypeIterator.next());

		GenotypeMatrix genotypeMatrix = new GenotypeMatrix(vcfIterator, 30, vcfSamples);
		VariantDataFactory f = new VariantDataFactory("TTTGCGCGGAGC:ATTGTTTAGGAG", allProbs, genotypeMatrix, 0.1);
		FindOptimalDonorMixture fodm = new FindOptimalDonorMixture(f);

		AllPairedSampleAssignmentsForCell result = fodm.findBestDonorPair("HUES53", vcfSamples, false);
		SamplePairAssignmentForCell best = result.getBestAssignment();
		// giving the optimizer a little bit of wiggle room because it seems instable between command line, running in eclipse, etc.
		// Wiggle room needed to be increased when moving to Java 8
		Assert.assertEquals(best.getSampleOneSingleLikelihood(),-145.64069867982204, 1.5);
		Assert.assertEquals(best.getSampleTwoSingleLikelihood(),-150.10190874973378, 1);
		Assert.assertEquals(best.getDoubletLikelihood(),-119.96098316288239, 1);
		Assert.assertEquals(best.getMixture(), 0.51, 0.01);

	}
	
	
}
