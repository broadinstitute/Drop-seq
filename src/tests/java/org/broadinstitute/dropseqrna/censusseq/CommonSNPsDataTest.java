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
package org.broadinstitute.dropseqrna.censusseq;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.variant.variantcontext.GenotypeType;

public class CommonSNPsDataTest {

	/*
	 * $readCounts
  	refCount altCount
        9        1
        8        2
        6        4

	$genotypeStates
     		sample1 sample2 sample3 sample4
	[1,]     0.5     0.0     0.0     0.0
	[2,]     0.5     0.5     0.0     0.5
	[3,]     0.5     1.0     0.5     0.0

	sampleMixture 0.25 0.25 0.25 0.25
	minorAllelefreqs 0.125 0.375 0.500
 	likeStart -6.920215

	sampleMixture 0.1 0.2 0.3 0.4
	minorAllelefreqs 0.05 0.35 0.40
	likeEnd -6.832927

	*/
	@Test
	public void testAddData () {
		String [] sampleNames = {"sample1", "sample2", "sample3", "sample4"};
		int [] [] gs = {{1,0,0,0}, {1,1,0,1},{1,2,1,0}};
		int [] [] alleleCounts={{9,1},{8,2},{6,4}};
		// [{2, 4, 6}, {4, 7}, {6}, {4}]
		CommonSNPsData d = new CommonSNPsData(Arrays.asList(sampleNames));

		for (int i=0; i<gs.length; i++) {
			int [] genos=gs[i];
			int [] ac = alleleCounts[i];
			d.addSNP(sampleNames, genos, ac[0], ac[1]);
		}
		int numVars = d.getNumVariants();
		Assert.assertEquals(gs.length, numVars);

		for (int i=0; i<numVars; i++) {
			int [] gsNew = d.getCountsAltAllele(i);
			Assert.assertEquals(gs[i], gsNew);

		}
	}

	@Test
	public void testGetAllelesFreqs () {
		double [] sampleMixtureStart = {0.25, 0.25, 0.25, 0.25};
		double [] sampleMixtureEnd = {0.1, 0.2, 0.3, 0.4};

		String [] sampleNames = {"sample1", "sample2", "sample3", "sample4"};
		int [] [] gs = {{1,0,0,0}, {1,1,0,1},{1,2,1,0}};
		int [] [] alleleCounts={{9,1},{8,2},{6,4}};
		// [{2, 4, 6}, {4, 7}, {6}, {4}]
		CommonSNPsData d = new CommonSNPsData(Arrays.asList(sampleNames));

		for (int i=0; i<gs.length; i++) {
			int [] genos=gs[i];
			int [] ac = alleleCounts[i];
			d.addSNP(sampleNames, genos, ac[0], ac[1]);
		}

		int totalNumReads=30;
		int tnr = d.getTotalSNPAlleleCounts();
		Assert.assertEquals(tnr, totalNumReads);



		double [] resultStart = d.getWeighedAlleleFrequencies(sampleMixtureStart);
		double [] resultEnd = d.getWeighedAlleleFrequencies(sampleMixtureEnd);

		double [] expectedStart = {0.125, 0.375, 0.500};
		double [] expectedEnd = {0.05, 0.35, 0.40};

		for (int i=0; i<resultStart.length; i++) {
			Assert.assertEquals(resultStart[i], expectedStart[i], 0.001);
			Assert.assertEquals(resultEnd[i], expectedEnd[i], 0.001);
		}
	}

	@Test
	public void testGetAllelesFreqsMissingData () {
		double [] sampleMixtureStart = {0.25, 0.25, 0.25, 0.25};
		double [] sampleMixtureEnd = {0.1, 0.2, 0.3, 0.4};

		String [] sampleNames = {"sample1", "sample2", "sample3", "sample4"};
		int [] [] gs = {{1,0,0,-1}, {1,1,0,1},{1,2,1,0}};
		int [] [] alleleCounts={{9,1},{8,2},{6,4}};
		// [{2, 4, 6}, {4, 7}, {6}, {4}]
		CommonSNPsData d = new CommonSNPsData(Arrays.asList(sampleNames));
		for (int i=0; i<gs.length; i++) {
			int [] genos=gs[i];
			int [] ac = alleleCounts[i];
			d.addSNP(sampleNames, genos, ac[0], ac[1]);
		}


		double [] resultStart = d.getWeighedAlleleFrequencies(sampleMixtureStart);
		double [] resultEnd = d.getWeighedAlleleFrequencies(sampleMixtureEnd);

		double [] expectedStart = {0.250, 0.375, 0.5};
		// example for the second result: weighted.mean(x=c(0.5, 0.5,0, 0.5), w=c(0.1, 0.2, 0.3, 0.4)); //  0.35
		// weighted.mean(x=c(0.5, 1,0.5, 0), w=c(0.1, 0.2, 0.3, 0.4));
		double [] expectedEnd = {0.25, 0.35, 0.4};

		for (int i=0; i<resultStart.length; i++) {
			Assert.assertEquals(resultStart[i], expectedStart[i], 0.001);
			Assert.assertEquals(resultEnd[i], expectedEnd[i], 0.001);
		}
	}

	@Test
	public void testGetAllelesFreqsMissingData2 () {
		double [] sampleMixture = {0.25, 0.25, 0.25, 0.25};


		String [] sampleNames = {"sample1", "sample2", "sample3", "sample4"};
		int [] [] gs = {{1,0,0,-1}, {1,1,0,1},{1,2,1,0}};
		int [] [] alleleCounts={{9,1},{8,2},{6,4}};
		// [{2, 4, 6}, {4, 7}, {6}, {4}]
		CommonSNPsData d = new CommonSNPsData(Arrays.asList(sampleNames));
		for (int i=0; i<gs.length; i++) {
			int [] genos=gs[i];
			int [] ac = alleleCounts[i];
			d.addSNP(sampleNames, genos, ac[0], ac[1]);
		}

		double missingAsHet = d.getWeighedAlleleFrequenciesOneSNPMissingAsHet(sampleMixture, 0);
		double missingRemoved = d.getWeighedAlleleFrequenciesOneSNPMissingRemoved(sampleMixture, 0);

		double expectedMissingAsHet = 0.25;
		double expectedMissingRemoved = 0.166;
		Assert.assertEquals(missingAsHet, expectedMissingAsHet, 0.001);
		Assert.assertEquals(missingRemoved, expectedMissingRemoved, 0.001);

	}




	@Test
	public void testPrintAlleleCountsPerDonor () {
		File snpReadCountFile = new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalallelecounts/sampleassignment/genomicpool/commonsnps/readCounts.txt.gz");
		File genotypeFile = new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalallelecounts/sampleassignment/genomicpool/commonsnps/sampleGenotypeStates.txt.gz");
		CommonSNPsData result = CommonSNPsData.parseFromFiles(snpReadCountFile, genotypeFile);
		Assert.assertTrue(result!=null);
		Assert.assertEquals(result.getNumVariants(), 50000);
		Assert.assertEquals(result.getSampleNames().size(), 20);

		List<String> samples = result.getSampleNames();
		for (String s: samples) {
			ObjectCounter<GenotypeType> g= result.getCountGenotypes(s);
			System.out.println(g);
		}

	}

	@Test
	public void testFileParser() {
		File snpReadCountFile = new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalallelecounts/sampleassignment/genomicpool/commonsnps/readCounts.txt.gz");
		File genotypeFile = new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalallelecounts/sampleassignment/genomicpool/commonsnps/sampleGenotypeStates.txt.gz");
		CommonSNPsData result = CommonSNPsData.parseFromFiles(snpReadCountFile, genotypeFile);
		Assert.assertTrue(result!=null);
		Assert.assertEquals(result.getNumVariants(), 50000);
		Assert.assertEquals(result.getSampleNames().size(), 20);
	}

	@Test
	public void testFileParserEmptyFile() throws IOException {
		File outFile = File.createTempFile("commonSNPsTest.", ".empty.txt");
		outFile.deleteOnExit();
		CommonSNPsData result = CommonSNPsData.parseFromFiles(outFile, outFile);
		Assert.assertNull(result);
	}

}
