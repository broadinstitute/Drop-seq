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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.broadinstitute.dropseqrna.censusseq.CommonSNPsData;
import org.junit.Assert;
import org.testng.annotations.Test;

import picard.util.BasicInputParser;

public class OptimizeSampleRatiosCommonSNPsTest {


	/*
	 * Going to need real data to test this.
	 * Generate 50000 variants and 10 donors worth of data.
	 * Write file parsers for genotype states and observed alleles
	 *
	*/

	// Changes to direct iteration to include an optimization step make this test fail.
	@Test(enabled=true)
	public void testOptimize () {
		List<double []> answerKey = getAnswerKey();
		File snpReadCountFile = new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalallelecounts/sampleassignment/genomicpool/commonsnps/readCounts.txt.gz");
		File genotypeFile = new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalallelecounts/sampleassignment/genomicpool/commonsnps/sampleGenotypeStates.txt.gz");
		CommonSNPsData d = CommonSNPsData.parseFromFiles(snpReadCountFile, genotypeFile);
		// data was calculated with an iteration factor of 0.05, need to maintain that.
		// changes to optimized iteration factor ruin this test.
		OptimizeSampleRatiosCommonSNPs optim = new OptimizeSampleRatiosCommonSNPs(d, 1, null);
		Map<String, Double> ratioMap = optim.directIteration().getResult();

		double [] answerDirectIteration=answerKey.get(2);
		List<String> sampleNames = d.getSampleNames();
		for (int i=0; i<answerDirectIteration.length; i++) {
			double actual = ratioMap.get(sampleNames.get(i));
			double expected = answerDirectIteration[i];
			Assert.assertEquals(expected, actual, 0.001);
		}

	}

	@Test
	// tests starting points random/non random, and logging
	public void testGetStartPoint () {
		File snpReadCountFile = new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalallelecounts/sampleassignment/genomicpool/commonsnps/readCounts.txt.gz");
		File genotypeFile = new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalallelecounts/sampleassignment/genomicpool/commonsnps/sampleGenotypeStates.txt.gz");
		CommonSNPsData d = CommonSNPsData.parseFromFiles(snpReadCountFile, genotypeFile);
		Random random = new Random(1);
		OptimizeSampleRatiosCommonSNPs optim = new OptimizeSampleRatiosCommonSNPs(d, 1, random);
		Map<String, Double> ratioMap = optim.directIteration().getResult();

		// values of the donors after optimization.
		double [] vals1 = ratioMap.values().stream().mapToDouble(x->x.doubleValue()).toArray();

		int numSamples=d.getSampleNames().size();
		// test the starting point for equal ratios.  1/ num samples.
		double [] starting = optim.getStartPoint(numSamples, 1, false);
		double expectedStart = (double) 1 / (double) starting.length;
		for (double element2 : starting)
			Assert.assertEquals(element2, expectedStart, 0.001);

		// test that the value log doesn't throw any large null-pointer type errors.
		List<double []> mixtureResults = new ArrayList<>();
		mixtureResults.add(starting);
		mixtureResults.add(vals1);
		optim.logGradientError(mixtureResults);

		// test the starting point for unequal ratios.  The values will not all be 1/n.
		starting = optim.getStartPoint(numSamples, 1, true);
		boolean allEqual=true;
		for (double element : starting)
			// if a result is not 1/n,
			if (Math.abs(element - (double) 1/starting.length) > 0.0001)
				allEqual=false;
		Assert.assertFalse(allEqual);
	}


	private List<double []>  getAnswerKey () {
		File answerKeyFile = new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalallelecounts/sampleassignment/genomicpool/commonsnps/answer_key.txt");
		BasicInputParser snpParser = new BasicInputParser(false, answerKeyFile);
		List<double []> result = new ArrayList<>();
		snpParser.next(); // skip header.
		while (snpParser.hasNext()){
			double [] mixtures = parseSNPReadCountLine(snpParser);
			result.add(mixtures);
		}
		return result;
	}

	private static double [] parseSNPReadCountLine (final BasicInputParser parser) {
		if (parser.hasNext()) {
			String [] line =parser.next();
			double [] l = new double [line.length-1];
			for (int i=0; i<l.length; i++)
				l[i]=Double.parseDouble(line[i+1]);
			return l;
		}
		return null;
	}
}
