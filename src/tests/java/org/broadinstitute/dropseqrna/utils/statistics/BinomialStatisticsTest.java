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
package org.broadinstitute.dropseqrna.utils.statistics;

import org.apache.commons.math3.stat.interval.ConfidenceInterval;
import org.apache.commons.math3.stat.interval.IntervalUtils;
import org.broadinstitute.dropseqrna.utils.statistics.BinomialStatistics;
import org.testng.Assert;
import org.testng.annotations.Test;

public class BinomialStatisticsTest {


	@Test(enabled=true)
	// Some of the methods (Clopper-Pearson method, Normal approximation) throws an exception because there are 0 trials.
	public void testData () {
		ConfidenceInterval i= IntervalUtils.getWilsonScoreInterval(5, 0, 0.95);
		ConfidenceInterval i2= IntervalUtils.getWilsonScoreInterval(5, 5, 0.95);
		Assert.assertNotNull(i);
		Assert.assertNotNull(2);
	}

	@Test(enabled=false)
	// this throws an exception because there are 0 trials.
	public void testNoData () {
		BinomialStatistics s1 = new BinomialStatistics(0, 0, 0.95, 0.5);
		Assert.assertNotNull(s1);
	}


	@Test(enabled=true)
	public void getBinomialConfidenceInterval() {
		//binom.test(x=10, n=20, p=0.5, alternative="two.sided", conf.level = 0.95)$conf.int[1:2]
		BinomialStatistics s1 = new BinomialStatistics(20, 10, 0.95, 0.5);
		Assert.assertEquals(s1.getBinomialConfidenceInterval().getLowerBound(), 0.299, 0.001);
		Assert.assertEquals(s1.getBinomialConfidenceInterval().getUpperBound(), 0.700, 0.001);

		//binom.test(x=15, n=20, p=0.5, alternative="two.sided", conf.level = 0.95)$conf.int[1:2]
		BinomialStatistics s2 = new BinomialStatistics(20, 15, 0.95, 0.5);
		Assert.assertEquals(s2.getBinomialConfidenceInterval().getLowerBound(), 0.53129, 0.001);
		Assert.assertEquals(s2.getBinomialConfidenceInterval().getUpperBound(), 0.88813, 0.001);

		//binom.test(x=5, n=20, p=0.5, alternative="two.sided", conf.level = 0.95)$conf.int[1:2]
		BinomialStatistics s3 = new BinomialStatistics(20, 5, 0.95, 0.5);
		Assert.assertEquals(s3.getBinomialConfidenceInterval().getLowerBound(), 0.11186, 0.001);
		Assert.assertEquals(s3.getBinomialConfidenceInterval().getUpperBound(), 0.46870, 0.001);
	}

	// R code to obtain pvalues:
	//	twoTailedUnbiasedTest<-function (numSuccesses, numTrials) {
	//	 	x=c(numSuccesses, numTrials-numSuccesses)
	//	 	dbinom(min(x), sum(x), 0.5) + ifelse( min(x)==0, 0, 2*sum(dbinom(0:(min(x)-1), sum(x), 0.5)))
	// 	}
	@Test(enabled=true)
	public void getTwoSidedPvalue() {

		// binom.test(x=10, n=20, p=0.5, alternative="two.sided", conf.level = 0.95)$p.value
		double pval = BinomialStatistics.getTwoSidedPvalue(20, 10, 0.5);
		Assert.assertEquals(pval, 1, 0.001);

		// binom.test(x=15, n=20, p=0.5, alternative="two.sided", conf.level = 0.95)$p.value
		pval = BinomialStatistics.getTwoSidedPvalue(20, 15, 0.5);
		Assert.assertEquals(pval, 0.0266037, 0.001);

		// binom.test(x=5, n=20, p=0.5, alternative="two.sided", conf.level = 0.95)$p.value
		pval = BinomialStatistics.getTwoSidedPvalue(20, 15, 0.5);
		Assert.assertEquals(pval, 0.0266037, 0.001);
	}


	@Test(enabled=true)
	public void getBinomialPvalueTruncated() {

		// binom.test(x=10, n=20, p=0.5, alternative="two.sided", conf.level = 0.95)$p.value
		double pval = BinomialStatistics.getTwoSidedPvalueApacheTruncated(20, 10, 0.5);
		Assert.assertEquals(pval, 1, 0.001);

		// binom.test(x=15, n=20, p=0.5, alternative="two.sided", conf.level = 0.95)$p.value
		pval = BinomialStatistics.getTwoSidedPvalueApacheTruncated(20, 15, 0.5);
		Assert.assertEquals(pval, 0.04138947, 0.001);

		// binom.test(x=5, n=20, p=0.5, alternative="two.sided", conf.level = 0.95)$p.value
		pval = BinomialStatistics.getTwoSidedPvalueApacheTruncated(20, 15, 0.5);
		Assert.assertEquals(pval, 0.04138947, 0.001);

		// binom.test(x=1, n=20, p=0.5, alternative="two.sided", conf.level = 0.95)$p.value
		pval = BinomialStatistics.getTwoSidedPvalueApacheTruncated(20, 1, 0.5);
		Assert.assertEquals(pval, 0.00004005432, 0.000001);
	}

	@Test(enabled=true)
	public void getRatio() {
		BinomialStatistics s1 = new BinomialStatistics(20, 10, 0.95, 0.5);
		Assert.assertEquals(s1.getRatio(), 0.5, 0.001);

		BinomialStatistics s2 = new BinomialStatistics(30, 10, 0.95, 0.5);
		Assert.assertEquals(s2.getRatio(), 0.333, 0.001);

	}

	@Test(enabled=true)
	public void other () {
		BinomialStatistics s1 = new BinomialStatistics(7, 7, 0.95, 0.5);
		double pval = s1.getBinomialPvalue();
		ConfidenceInterval i = s1.getBinomialConfidenceInterval();

	}
}
