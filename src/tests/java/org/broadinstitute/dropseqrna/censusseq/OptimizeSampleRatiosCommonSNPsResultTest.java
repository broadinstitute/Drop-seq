package org.broadinstitute.dropseqrna.censusseq;

import java.util.HashMap;
import java.util.Map;

import org.testng.annotations.Test;

import org.testng.Assert;

public class OptimizeSampleRatiosCommonSNPsResultTest {
	@Test
	public void test() {
		Map<String,Double> result = new HashMap<>();
		result.put("Donor1", 0.2);
		result.put("Donor2", 0.4);
		result.put("Donor3", 0.4);
		boolean converged = true;
		double bestLikelihood = 0.8;
		double secondBestLikelihood = 0.4;
		int numObservations = 4;
		boolean randomizedStart=false;
		OptimizeSampleRatiosCommonSNPsResult r = new OptimizeSampleRatiosCommonSNPsResult(result, converged, bestLikelihood, secondBestLikelihood, numObservations, randomizedStart);

		Map<String,Double> o = r.getResult();
		Assert.assertEquals(result, o);

		Assert.assertSame(r.isConverged(), converged);
		Assert.assertEquals(bestLikelihood, r.getBestLikelihood());
		Assert.assertEquals(secondBestLikelihood, r.getSecondBestLikelihood());
		Assert.assertEquals(numObservations, r.getNumObservations());
		Assert.assertSame(randomizedStart, r.isRandomizedStart());
		Assert.assertEquals(bestLikelihood/numObservations, r.getNormalizedLikelihood(), 0.001);
		Assert.assertEquals(bestLikelihood-secondBestLikelihood, r.getLikelihoodDelta(), 0.001);

	}
}
