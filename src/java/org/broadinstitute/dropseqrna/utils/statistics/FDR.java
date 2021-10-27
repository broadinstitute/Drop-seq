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

import java.util.stream.IntStream;

public class FDR {

	/**
	 * Benjamini-Hochberg pvalue adjustment.
	 * See http://www.biostathandbook.com/multiplecomparisons.html section on BH pvalue adjustment.
	 * @param pvalues A list of pvalues to adjust.
	 * @return
	 */
	public static double [] bhAdjustment (final double [] pvalues) {
		if (pvalues==null ||  pvalues.length==0) return pvalues;

		//sorts from smallest to largest pvalue.
		int [] idx = getSortOrder(pvalues);
		// rank of the pvalue is the count from the end of the array back, +1 (array is 0 based.)
		int currentRank=idx.length;

		double [] fdrCorrectedPvalues= new double [pvalues.length];
		// the last value adjusted is the same as the last pvalue.
		fdrCorrectedPvalues[idx[currentRank-1]]=pvalues[idx[currentRank-1]];
		currentRank--;
		int n=pvalues.length;

		// begin to process.
		while (currentRank>0) {
			double previousCorrectedPvalue=fdrCorrectedPvalues[idx[currentRank]];
			double uncorrectedPvalue=pvalues[idx[currentRank-1]];
			double correctedPvalue =  uncorrectedPvalue*((double)n/(double)currentRank);
			double bestValue= Math.min(previousCorrectedPvalue,correctedPvalue);
			fdrCorrectedPvalues[idx[currentRank-1]]=bestValue;
			currentRank--;
		}
		return fdrCorrectedPvalues;
	}

	/**
	 * Given a list of doubles, provide the index that would sort these numbers.
	 * IE: if you iterate through the array in the index rank ordering, you'd iterate from the smallest to largest value.
	 * @param values a bunch of doubles.
	 * @return a list of rankings (index) of the values from the smallest to largest value.
	 */
	public static int [] getSortOrder (final double [] values) {
		int[] sortedIndices = IntStream.range(0, values.length)
		                .boxed().sorted((i, j) -> Double.compare(values[i], values[j]))
		                .mapToInt(ele -> ele).toArray();
		return sortedIndices;
	}
}
