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

import java.util.Arrays;

/**
 * Implement the Shannon-Weaver Diversity calculations.
 * See: https://en.wikipedia.org/wiki/Diversity_index
 * @author nemesh
 *
 */
public class Diversity {

	/**
	 * Calculate the Shannon-Weaver diversity.
	 * Treat the data as proportions (each proportion data point Pi), calculate the sum of -Pi * log (Pi)
	 * In R for a vector of proportions this would be sum (sapply (proportions, function (x) -x*log(x)))
	 * @param data
	 * @return  The diversity of the data
	 */
	public static double diversity (final double [] data) {
		double [] proportions = getProportions(data);
		double diversity = Arrays.stream(proportions).map(x -> -x * Math.log(x)).sum();
		return diversity;
	}

	/**
	 * Calculate the Shannon-Weaver diversity.
	 * Treat the data as proportions (each proportion data point Pi), calculate the sum of Pi * log (Pi)
	 * @param data
	 * @return  The diversity of the data
	 */
	public static double diversity (final int [] counts) {
		double [] data = Arrays.stream(counts).asDoubleStream().toArray();
		return diversity(data);
	}

	/**
	 * Calculate the Shannon-Weaver equitability.  This is the diversity normalized to a score between 0 and 1.
	 * Treat the data as proportions (each proportion data point Pi), calculate the sum of Pi * log (Pi).  Divide by the log (number of data points).
	 * @param data
	 * @return  The diversity of the data
	 */
	public static double equitability (final double [] data) {
		double diversity = diversity(data);
		double equitability = diversity/ Math.log(data.length);
		return equitability;
	}

	/**
	 * Calculate the Shannon-Weaver equitability.  This is the diversity normalized to a score between 0 and 1.
	 * Treat the data as proportions (each proportion data point Pi), calculate the sum of Pi * log (Pi).  Divide by the log (number of data points).
	 * @param data
	 * @return  The diversity of the data
	 */
	public static double equitability (final int [] counts) {
		double [] data = Arrays.stream(counts).asDoubleStream().toArray();
		return equitability(data);
	}

	/**
	 * Convert the data to proportional representation.  Each value is divided by the sum of all values - the proportions add to one.
	 * @param data
	 * @return the data converted to proportions.
	 */
	public static double [] getProportions (final double [] data) {
		double total = Arrays.stream(data).sum();
		double [] proportions = Arrays.stream(data).map(x -> x/total).toArray();
		return (proportions);
	}




}
