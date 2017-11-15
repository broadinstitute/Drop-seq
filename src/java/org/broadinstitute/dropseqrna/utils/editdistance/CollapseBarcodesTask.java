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
package org.broadinstitute.dropseqrna.utils.editdistance;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.RecursiveTask;

public class CollapseBarcodesTask extends RecursiveTask<Set<String>> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1358460693439743534L;
	
	private String baseString;
	private List<String> comparisonStrings;
	private int splitSize=1000;
	private int editDistance;
	private int threshold;
	private boolean findIndels;
	
	/**
	 * Find <comparisonStrings> that are within <editDistance> of the <baseString>.
	 * @param baseString  The main string to find neighbors of
	 * @param comparisonStrings The List of strings to compare to the baseString
	 * @param splitSize How many strings should be compared in a block.  This controls how to split the data for the Fork / Join operation
	 * @param editDistance What is the maximum edit distance between the <baseString> and a compared string for the compared string to be returned  
	 * If the naive edit distance is greater than threshold, then <threshold> is returned as the edit distance.  This is saying "Don't waste time computing edit distances
	 * in the more costly way when the naive edit distance is already high." 
	 * @param findIndels Should indels with the sliding window technique be detected?
	 */
	public CollapseBarcodesTask(String baseString, List<String> comparisonStrings, int splitSize, int editDistance, boolean findIndels) {
		this.baseString=baseString;
		this.comparisonStrings=comparisonStrings;
		this.splitSize=splitSize;
		this.editDistance=editDistance;
		this.findIndels=findIndels;
	}
	
	@Override
	protected Set<String> compute() {
		
		int size=this.comparisonStrings.size();
		// check to see if the data should be split
		if (size>this.splitSize) {
			// split data in half.
			int midpoint=size/2;
			List<String> n1 = comparisonStrings.subList(0, midpoint);
			List<String> n2 = comparisonStrings.subList(midpoint, size-1);
			CollapseBarcodesTask t1 = new CollapseBarcodesTask(baseString, n1, this.splitSize, this.editDistance, this.findIndels);
			CollapseBarcodesTask t2 = new CollapseBarcodesTask(baseString, n2, this.splitSize, this.editDistance, this.findIndels);
			t1.fork();
            Set<String> result = t2.compute();
            result.addAll(t1.join());
            return (result);
            
		} 
		// compute task
		Set<String> result = new HashSet<String>();
		for (String b : this.comparisonStrings) {
			int ed;
			if (this.findIndels) {
				// EDUtils.getInstance().getStringsWithinEditDistanceWithIndel(b, comparisonStrings)
				ed = LevenshteinDistance.getIndelSlidingWindowEditDistance(this.baseString, b, threshold);
			} else {
				ed = HammingDistance.getHammingDistance(this.baseString, b);
			}
			
			if (ed <= editDistance)
				result.add(b);
		}
		return result;
	}
	
	
}
