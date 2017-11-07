/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils.editdistance;

import java.util.List;
import java.util.Set;
import java.util.concurrent.ForkJoinPool;

public class CollapseBarcodeThreaded {

	private ForkJoinPool threadpool;
	private final int blockSize;
	
	public CollapseBarcodeThreaded (int blockSize, Integer numThreads) {
		if (numThreads==null) {
			threadpool = new ForkJoinPool();
		} else {
			threadpool= new ForkJoinPool(numThreads);
		}
		this.blockSize = blockSize;
	}
	
	
	public Set<String> getStringsWithinEditDistanceWithIndel(String baseString,
			List<String> comparisonStrings, int editDistance, boolean findIndels) {
		return threadpool.invoke(new CollapseBarcodesTask(baseString, comparisonStrings, blockSize, editDistance, findIndels));
	}
	
	public int getNumThreads () {
		return this.threadpool.getParallelism();
	}

		
	
}
