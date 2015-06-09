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
