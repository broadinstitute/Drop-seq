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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import htsjdk.samtools.util.Log;

/**
 * A utility class that takes a list of strings (ordered by prevalence to determine which barcodes are merged into which)
 * and generates a map, containing the non merged barcodes as keys and the values of each key are a list of barcodes that are merged into a key.
 * @author nemesh
 *
 */
public class MapBarcodesByEditDistance {

	// private CollapseBarcodeThreaded cbt=null;
	private final int NUM_THREADS;
	private final Log log = Log.getInstance(MapBarcodesByEditDistance.class);
	private final int REPORT_PROGRESS_INTERVAL;
	// private int threadedBlockSize=20000;
	private final boolean verbose;
	private ForkJoinPool forkJoinPool;

    // https://blog.krecan.net/2014/03/18/how-to-specify-thread-pool-for-java-8-parallel-streams/
	public MapBarcodesByEditDistance (final boolean verbose, final int numThreads, final int reportProgressInterval) {
		this.verbose=verbose;
		this.NUM_THREADS=numThreads;
		this.REPORT_PROGRESS_INTERVAL=reportProgressInterval;
		// if (this.NUM_THREADS>1) cbt= new CollapseBarcodeThreaded(this.threadedBlockSize, this.NUM_THREADS);
		if (this.NUM_THREADS>1) forkJoinPool = new ForkJoinPool(numThreads);
	}

	public MapBarcodesByEditDistance (final boolean verbose, final int reportProgressInterval) {
		this(verbose, 1, reportProgressInterval);
	}

	public MapBarcodesByEditDistance (final boolean verbose) {
		this(verbose, 1, 100000);
	}

	/**
	 * Perform edit distance collapse where all barcodes are eligible to be collapse.
	 * Barcodes are ordered by total number of counts.
	 * @param barcodes
	 * @param findIndels
	 * @param editDistance
	 * @return
	 */
	public Map<String, List<String>> collapseBarcodes (final ObjectCounter<String> barcodes, final boolean findIndels, final int editDistance) {
		List<String> coreBarcodes = barcodes.getKeysOrderedByCount(true);
		return (collapseBarcodes(coreBarcodes, barcodes, findIndels, editDistance));
	}

	/**
	 * Collapses a list of barcodes.
	 * This works by iterating through every core barcode (ordered from largest to smallest in the barcodes object), and mapping any other barcode within editDistance
	 * to this barcode.  This means that the number of coreBarcodes at the end is less than or equal to the number of core barcodes submitted.  These barcodes are the keys of the output object.
	 * The barcodes that are not in the coreBarcodes list are eligible to be collapsed into a core barcode, but will never absorb other barcodes.
	 * We use coreBarcodes in order to limit the scope of the computational work, as the number of coreBarcodes can be small (the number of cells)
	 * compared to the total number of barcodes (number of beads + bead sequencing errors.)
	 * Smaller barcodes are always collapsed into larger ones.  Each barcode only exists once in the output - if a barcode A is edit distance 1 away
	 * from barcode B and C, it will be assigned to whichever of B and C is larger.
	 * @param coreBarcodes A list of barcode strings that are considered "core" or primary barcodes.
	 * @param barcodes An exhaustive list of all barcodes (both core and non-core) with assigned counts of observations of these barcodes.
	 * @param findIndels If true, we use Levenshtein indel sensitive collapse.  If false, use Hamming distance.
	 * @param editDistance The maximum edit distance two barcodes can be at to be collapsed.
	 * @return
	 */
	public Map<String, List<String>> collapseBarcodes(List<String> coreBarcodes, ObjectCounter<String> barcodes, final boolean findIndels, final int editDistance) {
		// don't allow side effects to modify input lists.
		coreBarcodes = new ArrayList<>(coreBarcodes);
		barcodes = new ObjectCounter<>(barcodes);

		Map<String, List<String>> result = new HashMap<>();
		int count = 0;
		int numBCCollapsed=0;

		List<String> barcodeList = barcodes.getKeysOrderedByCount(true);
		//List<BarcodeWithCount> barcodesWithCount=getBarcodesWithCounts(barcodes);
		//List<String> barcodeList = EDUtils.getInstance().getBarcodes(barcodesWithCount);

		//int totalCount = barcodes.getTotalCount();

		int coreBarcodeCount=coreBarcodes.size();
		long startTime = System.currentTimeMillis();
		while (coreBarcodes.isEmpty()==false) {
			String b = coreBarcodes.get(0);
			count++;
			coreBarcodes.remove(b);
			barcodeList.remove(b);

			Set<String> closeBC=processSingleBarcode(b, barcodeList, findIndels, editDistance);
			numBCCollapsed+=closeBC.size();


			if (result.containsKey(b))
				log.error("Result should never have core barcode");

			List<String> closeBCList = new ArrayList<>(closeBC);
			Collections.sort(closeBCList);
			result.put(b, closeBCList);

			barcodeList.removeAll(closeBC);
			coreBarcodes.removeAll(closeBC);
			if (this.REPORT_PROGRESS_INTERVAL!=0 && count % this.REPORT_PROGRESS_INTERVAL == 0) {
				if (barcodes.getSize()>10000) log.info("Processed [" + count + "] records, totals BC Space left [" + barcodeList.size() +"]", " # collapsed this set [" + numBCCollapsed+"]");
				numBCCollapsed=0;
			}
		}
		if (verbose) {
			long endTime = System.currentTimeMillis();
			long duration = (endTime - startTime)/1000;
			log.info("Collapse with [" + this.NUM_THREADS +"] threads took [" + duration + "] seconds to process");
			log.info("Started with core barcodes [" +coreBarcodeCount+  "] ended with [" + count + "] num collapsed [" +  (coreBarcodeCount-count) +"]");
		}
		return (result);
	}

	private Set<String> processSingleBarcode(final String barcode, final List<String> comparisonBarcodes, final boolean findIndels, final int editDistance) {
		Set<String> closeBarcodes =null;

		// Replaced with java 8 lambda method. woot?
		if (this.NUM_THREADS>1 ) {
			// closeBarcodes=cbt.getStringsWithinEditDistanceWithIndel(barcode, comparisonBarcodes, editDistance, findIndels);
			 closeBarcodes=processSingleBarcodeMultithreaded(barcode, comparisonBarcodes, findIndels, editDistance);
			 return closeBarcodes;
		}
		if (findIndels)
			closeBarcodes = EDUtils.getInstance().getStringsWithinEditDistanceWithIndel(barcode,comparisonBarcodes, editDistance);
		else
			closeBarcodes = EDUtils.getInstance().getStringsWithinEditDistance(barcode,comparisonBarcodes, editDistance);

		return (closeBarcodes);
	}

	/**
	 * The Java lambda way.
	 * @param barcode
	 * @param comparisonBarcodes
	 * @param findIndels
	 * @param editDistance
	 * @return
	 */
	private Set<String> processSingleBarcodeMultithreaded(final String barcode, final List<String> comparisonBarcodes, final boolean findIndels, final int editDistance) {
		Set<String> result = new HashSet<>();
		try {
			if (findIndels)
				result = forkJoinPool.submit(() -> comparisonBarcodes.parallelStream().filter(x -> LevenshteinDistance.getIndelSlidingWindowEditDistance(barcode, x) <= editDistance).collect(Collectors.toSet())).get();
			else
				result = forkJoinPool.submit(() -> comparisonBarcodes.stream().filter(x -> HammingDistance.getHammingDistance(barcode, x) <= editDistance).collect(Collectors.toSet())).get();
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}
		return result;
	}

	private List<BarcodeWithCount> getBarcodesWithCounts (final ObjectCounter<String> barcodes) {
		List<BarcodeWithCount> result = new ArrayList<>();
		List<String> keys = barcodes.getKeysOrderedByCount(true);
		for (String k: keys) {
			BarcodeWithCount b = new BarcodeWithCount(k, barcodes.getCountForKey(k));
			result.add(b);
		}
		return (result);
	}
}
