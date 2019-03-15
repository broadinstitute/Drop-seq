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
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;

import org.broadinstitute.dropseqrna.metrics.UmiSharingMetrics;
import org.broadinstitute.dropseqrna.metrics.umisharing.ParentEditDistanceMatcher;
import org.broadinstitute.dropseqrna.metrics.umisharing.ParentEditDistanceMatcher.TagValues;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import htsjdk.samtools.util.Log;

/**
 * A utility class that takes a list of strings (ordered by prevalence to determine which barcodes are merged into which)
 * and generates a map, containing the non merged barcodes as keys and the values of each key are a list of barcodes that are merged into a key.
 * @author nemesh
 *
 */
public class MapBarcodesByEditDistance {

	private final int NUM_THREADS;
	private final Log log = Log.getInstance(MapBarcodesByEditDistance.class);
	private final int REPORT_PROGRESS_INTERVAL;
	private final boolean verbose;
	private ForkJoinPool forkJoinPool;

	public MapBarcodesByEditDistance (final boolean verbose, final int numThreads, final int reportProgressInterval) {
		this.verbose=verbose;
		this.NUM_THREADS=numThreads;
		this.REPORT_PROGRESS_INTERVAL=reportProgressInterval;
		// https://blog.krecan.net/2014/03/18/how-to-specify-thread-pool-for-java-8-parallel-streams/
		forkJoinPool = new ForkJoinPool(numThreads);
	}

	public MapBarcodesByEditDistance (final boolean verbose, final int reportProgressInterval) {
		this(verbose, 1, reportProgressInterval);
	}

	public MapBarcodesByEditDistance (final boolean verbose) {
		this(verbose, 1, 100000);
	}

	/**
	 * Find if a repaired barcode [which have an N at the last base] can be matched up to one and ONLY one other barcode at the edit distance.
	 * This match up can be one of two ways:
	 * 1) An indel relationship at any position
	 * 2) A substitution relationship at base 12.
	 *
	 * @param repairedCellBarcodes cell barcodes that have been repaired, and we need to find their intended sequence.  They should all end in N.
	 * @param potentialIntendedSequences Other cell barcodes that do not exhibit UMI Bias, and are the search set to find the "original" barcodes that gave rise to the repaired barcodes.
	 * @return A map from the repaired sequence to the intended sequence.
	 */
	public Map<String, String> findIntendedIndelSequences (final Collection<String> repairedCellBarcodes, final List<String> potentialIntendedSequences, final int editDistance) {
		Map<String, String> result=new HashMap<>();
		long startTime = System.currentTimeMillis();

		for (String repairedBC: repairedCellBarcodes) {
			Set<String> possibleIntendedSequences = findIntendedIndelSequences(repairedBC, potentialIntendedSequences, editDistance);
			if (possibleIntendedSequences.size()==1)
				result.put(repairedBC, possibleIntendedSequences.iterator().next());
		}

		if (verbose) {
			long endTime = System.currentTimeMillis();
			long duration = (endTime - startTime)/1000;
			log.info("Collapsed ["+repairedCellBarcodes.size()+"] barcodes against [" + potentialIntendedSequences.size() + "] targets with [" + this.NUM_THREADS +"] threads took [" + duration + "] seconds to process");
		}
		return result;
	}

	/**
	 * An intended sequence is one of two things:
	 * 1) An indel relationship (not a substitution) at base 1-11
	 * 2) A substitution at base 12.
	 * At ED=1, because the repairedCellBarcode has an N at the end, it's automatically substitution >1 if the base isn't 12.
	 * @param repairedCellBarcode
	 * @param otherBarcodes
	 * @param editDistance
	 * @return
	 */
	Set<String> findIntendedIndelSequences (final String repairedCellBarcode, final List<String> otherBarcodes, final int editDistance) {
		Set<String> indelResult = processSingleBarcodeMultithreaded(repairedCellBarcode, otherBarcodes, true, editDistance);
		Set<String> hammingResult = processSingleBarcodeMultithreaded(repairedCellBarcode, otherBarcodes, false, editDistance);
		Set<String> hammingResultPos12 = new HashSet<>();

		// indel changes must be deletions in the repaired region (or insertions in the intended sequence)
		Set<String> indelResultFiltered = new HashSet<>();
		for (String s: indelResult) {
			LevenshteinDistanceResult r= LevenshteinDistance.computeLevenshteinDistanceResult(s, repairedCellBarcode, 1, 1, 2);
			String [] ops  = r.getOperations();
			// any position before the last is D, and last is I.
			for (int i=0; i<ops.length-2; i++)
				if (ops[i].equals("D") && ops[ops.length-1].equals("I")) {
					indelResultFiltered.add(s);
					break;
				}

		}

		// hamming result must be at the last position.
		for (String s: hammingResult) {
			int [] pos = HammingDistance.getHammingDistanceChangePositions(repairedCellBarcode, s);
			// if there's one change at the last base...
			if (pos.length==1 && pos[0]==s.length()-1)
				hammingResultPos12.add(s);
		}

		// add any hamming position 12 results to the indel results.
		indelResultFiltered.addAll(hammingResultPos12);
		return indelResultFiltered;

	}

	

	public BottomUpCollapseResult bottomUpCollapse (final ObjectCounter<String> barcodes, final int editDistance) {

		BottomUpCollapseResult result = new BottomUpCollapseResult();
		// if there are no barcodes to collapse, return an empty collapse result.
		if (barcodes.getSize()==0)
			return result;
		// ordered from smallest to largest.
		List<String> barcodeList = barcodes.getKeysOrderedByCount(false);
		List<char []> barcodeListArrays = barcodeList.stream().map(x-> x.toCharArray()).collect(Collectors.toList());

		// assert all the char [] are the same length to speed things up.
		if (!assertAllArraysSameLength(barcodeListArrays))
			throw new IllegalArgumentException("This collapse requires all strings to be the same length!");

		long startTime = System.currentTimeMillis();

		// process [i] vs [i+1:(end-1)]
		// can't collapse the last barcode with nothing...
		int len=barcodeListArrays.size();
		for (int i=0; i<(len-1); i++) {
			String smallBC = barcodeList.get(i);
			List<char [] > largerBarcodes= barcodeListArrays.subList(i+1, len);
			// get the small barcode as the char []
			Set<String> largerRelatedBarcodes = processHammingDistanceEqualSizedStrings(barcodeListArrays.get(i), largerBarcodes, editDistance);

			// if there's just 1 larger neighbor, the result is unambiguous.
			if (largerRelatedBarcodes.size()==1 ) {
				String largerBarcode=largerRelatedBarcodes.iterator().next();
				// only add if the barcode is larger.
				// this avoids ordering issues with equally sized barcodes vs alphanumeric barcode sorting.
				if (barcodes.getCountForKey(smallBC)<= barcodes.getCountForKey(largerBarcode))
					result.addPair(smallBC, largerBarcode);
			}
			if (largerRelatedBarcodes.size()>1)
				result.addAmbiguousBarcode(smallBC);
			if (this.REPORT_PROGRESS_INTERVAL!=0 && i % this.REPORT_PROGRESS_INTERVAL == 0)
				log.info("Processed [" + i + "] records of [" +len+"] barcodes");
		}

		if (verbose) {
			long endTime = System.currentTimeMillis();
			long duration = (endTime - startTime)/1000;
			log.info("Collapsed ["+len+"] barcodes with [" + this.NUM_THREADS +"] threads took [" + duration + "] seconds to process");
		}
		return result;
	}


	/**
	 * This is a particular optimization for large numbers of strings to avoid converting Strings into char []'s repeatedly.
	 *
	 * @param barcode
	 * @param comparisonBarcodes
	 * @param editDistance
	 * @return
	 */
	private Set<String> processHammingDistanceEqualSizedStrings (final char [] barcode, final List<char []> comparisonBarcodes, final int editDistance) {
		Set<char []> resultTemp = Collections.EMPTY_SET;
		if (this.NUM_THREADS==1)
			resultTemp = comparisonBarcodes.stream().filter(x -> HammingDistance.getHammingDistanceEqualSizedStrings(barcode, x) <= editDistance).collect(Collectors.toSet());
		else
			try {
				resultTemp = forkJoinPool.submit(() -> comparisonBarcodes.parallelStream().filter(x -> HammingDistance.getHammingDistanceEqualSizedStrings(barcode, x) <= editDistance).collect(Collectors.toSet())).get();
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		// convert the char [] back to strings
		Set<String> result = resultTemp.stream().map(x -> String.valueOf(x)).collect(Collectors.toSet());
		return result;
	}

	/**
	 * Make sure all character arrays are the same size so you can use the short-cut hamming distance that doesn't do this check.
	 * @param charArrays
	 * @return
	 */
	private boolean assertAllArraysSameLength (final Collection<char []> charArrays) {
		int length = charArrays.iterator().next().length;
		for (char [] s: charArrays){
			int l = s.length;
			if (l!=length) return false;
		}
		return true;
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
	 * Collapses barcodes by the given edit distance and indel settings.
	 * Returns an ObjectCounter of the barcodes, with the counts of the barcodes updated to reflect barcodes that were collapsed.
	 * @param barcodes
	 * @param findIndels
	 * @param editDistance
	 * @return
	 */
	public ObjectCounter<String> collapseAndMergeBarcodes (final ObjectCounter<String> barcodes, final boolean findIndels, final int editDistance) {
		if (editDistance==0) return barcodes;
		ObjectCounter <String> result = new ObjectCounter<>();
		Map<String, List<String>> collapseMap = this.collapseBarcodes(barcodes, findIndels, editDistance);

		for (String key: collapseMap.keySet()) {
			int totalCount = barcodes.getCountForKey(key);
			List<String> values = collapseMap.get(key);
			for (String bc: values) {
				int count = barcodes.getCountForKey(bc);
				totalCount+=count;
			}
			result.setCount(key, totalCount);
		}
		return (result);

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
	public Map<String, List<String>> collapseBarcodes (List<String> coreBarcodes, ObjectCounter<String> barcodes, final boolean findIndels, final int editDistance) {
		FindSimilarEntitiesByEditDistance function = new FindSimilarEntitiesByEditDistance(this, findIndels, editDistance);
		FindSimilarEntitiesResult<String, String> result = collapseBarcodesGeneric(coreBarcodes, barcodes, function, 0);
		return result.getEntityMap();
	}

	/**
	 * Collapses a list of barcodes.
	 * This works by iterating through every core barcode (ordered from largest to smallest in the barcodes object), and mapping any other barcode close
	 * to this barcode.
	 *
	 * In the case of adaptive collapse, "close" uses the edit distance distribution of all barcodes to the core barcode being examined, and assumes that this
	 * distribution is bimodal.  The edit distance of lowest density is selected as the maximal edit distance.
	 *
	 * This means that the number of coreBarcodes at the end is less than or equal to the number of core barcodes submitted.  These barcodes are the keys of the output object.
	 * The barcodes that are not in the coreBarcodes list are eligible to be collapsed into a core barcode, but will never absorb other barcodes.
	 * We use coreBarcodes in order to limit the scope of the computational work, as the number of coreBarcodes can be small (the number of cells)
	 * compared to the total number of barcodes (number of beads + bead sequencing errors.)
	 * Smaller barcodes are always collapsed into larger ones.  Each barcode only exists once in the output - if a barcode A is edit distance 1 away
	 * from barcode B and C, it will be assigned to whichever of B and C is larger.
	 * @param coreBarcodes A list of barcode strings that are considered "core" or primary barcodes.
	 * @param barcodes An exhaustive list of all barcodes (both core and non-core) with assigned counts of observations of these barcodes.
	 * @param findIndels If true, we use Levenshtein indel sensitive collapse.  If false, use Hamming distance.
	 * @param defaultEditDistance If the discovered edit distance threshold is less than this number, this is used instead.  Set to 0 to effectively ignore the parameter.
	 * @return
	 */
	public AdaptiveMappingResult collapseBarcodesAdaptive(List<String> coreBarcodes, ObjectCounter<String> barcodes, final boolean findIndels, final int defaultEditDistance, final int minEditDistance, final int maxEditDistance) {
		FindSimilarEntitiesByAdaptiveEditDistance function = new FindSimilarEntitiesByAdaptiveEditDistance(this, findIndels, defaultEditDistance, minEditDistance, maxEditDistance);
		FindSimilarEntitiesResult<String, EditDistanceMappingMetric> result = collapseBarcodesGeneric(coreBarcodes, barcodes, function, 0);
		AdaptiveMappingResult r = new AdaptiveMappingResult(result.getEntityMap(), new ArrayList <EditDistanceMappingMetric>(result.getCollapseMetric()));
		return (r);		
	}
	
	public AdaptiveMappingResult collapseBarcodesAdaptive (final ObjectCounter<String> barcodes, final boolean findIndels, final int defaultEditDistance, final int minEditDistance, final int maxEditDistance) {
		List<String> coreBarcodes = barcodes.getKeysOrderedByCount(true);
		return (collapseBarcodesAdaptive(coreBarcodes, barcodes, findIndels, defaultEditDistance, minEditDistance, maxEditDistance));
	}
	
	/**
	 * Perform edit distance collapse where all barcodes are eligible to be collapse.
	 * Barcodes are ordered by total number of counts.
	 * @param barcodes
	 * @param findIndels
	 * @param editDistance The maximum edit distance from the parent barcode to any child.
	 * @param pathStepSize The maximum edit distance from any child to another child.
	 * @return
	 */
	public Map<String, List<String>> collapseBarcodesByMutationalCollapse (final ObjectCounter<String> barcodes, final boolean findIndels, final int maxEditDistance, final int minSizeToCollapse, final int pathStepSize) {
		FindSimilarEntitiesByMutationalCollapse function = new FindSimilarEntitiesByMutationalCollapse(this, findIndels, maxEditDistance, pathStepSize);
		List<String> orderedBarcodes = barcodes.getKeysOrderedByCount(true);
		FindSimilarEntitiesResult<String, String> result = collapseBarcodesGeneric(orderedBarcodes, barcodes, function, minSizeToCollapse);		
		return (result.getEntityMap());		
	}
	
	public FindSimilarEntitiesResult<String, UmiSharingMetrics> collapseBarcodesByUmiSharing (final ObjectCounter<String> barcodes, ParentEditDistanceMatcher parentEditDistanceMatcher, final double sharingThreshold, Map<String, Set<TagValues>> umisPerBarcode) {
		FindSimilarEntitiesByUMISharing function = new FindSimilarEntitiesByUMISharing(this, parentEditDistanceMatcher, sharingThreshold, umisPerBarcode);
		List<String> orderedBarcodes = barcodes.getKeysOrderedByCount(true);
		FindSimilarEntitiesResult<String, UmiSharingMetrics> result = collapseBarcodesGeneric(orderedBarcodes, barcodes, function, 0);		
		return result;	
	}
	
	public Set<String> processSingleBarcode(final String barcode, final List<String> comparisonBarcodes, final boolean findIndels, final int editDistance) {
		Set<String> closeBarcodes =null;

		// Replaced with java 8 lambda method. woot?
		if (this.NUM_THREADS>1 ) {
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
	public Set<String> processSingleBarcodeMultithreaded(final String barcode, final List<String> comparisonBarcodes, final boolean findIndels, final int editDistance) {
		Set<String> result = Collections.EMPTY_SET;
		try {
			if (findIndels)
				result = forkJoinPool.submit(() -> comparisonBarcodes.parallelStream().filter(x -> LevenshteinDistance.getIndelSlidingWindowEditDistance(barcode, x) <= editDistance).collect(Collectors.toSet())).get();
			else
				result = forkJoinPool.submit(() -> comparisonBarcodes.parallelStream().filter(x -> HammingDistance.getHammingDistance(barcode, x) <= editDistance).collect(Collectors.toSet())).get();
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}
		return result;
	}

	/**
	 * The Java lambda way.
	 * @param barcode
	 * @param comparisonBarcodes
	 * @param findIndels
	 * @param editDistance
	 * @return
	 */
	public int [] getEditDistanceDistributioneMultithreaded(final String barcode, final Collection<String> comparisonBarcodes, final boolean findIndels) {
		int [] result=null;
		try {
			if (findIndels)
				result = forkJoinPool.submit(() -> comparisonBarcodes.parallelStream().mapToInt(x -> LevenshteinDistance.getIndelSlidingWindowEditDistance(barcode, x)).toArray()).get();
			else
				result = forkJoinPool.submit(() -> comparisonBarcodes.parallelStream().mapToInt(x -> HammingDistance.getHammingDistance(barcode, x)).toArray()).get();
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}
		return result;
	}


	public class AdaptiveMappingResult {
		private final Map<String, List<String>> barcodeCollapseResult;
		private final List<EditDistanceMappingMetric> metricResult;

		public AdaptiveMappingResult (final Map<String, List<String>> barcodeCollapseResult, final List<EditDistanceMappingMetric> metricResult) {
			this.barcodeCollapseResult=barcodeCollapseResult;
			this.metricResult=metricResult;
		}

		public Map<String, List<String>> getBarcodeCollapseResult() {
			return barcodeCollapseResult;
		}

		public List<EditDistanceMappingMetric> getMetricResult() {
			return metricResult;
		}

	}
		
	/**
	 * Generic method to do barcode collapse.  Orders barcodes by count from largest to smallest, and collapses "close" barcodes as described by the FindCloseEntities class.
	 * @param coreBarcodes A list of core barcodes that other barcodes can be collapsed into.
	 * @param entityCounts Counts of each core barcode [umis,reads,etc] which orders how barcodes are collapsed.
	 * @param function Defines how a parent barcode and its children are close enough to each other for collapse
	 * @param minSizeToCollapse Do not collapse core barcodes with counts less than this number.
	 * @return
	 */
	private <T extends Comparable<T>,M> FindSimilarEntitiesResult<T,M> collapseBarcodesGeneric(List<T> coreBarcodes, ObjectCounter<T> entityCounts, FindSimilarEntities<T,M> function, final int minSizeToCollapse) {
		// don't allow side effects to modify input lists.
		coreBarcodes = new ArrayList<T>(coreBarcodes);
		entityCounts = new ObjectCounter<T>(entityCounts);

		FindSimilarEntitiesResult<T,M> result = new FindSimilarEntitiesResult<>();
		int count = 0;
		int numBCCollapsed=0;

		List<T> barcodeList = entityCounts.getKeysOrderedByCount(true);
		int coreBarcodeCount=coreBarcodes.size();
		long startTime = System.currentTimeMillis();
		while (coreBarcodes.isEmpty()==false) {
			T b = coreBarcodes.get(0);
			int barcodeSize=entityCounts.getCountForKey(b);
			if (barcodeSize<minSizeToCollapse)
				break;  			
			count++;
			coreBarcodes.remove(b);
			barcodeList.remove(b);
			FindSimilarEntitiesResult<T,M> singleResult = function.find(b, barcodeList, entityCounts);
			
			List<T> closeBC=singleResult.getEntityMap().get(b);		
			numBCCollapsed+=closeBC.size();
			
			if (closeBC.contains(b))
				log.error("Related entities should not contain tested entity");

			// this sorting is for consistency in output.
			List<T> closeBCList = new ArrayList<>(closeBC);
			Collections.sort(closeBCList);
			result.addMapping(b, closeBCList);
			result.addMetrics(singleResult.getCollapseMetric());
			
			barcodeList.removeAll(closeBC);
			coreBarcodes.removeAll(closeBC);
			if (this.REPORT_PROGRESS_INTERVAL!=0 && count % this.REPORT_PROGRESS_INTERVAL == 0) {
				if (entityCounts.getSize()>10000) log.info("Processed [" + count + "] records, totals BC Space left [" + barcodeList.size() +"]", " # collapsed this set [" + numBCCollapsed+"]");
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


	/******************************************************
	 * TO BE DELETED
	 * ****************************************************
	 **/
	
	/**
	 * Collapses a list of barcodes.
	 * This works by iterating through every core barcode (ordered from largest to smallest in the barcodes object), and mapping any other barcode close
	 * to this barcode.
	 *
	 * In the case of adaptive collapse, "close" uses the edit distance distribution of all barcodes to the core barcode being examined, and assumes that this
	 * distribution is bimodal.  The edit distance of lowest density is selected as the maximal edit distance.
	 *
	 * This means that the number of coreBarcodes at the end is less than or equal to the number of core barcodes submitted.  These barcodes are the keys of the output object.
	 * The barcodes that are not in the coreBarcodes list are eligible to be collapsed into a core barcode, but will never absorb other barcodes.
	 * We use coreBarcodes in order to limit the scope of the computational work, as the number of coreBarcodes can be small (the number of cells)
	 * compared to the total number of barcodes (number of beads + bead sequencing errors.)
	 * Smaller barcodes are always collapsed into larger ones.  Each barcode only exists once in the output - if a barcode A is edit distance 1 away
	 * from barcode B and C, it will be assigned to whichever of B and C is larger.
	 * @param coreBarcodes A list of barcode strings that are considered "core" or primary barcodes.
	 * @param barcodes An exhaustive list of all barcodes (both core and non-core) with assigned counts of observations of these barcodes.
	 * @param findIndels If true, we use Levenshtein indel sensitive collapse.  If false, use Hamming distance.
	 * @param defaultEditDistance If the discovered edit distance threshold is less than this number, this is used instead.  Set to 0 to effectively ignore the parameter.
	 * @return
	 */
	/*
	public AdaptiveMappingResult collapseBarcodesAdaptiveOld(List<String> coreBarcodes, ObjectCounter<String> barcodes, final boolean findIndels, final int defaultEditDistance, final int minEditDistance, final int maxEditDistance) {
		// don't allow side effects to modify input lists.
		coreBarcodes = new ArrayList<>(coreBarcodes);
		barcodes = new ObjectCounter<>(barcodes);

		Set<String> allBarcodes = new HashSet<>();
		allBarcodes.addAll(coreBarcodes);
		allBarcodes.addAll(barcodes.getKeys());

		Map<String, List<String>> result = new HashMap<>();
		int count = 0;
		int numBCCollapsed=0;

		List<String> barcodeList = barcodes.getKeysOrderedByCount(true);

		int coreBarcodeCount=coreBarcodes.size();
		long startTime = System.currentTimeMillis();
		List<EditDistanceMappingMetric> metrics = new ArrayList<>();

		while (coreBarcodes.isEmpty()==false) {
			String b = coreBarcodes.get(0);
			count++;
			coreBarcodes.remove(b);
			barcodeList.remove(b);
			// find the edit distance threshold.
			int [] edList = getEditDistanceDistributioneMultithreaded(b, allBarcodes, findIndels);
			// filter out ed=0.
			edList = Arrays.stream( edList ).filter(x-> x>0).toArray();
			int editDistanceDiscovered=findEditDistanceThreshold(edList, minEditDistance, maxEditDistance);

			// constrain to min/max edit distances.  Retain the original discovered edit distance.
			int editDistance=editDistanceDiscovered;
			//if (editDistance > maxEditDistance) editDistance=defaultEditDistance;
			if (editDistance > maxEditDistance | editDistanceDiscovered==-1) editDistance=defaultEditDistance;

			Set<String> closeBC=processSingleBarcode(b, barcodeList, findIndels, editDistance);

			numBCCollapsed+=closeBC.size();
			// Steve reports all barcodes, not just collapsed ones.
			int mergedObservations=getMergedNumObservations(b, closeBC, barcodes);
			// add the edit distance distribution for this entitiy to all other entities here.
			EditDistanceMappingMetric metric = new EditDistanceMappingMetric(b, closeBC.size(), editDistance, editDistanceDiscovered, barcodes.getCountForKey(b), mergedObservations, edList);
			metrics.add(metric);
			
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
		AdaptiveMappingResult r = new AdaptiveMappingResult(result, metrics);
		return (r);
	}
	*/
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
	/*
	public Map<String, List<String>> collapseBarcodesOld(List<String> coreBarcodes, ObjectCounter<String> barcodes, final boolean findIndels, final int editDistance) {
		// don't allow side effects to modify input lists.
		coreBarcodes = new ArrayList<>(coreBarcodes);
		barcodes = new ObjectCounter<>(barcodes);

		Map<String, List<String>> result = new HashMap<>();
		int count = 0;
		int numBCCollapsed=0;

		List<String> barcodeList = barcodes.getKeysOrderedByCount(true);
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
	

	private int getMergedNumObservations (final String originalBC, final Set<String> closeBC, final ObjectCounter<String> barcodes) {
		int total=barcodes.getCountForKey(originalBC);
		for (String c: closeBC)
			total+=barcodes.getCountForKey(c);
		return total;
	}
	*/
	
	/**
	 * Perform edit distance collapse where all barcodes are eligible to be collapse.
	 * Barcodes are ordered by total number of counts.
	 * @param barcodes
	 * @param findIndels
	 * @param editDistance The maximum edit distance from the parent barcode to any child.
	 * @param pathStepSize The maximum edit distance from any child to another child.
	 * @return
	 */
	/*
	public Map<String, List<String>> collapseBarcodesByMutationalCollapseOld (final ObjectCounter<String> barcodes, final boolean findIndels, final int maxEditDistance, final int minSizeToCollapse, final int pathStepSize) {
		
		Map<String, List<String>> result = new HashMap<>();
		int originalBarcodeCount=barcodes.getSize();
		
		int count = 0;
		int numBCCollapsed=0;
		int totalCollapsed=0;
		
		List<String> barcodeList = barcodes.getKeysOrderedByCount(true);
		
		long startTime = System.currentTimeMillis();
		if (this.REPORT_PROGRESS_INTERVAL!=0)
			log.info("Start of mutational barcode collapse for [", barcodes.getKeys().size()+"] barcodes with minimum parent size [", minSizeToCollapse, "] and max edit distance [", maxEditDistance+"]");
		while (barcodeList.isEmpty()==false) {
			String b = barcodeList.get(0);
			int barcodeSize=barcodes.getCountForKey(b);
			// if you've iterated past the smallest barcode you're willing to consider, break the loop.
			if (barcodeSize<minSizeToCollapse)
				break;  
			count++;
			barcodeList.remove(b);
			
			Set<String> closeBC=findRelatedBarcodesByMutationalCollapse(b, barcodeList, findIndels, maxEditDistance, pathStepSize);
			numBCCollapsed+=closeBC.size();
			totalCollapsed+=closeBC.size();
			if (result.containsKey(b))
				log.error("Result should never have core barcode");

			List<String> closeBCList = new ArrayList<>(closeBC);
			Collections.sort(closeBCList);
			result.put(b, closeBCList);

			barcodeList.removeAll(closeBC);			
			if (this.REPORT_PROGRESS_INTERVAL!=0 && count % this.REPORT_PROGRESS_INTERVAL == 0) {
				if (barcodes.getSize()>10000) log.info("Processed [" + count + "] records, totals BC Space left [" + barcodeList.size() +"]", " # collapsed this set [" + numBCCollapsed+"]");				
				numBCCollapsed=0;
				
			}
		}
		if (verbose) {
			long endTime = System.currentTimeMillis();
			long duration = (endTime - startTime)/1000;
			log.info("Collapse with [" + this.NUM_THREADS +"] threads took [" + duration + "] seconds to process");
			log.info("Started with core barcodes [" +originalBarcodeCount+  "] ended with [" + count + "] num collapsed [" +  (totalCollapsed) +"]");
		}
		return (result);
	}
	*/
	
	/**
	 * Find barcodes related to the core barcode.
	 * For minEditDistance=1:
	 * First find all barcodes related by ED=1.  Then find any barcodes that are related to the ED=1 barcodes by ED=1 AND related to the core barcode at ED=2.
	 * Iterate to the max edit distance.
	 *
	 * @param barcode The barcode to find neighbors for
	 * @param allBarcodes All barcodes to search
	 * @param findIndels Should we find indels, or hamming distance only?
	 * @param maxEditDistance The maximum edit distance to search
	 * @param pathStepSize The maximum edit distance from any child to another child.
	 * @return
	 */
	/*
	Set<String> findRelatedBarcodesByMutationalCollapse (final String barcode, final List<String> allBarcodes, final boolean findIndels, final int maxEditDistance, final int pathStepSize) {
		// parameterizing minEditDistance could lead to complications - how far apart are subseqeunt jumps from the first set of barcodes found?  ED=1 or ED=minEditDistance?		
		// store results for each edit distance here.
		Map<Integer, List<String>> validBarcodes = new HashMap<> ();

		// map from the initial barcode to all the other barcodes, group by edit distance.
		Map<Integer, Set<String>> barcodesAtED = new HashMap<>();
		int [] edList = getEditDistanceDistributioneMultithreaded(barcode, allBarcodes, findIndels);
		for (int i=0; i<allBarcodes.size(); i++) {
			int ed = edList[i];
			if (ed<=maxEditDistance) {
				Set<String> bcSet = barcodesAtED.get(ed);
				if (bcSet==null) {
					bcSet=new HashSet<>();
					barcodesAtED.put(ed, bcSet);
				}
				bcSet.add(allBarcodes.get(i));
			}
		}

		for (int editDistance=pathStepSize; editDistance<=maxEditDistance; editDistance+=pathStepSize) {
			// short circuit this edit distance if there are no barcodes at the edit distance.
			if (!barcodesAtED.containsKey(editDistance)) continue;

			// test barcodes in this iteration.  Either the starting barcode, or all neighbors that were editDistance-1 of the current distance.
			List<String> validBarcodesLastIteration = null;
			List<String>  barcodesToTest=new ArrayList<>(barcodesAtED.get(editDistance));

			// if we're on the first iteration, then all results at this edit distance hop are "valid" without further checks.
			if (editDistance==pathStepSize) {
				validBarcodes.put(editDistance, barcodesToTest);
				continue; // break out of loop, you're done.
			}

			// we're on some other iteration, we use the results of the last iteration as one of our two barcode lists.
			validBarcodesLastIteration=validBarcodes.get(editDistance-pathStepSize);
			if (validBarcodesLastIteration==null) validBarcodesLastIteration=Collections.EMPTY_LIST;

			// now test the valid barcodes against the barcodes to test.
			Set<String> newBarcodesThisIter=new HashSet<>();

			for (String bc: validBarcodesLastIteration) {
				// find barcodes ED=pathStepSize away from the last iteration results.
				Set<String> resultOneIter= processSingleBarcode(bc, barcodesToTest, findIndels, pathStepSize);
				newBarcodesThisIter.addAll(resultOneIter);
			}
			// iteration finished
			validBarcodes.put(editDistance, new ArrayList<>(newBarcodesThisIter));
		}
		// flatten out results from each edit distance iteration.
		Set<String> result = new HashSet<>();
		for (List<String> l: validBarcodes.values())
			result.addAll(l);

		return result;
	}
	*/

	 
	 
}
