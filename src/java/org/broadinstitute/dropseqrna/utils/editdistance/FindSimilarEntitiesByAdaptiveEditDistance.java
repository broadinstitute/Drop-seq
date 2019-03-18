package org.broadinstitute.dropseqrna.utils.editdistance;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;

public class FindSimilarEntitiesByAdaptiveEditDistance implements FindSimilarEntities<String, EditDistanceMappingMetric> {

	private final MapBarcodesByEditDistance mbed;
	private final boolean findIndels;
	private final int defaultEditDistance;
	private final int minEditDistance;
	private final int maxEditDistance;
	
	
	public FindSimilarEntitiesByAdaptiveEditDistance(final MapBarcodesByEditDistance mbed, final boolean findIndels, final int defaultEditDistance, final int minEditDistance, final int maxEditDistance) {
		this.mbed=mbed;
		this.findIndels=findIndels;
		this.defaultEditDistance=defaultEditDistance;	
		this.minEditDistance=minEditDistance;
		this.maxEditDistance=maxEditDistance;
	}


	@Override
	/**
	 * Adaptive collapse. The definition of similar in this strategy uses the edit distance distribution of all barcodes to the core barcode being examined, and assumes that this
	 * distribution is bimodal.  The edit distance of lowest density is selected as the maximal edit distance, to separate some group of closely related barcodes from another set that 
	 * is much further away.
	 * 
	 * @param entity A barcode that may have one or more children in the search space.
	 * @param searchSpace A list of other barcodes, a subset of which may be related to entity.
	 * @param counts The number of observations of each barcode (usually a UMI count).  This orders barcodes from largest to smallest for collapse.
	 */
	public FindSimilarEntitiesResult<String, EditDistanceMappingMetric> find(String entity, List<String> searchSpace, ObjectCounter<String> counts) {
		int [] edList = mbed.getEditDistanceDistributioneMultithreaded(entity, searchSpace, findIndels);
		// filter out ed=0.
		edList = Arrays.stream( edList ).filter(x-> x>0).toArray();
		int editDistanceDiscovered=findEditDistanceThreshold(edList);

		// constrain to min/max edit distances.  Retain the original discovered edit distance.
		int editDistance=editDistanceDiscovered;
		if (editDistance > maxEditDistance | editDistanceDiscovered==-1) editDistance=defaultEditDistance;

		Set<String> closeBC=mbed.processSingleBarcode(entity, searchSpace, findIndels, editDistance);		
		// Steve reports all barcodes, not just collapsed ones.
		int mergedObservations=getMergedNumObservations(entity, closeBC, counts);
		// add the edit distance distribution for this entity to all other entities here.
		EditDistanceMappingMetric metric = new EditDistanceMappingMetric(entity, closeBC.size(), editDistance, editDistanceDiscovered, counts.getCountForKey(entity), mergedObservations, edList);
		FindSimilarEntitiesResult<String, EditDistanceMappingMetric> result = new FindSimilarEntitiesResult<>();
		result.addMapping(entity, closeBC);
		result.addMetrics(metric);
		return (result);		
	}
	
	private int getMergedNumObservations (final String originalBC, final Set<String> closeBC, final ObjectCounter<String> barcodes) {
		int total=barcodes.getCountForKey(originalBC);
		for (String c: closeBC)
			total+=barcodes.getCountForKey(c);
		return total;
	}
	
	/**
	 * 
	 * @param targetBarcode
	 * @param allBarcodes
	 * @return
	 */
	public Integer findEditDistanceThreshold (final String targetBarcode, final Set<String> allBarcodes) {		
		int [] edList = mbed.getEditDistanceDistributioneMultithreaded(targetBarcode, allBarcodes, findIndels);
		// build and populate an object counter.
		ObjectCounter<Integer> counts = new ObjectCounter<>();
		Arrays.stream(edList).forEach(x -> counts.increment(x));
		// System.out.println(Arrays.toString(edList));

		// loop through to find the minimum count point, from the lowest edit distance to the maximum ED.
		// int maxED = Arrays.stream(edList).max().getAsInt();
		List<Integer> result = new ArrayList<>();
		int resultCount=Integer.MAX_VALUE;
		// find the minimum position in the object counter.  This may include keys that don't exist (ie: they are 0.)
		// if there are multiple minimums, then select the lowest.
		for (int i=minEditDistance; i<=maxEditDistance; i++) {
			int count = counts.getCountForKey(i);
			if (count==0) {
				resultCount=count;
				result.add(i);
			}
		}
		// no results found.
		if (result.size()==0) return -1;
		// special case: no entries within all of scanned results.  return -1.
		int numPosScanned=(maxEditDistance-minEditDistance)+1;
		if (result.size()==numPosScanned && resultCount==0) return -1;

		// there may be ties.  select the first entry (smallest ED)
		Collections.sort(result);
		int threshold=result.get(0);
		threshold=threshold-1; // steve reports the last filled position, not the empty one.
		return threshold;
	}

	public Integer findEditDistanceThreshold (final int [] edList) {
		// build and populate an object counter.
		ObjectCounter<Integer> counts = new ObjectCounter<>();
		Arrays.stream(edList).forEach(x -> counts.increment(x));
		// System.out.println(Arrays.toString(edList));

		// loop through to find the minimum count point, from the lowest edit distance to the maximum ED.
		// int maxED = Arrays.stream(edList).max().getAsInt();
		List<Integer> result = new ArrayList<>();
		int resultCount=Integer.MAX_VALUE;
		// find the minimum positions in the object counter.  This may include keys that don't exist (ie: they are 0.)
		// if there are multiple minimums, then select the lowest.
		for (int i=minEditDistance; i<=maxEditDistance; i++) {
			int count = counts.getCountForKey(i);
			// a lower or equal count.
			if (count<=resultCount) {
				// if this is a lower count, then it's the only result.
				// if this is the same count, then add it to the result.
				if (count<resultCount) result.clear();
				result.add(i);
				resultCount=count;

			}
		}
		// no results found.
		if (result.size()==0) return -1;
		// special case: no entries within all of scanned results.  return -1.
		int numPosScanned=(maxEditDistance-minEditDistance)+1;
		if (result.size()==numPosScanned && resultCount==0) return -1;

		// there may be ties.  select the first entry (smallest ED)
		Collections.sort(result);
		int threshold=result.get(0);
		threshold=threshold-1; // steve reports the last filled position, not the empty one.
		return threshold;
	}

			
}
