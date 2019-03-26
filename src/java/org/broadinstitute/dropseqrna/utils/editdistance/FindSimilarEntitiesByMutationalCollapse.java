package org.broadinstitute.dropseqrna.utils.editdistance;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;

public class FindSimilarEntitiesByMutationalCollapse implements FindSimilarEntities<String, String>  {
	private final MapBarcodesByEditDistance mbed;
	private boolean findIndels;
	private int maxEditDistance;
	private int pathStepSize;
	
	/**
	 * First find all barcodes related by edit distance<=pathStepSize.  Then find any barcodes that are related to that set of barcodes by ED=pathStepSize AND related to the core barcode at ED=2*path step size..
	 * Iterate to the max edit distance.
	 * 
	 * @param mbed A mapping object to perform string functions
	 * @param findIndels Uses levenstein distance when set to true, or hamming distance when set to false
	 * @param maxEditDistance The maximum edit distance from the original core barcode to the most distant barcode.
	 * @param pathStepSize The edit distance between barcodes.
	 */
	public FindSimilarEntitiesByMutationalCollapse (final MapBarcodesByEditDistance mbed, final boolean findIndels, final int maxEditDistance, final int pathStepSize) {
		this.mbed=mbed;
		this.findIndels=findIndels;
		this.maxEditDistance=maxEditDistance;
		this.pathStepSize=pathStepSize;
	}
	

	@Override
	/**
	 * 
	 */
	public FindSimilarEntitiesResult<String, String> find(String entity, List<String> searchSpace, ObjectCounter<String> counts) {
		// parameterizing minEditDistance could lead to complications - how far apart are subseqeunt jumps from the first set of barcodes found?  ED=1 or ED=minEditDistance?		
		// store results for each edit distance here.
		Map<Integer, List<String>> validBarcodes = new HashMap<> ();

		// map from the initial barcode to all the other barcodes, group by edit distance.
		Map<Integer, Set<String>> barcodesAtED = new HashMap<>();
		int [] edList = mbed.getEditDistanceDistributioneMultithreaded(entity, searchSpace, findIndels);
		for (int i=0; i<searchSpace.size(); i++) {
			int ed = edList[i];
			if (ed<=maxEditDistance) {
				Set<String> bcSet = barcodesAtED.get(ed);
				if (bcSet==null) {
					bcSet=new HashSet<>();
					barcodesAtED.put(ed, bcSet);
				}
				bcSet.add(searchSpace.get(i));
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
				Set<String> resultOneIter= mbed.processSingleBarcode(bc, barcodesToTest, findIndels, pathStepSize);
				newBarcodesThisIter.addAll(resultOneIter);
			}
			// iteration finished
			validBarcodes.put(editDistance, new ArrayList<>(newBarcodesThisIter));
		}
		// flatten out results from each edit distance iteration.
		Set<String> result = new HashSet<>();
		for (List<String> l: validBarcodes.values())
			result.addAll(l);

		FindSimilarEntitiesResult<String,String> rr = new FindSimilarEntitiesResult<>();
		rr.addMapping(entity, result);				
		return rr;
	}

	
}
