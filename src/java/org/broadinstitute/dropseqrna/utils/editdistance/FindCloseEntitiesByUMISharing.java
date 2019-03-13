package org.broadinstitute.dropseqrna.utils.editdistance;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.dropseqrna.metrics.UmiSharingMetrics;
import org.broadinstitute.dropseqrna.metrics.umisharing.ParentEditDistanceMatcher;
import org.broadinstitute.dropseqrna.metrics.umisharing.ParentEditDistanceMatcher.TagValues;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

public class FindCloseEntitiesByUMISharing implements FindCloseEntities<String, UmiSharingMetrics>{

	private final MapBarcodesByEditDistance mbed;
	private final ParentEditDistanceMatcher parentEditDistanceMatcher;
	private final double sharingThreshold;
	private final Map<String, Set<TagValues>> umisPerBarcode;
	
	public FindCloseEntitiesByUMISharing (MapBarcodesByEditDistance mbed, ParentEditDistanceMatcher parentEditDistanceMatcher, final double sharingThreshold, Map<String, Set<TagValues>> umisPerBarcode) {
		this.mbed=mbed;
		this.parentEditDistanceMatcher=parentEditDistanceMatcher;
		this.sharingThreshold=sharingThreshold;
		this.umisPerBarcode=umisPerBarcode;
	}
	
	@Override
	public FindCloseEntitiesResult<String, UmiSharingMetrics> find(String entity, List<String> searchSpace, ObjectCounter<String> counts) {
		FindCloseEntitiesResult<String, UmiSharingMetrics> result = new FindCloseEntitiesResult();
		
		Set<TagValues> parentTuples=umisPerBarcode.get(entity);
		Set<String> resultBarcodes = new HashSet<String>();
		
		for (String child: searchSpace) {
			UmiSharingMetrics metrics = new UmiSharingMetrics();
			Set<TagValues> childTuples=umisPerBarcode.get(child);						
			metrics.PARENT = entity;
            metrics.CHILD = child;
            metrics.NUM_PARENT = parentTuples.size();
            metrics.NUM_CHILD = childTuples.size();
            metrics.NUM_SHARED = parentEditDistanceMatcher.computeNumShared(parentTuples, childTuples);
            metrics.FRAC_SHARED = metrics.NUM_SHARED/(double)metrics.NUM_CHILD;
            if (metrics.FRAC_SHARED>=this.sharingThreshold)
            	resultBarcodes.add(child);
            result.addMetrics(metrics);
		}
		result.addMapping(entity, resultBarcodes);
		return (result);
	}

}
