package org.broadinstitute.dropseqrna.utils.editdistance;

import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;

public class FindSimilarEntitiesByEditDistance implements FindSimilarEntities<String, String> {

	private final MapBarcodesByEditDistance mbed;
	private final boolean findIndels;
	private final int editDistance;
	
	public FindSimilarEntitiesByEditDistance(final MapBarcodesByEditDistance mbed, final boolean findIndels, final int editDistance) {
		this.mbed=mbed;
		this.findIndels=findIndels;
		this.editDistance=editDistance;		
	}
	
	@Override
	public FindSimilarEntitiesResult<String,String> find (String entity, List<String> searchSpace, ObjectCounter<String> counts) {
		Set<String> r= mbed.processSingleBarcode(entity, searchSpace, findIndels, editDistance);
		FindSimilarEntitiesResult<String,String> rr = new FindSimilarEntitiesResult<>();
		rr.addMapping(entity, r);				
		return rr;		
	}

}
