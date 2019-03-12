package org.broadinstitute.dropseqrna.utils.editdistance;

import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;

public class FindCloseEntitiesByEditDistance implements FindCloseEntities<String, String> {

	private final MapBarcodesByEditDistance mbed;
	private final boolean findIndels;
	private final int editDistance;
	
	public FindCloseEntitiesByEditDistance(final MapBarcodesByEditDistance mbed, final boolean findIndels, final int editDistance) {
		this.mbed=mbed;
		this.findIndels=findIndels;
		this.editDistance=editDistance;		
	}
	
	@Override
	public FindCloseEntitiesResult<String,String> find (String entity, List<String> searchSpace, ObjectCounter<String> counts) {
		Set<String> r= mbed.processSingleBarcode(entity, searchSpace, findIndels, editDistance);
		FindCloseEntitiesResult<String,String> rr = new FindCloseEntitiesResult<>();
		rr.addMapping(entity, r);				
		return rr;		
	}

}
