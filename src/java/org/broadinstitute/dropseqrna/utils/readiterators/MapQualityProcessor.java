package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class MapQualityProcessor implements SAMReadProcessorI {

	private Integer mapQuality;
	private boolean rejectNonPrimaryReads;
	
	public MapQualityProcessor (Integer mapQuality, boolean rejectNonPrimaryReads) {
		this.mapQuality=mapQuality;
		this.rejectNonPrimaryReads=rejectNonPrimaryReads;
	}
	
	@Override
	public Collection<SAMRecord> processRead(SAMRecord r) {
		List<SAMRecord> result = new ArrayList<SAMRecord>(1);
		
		// drop read if we filter on non-primary and we see a non-primary read.
		if (rejectNonPrimaryReads && r.isSecondaryOrSupplementary()) return (result);
		// reject on map quality if quality score is set.
		if (this.mapQuality!=null && this.mapQuality < r.getMappingQuality()) return (result);
		
		// you're all good
		result.add(r);
		return result;
	}

}
