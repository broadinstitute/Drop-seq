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
	public Collection<SAMRecord> processRead(SAMRecord r, Collection<SAMRecord> outList) {
		
		
		// drop read if we filter on non-primary and we see a non-primary read.
		if (rejectNonPrimaryReads && r.isSecondaryOrSupplementary()) return (outList);
		// reject on map quality if quality score is set.
		if (this.mapQuality!=null && r.getMappingQuality() < this.mapQuality) return (outList);
		
		// you're all good
		outList.add(r);
		return outList;
	}

}
