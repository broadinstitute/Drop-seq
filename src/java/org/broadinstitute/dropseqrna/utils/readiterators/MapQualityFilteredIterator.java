package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

public class MapQualityFilteredIterator extends FilteredIterator<SAMRecord> {

	private final Integer mapQuality;
	private final boolean rejectNonPrimaryReads;
	
	public MapQualityFilteredIterator(final Iterator<SAMRecord> underlyingIterator, final Integer mapQuality, final boolean rejectNonPrimaryReads) {
        super(underlyingIterator);
		this.mapQuality=mapQuality;
		this.rejectNonPrimaryReads=rejectNonPrimaryReads;
	}

    @Override
    protected boolean filterOut(SAMRecord r) {
        return (rejectNonPrimaryReads && r.isSecondaryOrSupplementary()) ||
                (this.mapQuality!=null && r.getMappingQuality() < this.mapQuality);
    }
}
