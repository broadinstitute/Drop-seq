package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;

import java.io.File;
import java.util.Collection;
import java.util.List;

import org.broadinstitute.dropseqrna.barnyard.Utils;

public class SimpleTagIterator extends AbstractBAMTagIterator <Collection<SAMRecord>> {
	
	private static Log log = Log.getInstance(SimpleTagIterator.class);
	
	private List<String>sortingTags;
	
	public SimpleTagIterator (File bamFile, List<String>sortingTags, Integer mapQuality, boolean rejectNonPrimaryReads) {
		super();
		this.sortingTags=sortingTags;
		MapQualityProcessor f = new MapQualityProcessor(mapQuality, rejectNonPrimaryReads);
		super.initialize(bamFile, sortingTags, f);
		
	}
	
	public SimpleTagIterator (File bamFile, List<String>sortingTags, ReadProcessorCollection rpc) {
		this.sortingTags=sortingTags;
		super.initialize(bamFile, sortingTags, rpc);
		
	}
	
	public Collection<SAMRecord> next() {
		Collection<SAMRecord> records = iterUtils.getRecordCollection(this.iter, this.sortingTags);
		return (records);
	}
	
	
}
