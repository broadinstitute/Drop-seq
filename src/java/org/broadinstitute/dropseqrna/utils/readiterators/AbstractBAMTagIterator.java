package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;

import java.io.File;
import java.util.Iterator;
import java.util.List;

public abstract class AbstractBAMTagIterator<T> implements Iterator<T>{
	
	public DEIteratorUtils iterUtils;
	public PeekableIterator<SAMRecord> iter;
	private static Log log = Log.getInstance(AbstractBAMTagIterator.class);
	
	public AbstractBAMTagIterator (int maxRecordsInRAM) {
		this.iterUtils = new DEIteratorUtils(maxRecordsInRAM);		
	}
	
	public AbstractBAMTagIterator () {
		this.iterUtils = new DEIteratorUtils();		
	}
	
	public boolean hasNext() {
		return this.iter.hasNext();
	}
	
	public void remove () {
		throw new UnsupportedOperationException("Remove not supported");
	}
	
	
	private void initialize (CloseableIterator<SAMRecord> cIter, List<String> sortingTags) {
		this.iter = new PeekableIterator<SAMRecord>(cIter);
		
		iter = iterUtils.primeIterator (iter, sortingTags);
		if (iter.hasNext()==false) {
			log.error("No records found to iterate on.  Did you select the right BAM tags for your cell barcode, and gene/exon tag?");
		}
	}
	
	void initialize (File bamFile, List<String> sortingTags, SAMReadProcessorI processor) {
		CloseableIterator<SAMRecord> cIter = iterUtils.getReadsInTagOrder(bamFile, sortingTags, processor);
		initialize(cIter, sortingTags);		
	}
	
	void initialize (File bamFile, List<String> sortingTags, ReadProcessorCollection rpc) {
		CloseableIterator<SAMRecord> cIter = iterUtils.getReadsInTagOrder(bamFile, sortingTags, rpc);
		initialize(cIter, sortingTags);
	}
	
	public abstract T next ();
	
	
	
}
