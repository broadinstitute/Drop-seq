package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

public class ReadProcessorCollection {

	List<SAMReadProcessorI> processorList = null;
	
	public ReadProcessorCollection () {
		this.processorList = new ArrayList<SAMReadProcessorI>();
	}
	
	public ReadProcessorCollection (Collection <SAMReadProcessorI> processors) {
		this.processorList = new ArrayList<SAMReadProcessorI>();
		this.processorList.addAll(processors);
	}

	public void addFilter (SAMReadProcessorI processor) {
		this.processorList.add(processor);
	}
	
	// Pass in a collection of samrecords that's empty.
	// shove records into the temp list.
	// swap the result with the temp list
	// empty the temp list.
	public Collection<SAMRecord> processRead (SAMRecord r) {
		Collection<SAMRecord> result = new ArrayList<SAMRecord>();
		result.add(r);
		Collection<SAMRecord> temp = new ArrayList<SAMRecord>();
		
		for (SAMReadProcessorI p: processorList) {
			Collection<SAMRecord> processorResult = new ArrayList<SAMRecord>();
			for (SAMRecord rec: result) {
				temp = p.processRead(rec, temp);
				processorResult.addAll(temp);
				temp.clear();
			}
			// a processor is done - push result into current and loop again.
			result=processorResult;
		}
		
		return (result);
	}
	
	
	
	
	
	
}
