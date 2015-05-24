package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;

import java.util.Collection;

/**
 * Classes implementing this interface examine the read and decide if the read should be returned 0 times, once, or many times.
 * Alterations to the read can be made through the implementing class.
 * @author nemesh
 *
 */
public interface SAMReadProcessorI {

	public Collection<SAMRecord> processRead (SAMRecord r);
	
}
