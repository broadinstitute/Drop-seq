package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Iterator;

import org.broadinstitute.dropseqrna.utils.TransformingIterator;

import htsjdk.samtools.SAMRecord;

/**
 * All reads receive the default value for the outTag, which is the collapse tag's current value
 * @author nemesh
 *
 */
public class DefaultTaggingIterator extends TransformingIterator<SAMRecord,SAMRecord> {

	private final String collapseTag;
	private final String outTag;
	
	public DefaultTaggingIterator(Iterator<SAMRecord> underlyingIterator, final String collapseTag, final String outTag) {
		super(underlyingIterator);
		this.collapseTag=collapseTag;
		this.outTag=outTag;
	}

	@Override
	public SAMRecord next() {
		SAMRecord r= this.underlyingIterator.next();
		String v = r.getStringAttribute(collapseTag);
		if (v!=null)
			r.setAttribute(outTag, v);
		return (r);
	}
	
}
