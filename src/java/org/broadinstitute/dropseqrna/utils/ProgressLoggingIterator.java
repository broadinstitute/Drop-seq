package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ProgressLogger;

import java.util.Iterator;

/**
 * A simple iterator that logs each read that passes through it.
 * @author nemesh
 *
 */
public class ProgressLoggingIterator extends FilteredIterator<SAMRecord> {

	private final ProgressLogger progressLogger;

	public ProgressLoggingIterator (final Iterator<SAMRecord> underlyingIterator, final ProgressLogger progressLogger) {
		super(underlyingIterator);
		this.progressLogger=progressLogger;
	}

	@Override
	protected boolean filterOut(final SAMRecord rec) {
		this.progressLogger.record(rec);
		return false;
	}






}
