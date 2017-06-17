package org.broadinstitute.dropseqrna.utils;

import java.util.Iterator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ProgressLogger;

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
	public boolean filterOut(final SAMRecord rec) {
		this.progressLogger.record(rec);
		return false;
	}






}
