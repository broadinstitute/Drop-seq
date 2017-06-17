package org.broadinstitute.dropseqrna.utils;

import java.util.Iterator;

import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantContextProgressLoggerIterator extends FilteredIterator<VariantContext>{
	private final ProgressLogger progressLogger;

	public VariantContextProgressLoggerIterator (final Iterator<VariantContext> underlyingIterator, final ProgressLogger progressLogger) {
		super(underlyingIterator);
		this.progressLogger=progressLogger;
	}

	@Override
	public boolean filterOut(final VariantContext rec) {
		this.progressLogger.record(rec.getContig(), rec.getStart());
		return false;
	}
}
