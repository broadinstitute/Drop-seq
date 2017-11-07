/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
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
