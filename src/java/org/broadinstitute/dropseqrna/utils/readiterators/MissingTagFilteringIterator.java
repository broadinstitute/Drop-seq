/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2015 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Iterator;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTagUtil;

/**
 * Iterator wrapper that emits a SAMRecord only if *all* the required tags are present.
 */
public class MissingTagFilteringIterator extends FilteredIterator<SAMRecord> {
    final short[] requiredTags;

    public MissingTagFilteringIterator(final Iterator<SAMRecord> underlyingIterator, final String... requiredTags) {
        super(underlyingIterator);
        this.requiredTags = new short[requiredTags.length];
        for (int i = 0; i < requiredTags.length; ++i)
			this.requiredTags[i] = SAMTagUtil.getSingleton().makeBinaryTag(requiredTags[i]);
    }

    @Override
    public boolean filterOut(final SAMRecord rec) {
        for (final short tag : requiredTags)
			if (rec.getAttribute(tag) == null)
				return true;
        return false;
    }
}
