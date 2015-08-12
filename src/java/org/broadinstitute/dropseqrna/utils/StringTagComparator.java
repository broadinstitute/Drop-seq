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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTagUtil;

import java.util.Comparator;

/**
 * Comparator for SAMRecord that orders by the value for the given tag, which must be a string value.
 * Not-present value for the tag sorts before any non-null value.
 */
public class StringTagComparator implements Comparator<SAMRecord> {
    private final short tag;

    public StringTagComparator(String tag) {
        this.tag = SAMTagUtil.getSingleton().makeBinaryTag(tag);
    }

    @Override
    public int compare(SAMRecord rec1, SAMRecord rec2) {
        final String s1 = (String)rec1.getAttribute(tag);
        final String s2 = (String)rec2.getAttribute(tag);

        if (s1 != null) {
            if (s2 == null)
                return 1;
            else {
                return s1.compareTo(s2);
            }
        } else if (s2 != null) {
            return -1;
        } else {
            return 0;
        }
    }
}
