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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

public class SamHeaderAndIterator {
    public final SAMFileHeader header;
    public final CloseableIterator<SAMRecord> iterator;

    public SamHeaderAndIterator(SAMFileHeader header, CloseableIterator<SAMRecord> iterator) {
        this.header = header;
        this.iterator = iterator;
    }
}
