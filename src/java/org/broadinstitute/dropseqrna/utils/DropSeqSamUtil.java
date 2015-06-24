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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;

import java.io.File;

public class DropSeqSamUtil {

    public static SAMSequenceDictionary loadSequenceDictionary(final File file) {
        final SamReader r = SamReaderFactory.makeDefault().open(file);
        SAMSequenceDictionary dict = r.getFileHeader().getSequenceDictionary();
        CloserUtil.close(r);
        return dict;
    }
}
