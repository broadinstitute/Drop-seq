/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2016 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import picard.cmdline.CommandLineProgram;

public class SamHeaderUtil {

    public static void addPgRecord(final SAMFileHeader header, final CommandLineProgram clp) {
        int highestIdSeen = -1;
        for (final SAMProgramRecord pg : header.getProgramRecords()) {
            try {
                int id = Integer.parseInt(pg.getId());
                if (id > highestIdSeen) {
                    highestIdSeen = id;
                }
            } catch (NumberFormatException e) {
                // ignore any non-integer PG IDs
            }
        }
        final SAMProgramRecord pg = new SAMProgramRecord(Integer.toString(highestIdSeen + 1));
        pg.setProgramName(clp.getClass().getSimpleName());
        pg.setCommandLine(clp.getCommandLine());
        pg.setProgramVersion(clp.getVersion());
        header.addProgramRecord(pg);
    }
}
