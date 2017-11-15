/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
