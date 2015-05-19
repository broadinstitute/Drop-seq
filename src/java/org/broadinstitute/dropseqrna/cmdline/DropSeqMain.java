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
package org.broadinstitute.dropseqrna.cmdline;

import picard.cmdline.PicardCommandLine;

import java.util.ArrayList;
import java.util.List;

public class DropSeqMain extends PicardCommandLine {
    /** The name of this unified command line program **/
    private final static String COMMAND_LINE_NAME = DropSeqMain.class.getSimpleName();

    /** The packages we wish to include in our command line **/
    protected static List<String> getPackageList() {
        final List<String> packageList = new ArrayList<String>();
        packageList.add("org.broadinstitute.dropseqrna");
        return packageList;
    }
    public static void main(final String[] args) {
        System.exit(new DropSeqMain().instanceMain(args, getPackageList(), COMMAND_LINE_NAME));
    }

}
