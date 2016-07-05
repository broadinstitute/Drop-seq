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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

public class DgeHeaderCommand {
    private String commandLine;

    public DgeHeaderCommand(String commandLine) {
        if (commandLine.contains("\t")) {
            throw new IllegalArgumentException("Command string may not contain tab character");
        }
        this.commandLine = commandLine;
    }

    public String getCommandLine() {
        return commandLine;
    }

    public void setCommandLine(String commandLine) {
        if (commandLine.contains("\t")) {
            throw new IllegalArgumentException("Command string may not contain tab character");
        }
        this.commandLine = commandLine;
    }
}
