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
package org.broadinstitute.dropseqrna.cmdline;

import picard.cmdline.CommandLineProgramGroup;

public class DropSeq implements CommandLineProgramGroup {
	@Override
	public String getName() {
		return "DropSeq Tools";
	}

	@Override
	public String getDescription() {
		return "Tools for aligning or analyzing DropSeq experiments.";
	}
}
