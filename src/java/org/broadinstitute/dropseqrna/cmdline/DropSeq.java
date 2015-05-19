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
