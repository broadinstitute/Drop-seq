package org.broadinstitute.dropseqrna.cmdline;

import picard.cmdline.CommandLineProgramGroup;

public class DropNet implements CommandLineProgramGroup {
	@Override
	public String getName() {
		return "DropNet Tools";
	}

	@Override
	public String getDescription() {
		return "Tools for aligning or analyzing DropNet experiments.";
	}
}



