package org.broadinstitute.dropseqrna.cmdline;

import picard.cmdline.CommandLineProgramGroup;

public class MetaData implements CommandLineProgramGroup {
	@Override
	public String getName() {
		return "MetaData Tools";
	}

	@Override
	public String getDescription() {
		return "Tools to generate or manipulate various kinds of meta data";
	}
}
