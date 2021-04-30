package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import htsjdk.samtools.util.Interval;

public interface SNPIntervalRecordI {

	public Interval getSNPInterval ();
	public int getNumBases();
}
