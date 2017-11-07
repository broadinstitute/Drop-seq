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
package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

public class Utils {
	private static Log log = Log.getInstance(Utils.class);
	private ProgressLogger progress = new ProgressLogger(log, 1000000);
	private static String DEFAULT_CELL_BARCODE="DEFAULT";

	public static SAMRecord getClone (final SAMRecord r) {
		SAMRecord rr=null;
		try {
			rr = (SAMRecord) r.clone();
		} catch (CloneNotSupportedException e) {
			log.info("This should never happen.  SAMRecord can't be cloned?");
		}
		return (rr);
	}

	public static String getCellBC (final SAMRecord r, final String cellBCTag) {
		String currentCell = r.getStringAttribute(cellBCTag);
		if (currentCell==null) return (DEFAULT_CELL_BARCODE);
		return (currentCell);
	}

	public static String strandToString(final boolean strand) {
		if (strand) return "+";
		return "-";
	}

}
