package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;

import java.util.ArrayList;
import java.util.List;

public class Utils {
	private static Log log = Log.getInstance(Utils.class);
	private ProgressLogger progress = new ProgressLogger(log, 1000000);
	private static String DEFAULT_CELL_BARCODE="DEFAULT";
	
	public static SAMRecord getClone (SAMRecord r) {
		SAMRecord rr=null;
		try {
			rr = (SAMRecord) r.clone();
		} catch (CloneNotSupportedException e) {
			log.info("This should never happen.  SAMRecord can't be cloned?");
		}
		return (rr);
	}
	
	public static String getCellBC (SAMRecord r, String cellBCTag) {
		String currentCell = r.getStringAttribute(cellBCTag);
		if (currentCell==null) return (DEFAULT_CELL_BARCODE);
		return (currentCell);
	}
	
	public static String strandToString(boolean strand) {
		if (strand) return "+";
		return "-";
	}
		
}
