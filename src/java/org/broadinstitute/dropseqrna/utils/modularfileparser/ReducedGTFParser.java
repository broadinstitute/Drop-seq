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
package org.broadinstitute.dropseqrna.utils.modularfileparser;


public class ReducedGTFParser implements Parser {

	private final String delimiter="\t";
	
	public ReducedGTFParser () {
	}
	
	@Override
	/**
	 * Bed files are tab delimited by default.
	 */
	public String[] parseLine(String nextLine) {
		String [] result = nextLine.split(this.delimiter);
		return result;
	}

	
	
}
