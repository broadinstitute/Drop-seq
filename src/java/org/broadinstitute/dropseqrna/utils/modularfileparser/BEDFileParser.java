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

import java.util.Arrays;

public class BEDFileParser implements Parser {

	private String delimiter;
	private boolean firstBodyLineRetrieved=false;
	int numCols=0;
	private String [] fullHeader = {"chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"};
	
	public BEDFileParser (String delimiter) {
		this.delimiter= delimiter;
	}
	
	@Override
	/**
	 * Bed files are tab delimited by default.
	 */
	public String[] parseLine(String nextLine) {
		String [] result = parse(nextLine); 
		if (firstBodyLineRetrieved==false) {
			firstBodyLineRetrieved=true;
			this.numCols=result.length;
		}
		return result;
	}

	
	/**
	 * Get the default header for a BEDFile, since the bedFile typically doesn't have a set of column definitions.
	 * Because the number of columns is variable, a set list of all the available columns is restricted to the number of results in the body of the data.
	 * After the first line of the file is parsed, you can get the header subsetted correctly.
	 */
	public String[] getDefaultHeader() {
		String [] result = Arrays.copyOfRange(this.fullHeader, 0, this.numCols);
		return (result);
	}
		
	private String [] parse(String nextLine) {
		String [] result = nextLine.split(this.delimiter);
		return result;
	}

	
}
