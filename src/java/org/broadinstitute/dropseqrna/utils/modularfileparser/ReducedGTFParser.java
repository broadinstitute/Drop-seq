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
