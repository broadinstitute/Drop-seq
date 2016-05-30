package org.broadinstitute.dropseqrna.utils.modularfileparser;

public class DelimiterParser implements Parser {

	private String delimiter;
	
	public DelimiterParser(String delimiter) {
		this.delimiter=delimiter;
	}
	
	@Override
	public String[] parseLine(String nextLine) {
		return nextLine.split(delimiter);
	}
	
}
