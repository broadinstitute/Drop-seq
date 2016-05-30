package org.broadinstitute.dropseqrna.utils.modularfileparser;




/**
 * Produces a modular file parser with the parser class and options set.
 * @author nemesh
 *
 */
/*
public class ParserFactory {

	private enum RegisteredParsers {
		DELIMITED("delimited", "delimited (csv, tab, etc) file parser");
		
		private final String name;
		private final String description;
		
		RegisteredParsers(String n, String d) {
			this.name = n;
			this.description=d;
		}
		
		public String getName () {
			return this.name;
		}
		
		public String getDescription () {
			return this.description;
		}
	}
	
	private static class SingletonHolder {
		public static final ParserFactory INSTANCE = new ParserFactory();
	
	}
	
	private ParserFactory() {}

	public static ParserFactory getInstance() {
		return SingletonHolder.INSTANCE;
	}
	
	public Parser getDelimitedParser(String delimiter) {
		DelimiterParser p = new DelimiterParser(delimiter);
		return (p);
	}
	
	public List<String> getFileParserList() {
		List<String> parsers=new ArrayList();
		for (RegisteredParsers p: RegisteredParsers.values()) {
			parsers.add(p.getName());
		}
		return parsers;
		
	}
}
*/
