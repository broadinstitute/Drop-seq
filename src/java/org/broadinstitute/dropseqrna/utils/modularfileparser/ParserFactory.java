/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
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
