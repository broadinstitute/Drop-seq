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
