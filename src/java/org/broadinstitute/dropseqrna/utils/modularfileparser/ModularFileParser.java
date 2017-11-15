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

import htsjdk.samtools.util.Log;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class ModularFileParser  {
	private final Log log = Log.getInstance(ModularFileParser.class);

	private Parser parser;
	private BufferedReader in;
	private boolean hasNext = true;
    private int linesToSkip;
    private boolean linesSkiped;
	
	public ModularFileParser(Parser parser, File inFile, int linesToSkip)  {
		this.parser=parser;
		try {
			in = new BufferedReader(new FileReader(inFile));
		} catch (FileNotFoundException e) {
			throw new ModularFileParserException("File not found: " + inFile.getAbsolutePath());
		}
		this.linesToSkip=linesToSkip;
	}
	
	
	public String[] readNextLine()  {
		String rawLine = getNextRawLine();
		if (rawLine==null) return (null);
		
    	String [] result = parser.parseLine(rawLine);
    	return (result);
    	
	}
	
	
	/**
	 * Get the next raw line from the buffered reader.  If there are no more lines, return null and close the buffered reader.  
	 * Subsequent calls will return null with no other effects.
	 * @return
	 * @throws IOException
	 */
	private String getNextRawLine()  {
		if (hasNext == false) return null;
		String nextLine=null;
		try {
			if (!this.linesSkiped) {
				for (int i = 0; i < linesToSkip; i++) {
					in.readLine();
				}
				this.linesSkiped = true;
				log.info("line skipping complete");
			}
			nextLine = in.readLine();
			if (nextLine == null) {
				hasNext = false;
				in.close();
				log.info("End of file reached, closing buffered reader.");
			}
		}
		catch (IOException e) {
			throw new ModularFileParserException("Can't read from file: ");
		}
        return nextLine;
    }

    /**
     * Closes the underlying reader.
     * 
     * @throws IOException if the close fails
     */
    public void close() {
    	try {
			in.close();
		} catch (IOException e) {
			throw new ModularFileParserException("Can't close file after reading finished");
		}
    }


	
    
	
}
