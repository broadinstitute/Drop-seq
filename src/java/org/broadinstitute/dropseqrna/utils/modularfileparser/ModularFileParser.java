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
