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
package org.broadinstitute.dropseqrna.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import picard.PicardException;

/**
 * Hold a buffered file writer and handle exceptions.
 * @author nemesh
 *
 */
public class OutputWriterUtil {

	public static BufferedWriter getWriter (File output) {
		BufferedWriter writer = null;
		try {
            writer = new BufferedWriter(new FileWriter(output));
        }
        catch (IOException ioe) {
            throw new PicardException("Error creating BufferedWriter " + output.getAbsolutePath() +
                    ": " + ioe.getMessage(), ioe);
        }
        return (writer);
	}
	
	public static void closeWriter (BufferedWriter writer) {
		
		try {
			writer.close();
		}
        catch (IOException ioe) {
            throw new PicardException("Error closing BufferedWriter "+
                    ": " + ioe.getMessage(), ioe);
        }
	}
	
	public static void writeResult (String result, BufferedWriter out) {
		try {
			out.write(result+"\n");
		} catch (IOException ioe) {
			throw new PicardException("Error writing to result file: " + ioe.getMessage(), ioe);
		}
	}
}
