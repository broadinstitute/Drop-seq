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
