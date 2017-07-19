/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2015 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.util.IOUtil;
import picard.PicardException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Utility to write Fasta files.
 * Copied from Picard-private in order to eliminate dependence on picard-private.jar
 *
 * @author ktibbett@broadinstitute.org
 */
public class FastaSequenceFileWriter {

    private int lineLength = 50;    // Length of lines of sequence in the fasta file
    private BufferedWriter writer = null;

    /**
     * Constructor that uses the default line length.  Checks that the file is
     * writeable and creates a BufferedWriter to do the writing
     */
    public FastaSequenceFileWriter(File fastaFile) {
        IOUtil.assertFileIsWritable(fastaFile);
        try {
            writer = new BufferedWriter(new FileWriter(fastaFile));
        }
        catch (IOException ioe) {
            throw new PicardException("Error creating BufferedWriter " + fastaFile.getAbsolutePath() +
                    ": " + ioe.getMessage(), ioe);
        }
    }

    /**
     * Constructor that uses a user-provided line length
     */
    public FastaSequenceFileWriter(File fastaFile, int lineLength) {
        this(fastaFile);
        this.lineLength = lineLength;
    }


    /**
     * Writes a sequence to the file.  Prefaces the name with ">" and writes
     * the sequence itself in lines whose length is specified by <code>lineLength</code>
     */
    public void writeSequence(String name, String sequence) {
        try {
            writer.write(">" + name);
            writer.newLine();
            int startPos = 0;
            do {
                int endPos = Math.min(startPos + lineLength, sequence.length());
                writer.write(sequence.substring(startPos, endPos));
                writer.newLine();
                startPos += lineLength;
            }
            while (startPos < sequence.length());
        }
        catch (IOException ioe) {
            throw new PicardException("Error writing to fasta file: " + ioe.getMessage(), ioe);
        }
    }


    /**
     * Closes the BufferedWriter
     */
    public void close() {
        if (writer != null) {
            try {
                writer.close();
            }
            catch (IOException ioe) {
                throw new RuntimeException(ioe);
            }
        }
    }
}
