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

import htsjdk.samtools.reference.FastaReferenceWriter;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.IOException;

/**
 * Thin wrapper around htsjdk.samtools.reference.FastaReferenceWriter
 */
public class FastaSequenceFileWriter {

    private static final int DEFAULT_LINE_LENGTH = 50;
    private final FastaReferenceWriter writer;

    public FastaSequenceFileWriter(File fastaFile) {
        this(fastaFile, DEFAULT_LINE_LENGTH);
    }

    public FastaSequenceFileWriter(File fastaFile, int lineLength) {
        IOUtil.assertFileIsWritable(fastaFile);
        try {
            writer = new FastaReferenceWriterBuilder().setBasesPerLine(lineLength).setFastaFile(fastaFile.toPath()).build();
        }
        catch (IOException ioe) {
            throw new RuntimeIOException("Error creating FastaReferenceWriter " + fastaFile.getAbsolutePath(), ioe);
        }
    }


    public void writeSequence(String name, String sequence) {
        try {
            writer.startSequence(name).appendBases(sequence);
        }
        catch (IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }

    public void writeSequence(String name, String sequence, String description) {
        try {
            writer.startSequence(name, description).appendBases(sequence);
        }
        catch (IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }


    public void close() {
        try {
            writer.close();
        }
        catch (IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }
}
