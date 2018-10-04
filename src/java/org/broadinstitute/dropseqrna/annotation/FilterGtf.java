/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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
package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.MetaData;
import org.broadinstitute.dropseqrna.utils.DropSeqSamUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@CommandLineProgramProperties(
        summary = "Remove GTF records with undesired gene_biotypes from GTF",
        oneLineSummary = "Filter a GTF",
        programGroup = MetaData.class
)

public class FilterGtf extends CommandLineProgram {
    @Argument(doc="Input GTF file to be filtered.")
    public File GTF;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output filtered GTF file")
    public File OUTPUT;

    @Argument(shortName = "G", doc="gene_biotype value that flags a GTF record as undesired")
    public List<String> UNDESIRED_GENE_BIOTYPE;

    @Argument(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, optional = true,
            doc="The reference sequence dictionary. If specified, GTF records for sequences not in the dictionary will be discarded.")
    public File SEQUENCE_DICTIONARY;


    private Set<String> geneBiotypesToFilter = new HashSet<>();
    private SAMSequenceDictionary sequenceDictionary;

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(this.GTF);
        IOUtil.assertFileIsWritable(this.OUTPUT);

        geneBiotypesToFilter.addAll(UNDESIRED_GENE_BIOTYPE);
        if (SEQUENCE_DICTIONARY != null) {
            sequenceDictionary = DropSeqSamUtil.loadSequenceDictionary(SEQUENCE_DICTIONARY);
        }

        // Don't use TabbedInputParser because we want to preserve comments
        // Don't use GTFReader because 1) we want to preserve comments; and 2) what is needed at this point is not
        // complicated and custom parsing will be more straightforward.
        final BufferedReader reader = IOUtil.openFileForBufferedReading(GTF);
        final BufferedWriter writer = IOUtil.openFileForBufferedWriting(OUTPUT);
        String inputLine;
        try {
            while ((inputLine = reader.readLine()) != null ) {
                if (inputLine.startsWith("#") || goodLine(inputLine)) {
                    writer.write(inputLine);
                    writer.newLine();
                }
            }
            CloserUtil.close(reader);
            writer.close();
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
        return 0;
    }

    /**
     * See if GTF record contains one of the undesired gene_biotypes
     * @param inputLine
     * @return true if the GTF record has no gene_biotype, or the gene_biotype is not one of the undesired ones.
     */
    private boolean goodLine(String inputLine) {
        final String[] fields = inputLine.split("\t");
        if (fields.length != 9) {
            throw new IllegalArgumentException("GTF line does not have 9 fields: " + inputLine);
        }
        if (sequenceDictionary != null && sequenceDictionary.getSequence(fields[0]) == null) {
            return false;
        }
        final String attributesString = fields[8].replaceAll("\"", "");
        for (final String attribute : attributesString.split(" *; *")) {
            final String[] keyValue = attribute.split(" ", 2);
            if (keyValue.length < 2) {
                throw new IllegalArgumentException("GTF attribute does not have space: " + inputLine);
            }
            if (keyValue[0].equals("gene_biotype") && geneBiotypesToFilter.contains(keyValue[1])) {
                return false;
            }
        }
        return true;
    }

    public static void main(final String[]args) {
        new FilterGtf().instanceMainWithExit(args);
    }
}
