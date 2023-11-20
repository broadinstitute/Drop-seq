/*
 * MIT License
 *
 * Copyright 2023 Broad Institute
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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintWriter;
import picard.cmdline.CommandLineProgram;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

public abstract class AbstractTripletDgeWriterClp
extends CommandLineProgram {
    private final Log log = Log.getInstance(AbstractTripletDgeWriterClp.class);
    @Argument(shortName = "G", doc="Unheadered text file with one gene per line.  If filename ends with .gz, will be compressed.", optional = true)
    public File OUTPUT_GENES;
    @Argument(shortName = "C", doc="Unheadered text file with one cell barcode per line.  If filename ends with .gz, will be compressed.", optional = true)
    public File OUTPUT_CELLS;
    @Argument(shortName = "F", doc="Unheadered tab-separated text file with each line containing {gene,gene,feature-type}.  If filename ends with .gz, will be compressed.", optional = true)
    public File OUTPUT_FEATURES;
    @Argument(doc = "When writing OUTPUT_FEATURES, the value to put in the 3rd column")
    public String FEATURE_TYPE = "Gene Expression";

    protected void writeListFile(final File outputFile, final List<String> values) {
        log.info("Writing ", outputFile);
        IOUtil.assertFileIsWritable(outputFile);
        ErrorCheckingPrintWriter out = new ErrorCheckingPrintWriter(IOUtil.openFileForWriting(outputFile));
        for (final String value: values) {
            out.println(value);
        }
        out.close();
    }

    protected void maybeWriteNamesFiles(final List<String> genes, final List<String> cellBarcodes) {
        if (OUTPUT_GENES != null) {
            writeListFile(OUTPUT_GENES, genes);
        }
        if (OUTPUT_FEATURES != null) {
            final List<String> features = genes.stream().map(
                    gene -> StringUtil.join("\t", gene, gene, FEATURE_TYPE)).collect(Collectors.toList());
            writeListFile(OUTPUT_FEATURES, features);
        }
        if (OUTPUT_CELLS != null) {
            writeListFile(OUTPUT_CELLS, cellBarcodes);
        }

    }
}
