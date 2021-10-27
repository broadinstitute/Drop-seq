/*
 * MIT License
 *
 * Copyright 2020 Broad Institute
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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.ReportFileUtil;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

/**
 * @author skashin
 */
@CommandLineProgramProperties(summary = "Merges cells to samples assignment files.",
        oneLineSummary = "Merges cells to samples assignment files.",
        programGroup = DropSeq.class)
public class MergeCellToSampleAssignments extends CommandLineProgram {

    private static final Log log = Log.getInstance(MergeCellToSampleAssignments.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input assignment files to be merged.", minElements = 1)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The merged assignment file.")
    public File OUTPUT;

    @Argument (doc="Set to true to emit a column that contains a population average likelihood score aggregated across all SNPs.  Setting this to false is purely for backwards compatability", optional=true)
	public Boolean EMIT_POPULATION_LIKELIHOOD=true;
	
    @Override
    protected int doWork() {
        IOUtil.assertFileIsWritable(OUTPUT);

        INPUT = FileListParsingUtils.expandFileList(INPUT);

        List < String > samples = null;
        CellCollectionSampleLikelihoodCollection mergedCollection = null;
        for (File assignmentFile : INPUT) {
            CellCollectionSampleLikelihoodCollection collection = CellCollectionSampleLikelihoodCollection.parseFile(assignmentFile);
            if (mergedCollection == null) {
                mergedCollection = collection;
                samples = collection.getSamples();
            } else {
                if (!collection.getSamples().equals(samples)) {
                    throw new RuntimeException("Make sure all assignment files contain the same list of samples");
                }

                for (String cellBarcode : collection.getCellBarcodes()) {
                    if (mergedCollection.getCellBarcodes().contains(cellBarcode)) {
                        throw new RuntimeException("Cell barcode " + cellBarcode + " is in multiple assignment files");
                    }
                    mergedCollection.add(collection.getLikelihoodCollection(cellBarcode));
                }
            }
        }

        ErrorCheckingPrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.OUTPUT));
        for (String comment : ReportFileUtil.readComments(INPUT))
            out.println(comment);
        for (Header header : getDefaultHeaders())
            out.println(ReportFileUtil.COMMENT_PREFIX + header);

        AssignCellsToSamples.writeBestLikelihoods(samples, mergedCollection, null, out, this.EMIT_POPULATION_LIKELIHOOD);
        out.close();

        return 0;
    }
}
