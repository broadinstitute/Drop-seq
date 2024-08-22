/*
 * MIT License
 *
 * Copyright 2022 Broad Institute
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
package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@CommandLineProgramProperties(summary = "Perform edit distance UMI collapse, CBC filtering and chimeric read filtering " +
        "on a chimeric transcripts report (columns CELL_BARCODE, GENE, MOLECULAR_BARCODE, NUM_OBS, CHIMERIC) " +
        "produced by MarkChimericReads in order to create a molBC.txt file (columns CELL_BARCODE, GENE, MOLECULAR_BARCODE, NUM_OBS).",
        oneLineSummary = "Create a molBC file from a chimeric_transcripts file.",
        programGroup = DropSeq.class)
public class ChimericReportEditDistanceCollapse
        extends CommandLineProgram {
    private static final Log log = Log.getInstance(ChimericReportEditDistanceCollapse.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "chimeric_transcript file produced by MarkChimericReads. (columns CELL_BARCODE, GENE, MOLECULAR_BARCODE, NUM_OBS, CHIMERIC)")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "molBC file (columns CELL_BARCODE, GENE, MOLECULAR_BARCODE, NUM_OBS)")
    public File OUTPUT;

    @Argument(doc="Edit distance for UMI collapse")
    public int EDIT_DISTANCE = 1;

    @Argument(doc="List of cell barcodes to include.  Default: include all cell barcodes.", optional = true)
    public File CELL_BC_FILE;

    @Argument(doc="If true, skip lines in INPUT with CHIMERIC=true")
    public boolean IGNORE_CHIMERIC = true;

    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> list = new ArrayList<>(1);
        if (EDIT_DISTANCE < 0) {
            list.add("EDIT_DISTANCE must be non-negative.");
        }
        return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        final Set<String> selectedCellBarcodes;
        if (CELL_BC_FILE == null) {
            selectedCellBarcodes = null;
        } else {
            IOUtil.assertFileIsReadable(CELL_BC_FILE);
            selectedCellBarcodes = new HashSet<>( ParseBarcodeFile.readCellBarcodeFile(CELL_BC_FILE));
            log.info("Found " + selectedCellBarcodes.size()+ " cell barcodes in file " + CELL_BC_FILE);
        }

        BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT);
        GatherMolecularBarcodeDistributionByGene.writePerTranscriptHeader(out, false);

        final UMICollectionByCellParser parser = new UMICollectionByCellParser(INPUT, IGNORE_CHIMERIC);
        for (final List<UMICollection> umiCollections : parser) {
            if (selectedCellBarcodes != null && !selectedCellBarcodes.contains(umiCollections.getFirst().getCellBarcode())) {
                continue;
            }
            for (final UMICollection umiCollection : umiCollections) {
                GatherMolecularBarcodeDistributionByGene.writePerTranscriptStats(
                        umiCollection.getGeneName(),
                        umiCollection.getCellBarcode(),
                        umiCollection.getMolecularBarcodeCountsCollapsed(EDIT_DISTANCE),
                        out);
            }
        }
        try {
            out.close();
        } catch (IOException e) {
            throw new RuntimeIOException("Exception writing " + OUTPUT, e);
        }
        return 0;
    }
}
