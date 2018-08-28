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
package org.broadinstitute.dropseqrna.barnyard;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BarcodeListRetrievalTest {

    private static final String[] BARCODES = {
            "AGTGAGACAAGG",
            "AGTGCGACAAGG",
            "ACTGAGACAAGG",
            "GGTGAGACAAGG"
    };

    private static final File BAM = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.bam");
    private static final String CELL_BARCODE_TAG = "XC";
    private static final String MOLECULAR_BARCODE_TAG = "XM";
    private static final int READ_QUALITY = 10;

    private final BarcodeListRetrieval blr = new BarcodeListRetrieval();

    @Test
    public void testGetCellBarcodesFromFile() throws IOException {
        final File barcodeListFile = File.createTempFile("BarcodeListRetrievalTest.", ".barcode_list");
        barcodeListFile.deleteOnExit();
        final FileWriter writer = new FileWriter(barcodeListFile);
        for (final String barcode: BARCODES) {
            writer.write(barcode);
            writer.write("\n");
        }
        writer.close();
        final List<String> barcodes = blr.getCellBarcodes(null, null, null,
                null, null, null, null, null,
                barcodeListFile,
                null, null, null, null,
                null, null, null);
        Assert.assertEquals(barcodes, Arrays.asList(BARCODES));
    }

    @Test
    public void testGetCellBarcodesByGeneCount() {
        List<String> barcodes = blr.getCellBarcodes(BAM, CELL_BARCODE_TAG, MOLECULAR_BARCODE_TAG,
                GeneFunctionCommandLineBase.DEFAULT_GENE_NAME_TAG,
        null, null, null, null, null, READ_QUALITY,
        null, 3, null, null,
                null, null);
        Assert.assertEquals(barcodes.size(), 14);
        barcodes = blr.getCellBarcodes(BAM, CELL_BARCODE_TAG, MOLECULAR_BARCODE_TAG,
                GeneFunctionCommandLineBase.DEFAULT_GENE_NAME_TAG,
                null, null, null, null, null, READ_QUALITY,
                null, 6, null, null,
                null, null);
        Assert.assertEquals(barcodes.size(), 0);
    }

    @Test
    public void testGetCellBarcodesByTranscriptCount() {
        List<String> barcodes = blr.getCellBarcodes(BAM, CELL_BARCODE_TAG, MOLECULAR_BARCODE_TAG,
                GeneFunctionCommandLineBase.DEFAULT_GENE_NAME_TAG,
                GeneFunctionCommandLineBase.DEFAULT_GENE_STRAND_TAG, GeneFunctionCommandLineBase.DEFAULT_GENE_FUNCTION_TAG,
                GeneFunctionCommandLineBase.DEFAULT_STRAND_STRATEGY, GeneFunctionCommandLineBase.DEFAULT_LOCUS_FUNCTION_LIST,
                null, READ_QUALITY,
                10, null, 10, null,
                1, 10);
        Assert.assertEquals(barcodes.size(), 8);

        barcodes = blr.getCellBarcodes(BAM, CELL_BARCODE_TAG, MOLECULAR_BARCODE_TAG,
                GeneFunctionCommandLineBase.DEFAULT_GENE_NAME_TAG,
                GeneFunctionCommandLineBase.DEFAULT_GENE_STRAND_TAG, GeneFunctionCommandLineBase.DEFAULT_GENE_FUNCTION_TAG,
                GeneFunctionCommandLineBase.DEFAULT_STRAND_STRATEGY, GeneFunctionCommandLineBase.DEFAULT_LOCUS_FUNCTION_LIST,
                null, READ_QUALITY,
                200, null, 10, null,
                1, 10);
        Assert.assertEquals(barcodes.size(), 0);
    }

    @Test
    public void testGetCellBarcodesByReadCount() {
        List<String> barcodes = blr.getCellBarcodes(BAM, CELL_BARCODE_TAG, null, null,
                null, null, null, null,null, READ_QUALITY,
                null, null, 10, null,
                null, null);
        Assert.assertEquals(barcodes.size(), 16);
        barcodes = blr.getCellBarcodes(BAM, CELL_BARCODE_TAG, null, null,
                null, null, null, null,null, READ_QUALITY,
                null, null, 5000, null,
                null, null);
        Assert.assertEquals(barcodes.size(), 0);
    }

    @Test
    public void testGetCellBarcodesByNumCoreBarcodes() {
        List<String> barcodes = blr.getCellBarcodes(BAM, CELL_BARCODE_TAG, null, null,
                null, null, null, null,null, READ_QUALITY,
                null, null, null, 100,
                null, null);
        Assert.assertEquals(barcodes.size(), 59);
        barcodes = blr.getCellBarcodes(BAM, CELL_BARCODE_TAG, null, null,
                null, null, null, null,null, 1,
                null, null, null, 100,
                null, null);
        Assert.assertEquals(barcodes.size(), 64);
        barcodes = blr.getCellBarcodes(BAM, CELL_BARCODE_TAG, null, null,
                null, null, null, null,null, READ_QUALITY,
                null, null, null, 10,
                null, null);
        Assert.assertEquals(barcodes.size(), 10);
    }
}
