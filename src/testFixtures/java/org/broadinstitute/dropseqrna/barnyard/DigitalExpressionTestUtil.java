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
package org.broadinstitute.dropseqrna.barnyard;

import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintWriter;
import org.testng.Assert;

import java.io.File;
import java.io.IOException;
import java.util.Random;

public class DigitalExpressionTestUtil {
    public static final String DIGITAL_EXPRESSION_EXTENSION = ".digital_expression.txt.gz";
    public static final String DIGITAL_EXPRESSION_SUMMARY_EXTENSION = ".digital_expression_summary.txt";
    /**
     * For the integration test, make a relatively small BAM file and read it in.
     * Manually test the counts of the reads across UMIs - umi testing to collapse molecular barcodes is in its own test.
     *
     * Data taken from 9-27-14 100 cells data set.
     *
     * cells
     * ATCAGGGACAGA
     * AGGGAAAATTGA
     * TTGCCTTACGCG
     * TGGCGAAGAGAT
     * TACAATTAAGGC
     *
     * genes
     * HUMAN_10:101948055-101989376:CHUK
     * HUMAN_15:101821715-101835487:SNRPA1
     * HUMAN_3:42642106-42690227:NKTR
     *
     * /fg/software/gap/gap_analysis/FilterBamByTag I=100cells_star_bq10_noPseudoGenes.bam O=test.bam TAG=ZC TAG_VALUES_FILE=bc.txt ACCEPT_TAG=true
     * samtools view -H test.bam > 5cell3gene.sam
     * samtools view test.bam |grep -f genes.txt >> 5cell3gene.sam
     * samtools view -Sb 5cell3gene.sam > 5cell3gene.bam
     *
     */

    static final File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.bam");
    static final String [] barcodes ={"ATCAGGGACAGA", "AGGGAAAATTGA", "TTGCCTTACGCG", "TGGCGAAGAGAT", "TACAATTAAGGC"};

    /**
     * @return A digital expression file, which is marked as deleteOnExit
*/
    public static File makeDigitalExpressionFile(final boolean outputHeader) throws IOException {
        return makeDigitalExpressionFile(outputHeader, new File("."));
    }

    /** Use this version when calling from private tests */
public static File makeDigitalExpressionFile(final boolean outputHeader, final File basedir) throws IOException {
        final File outFile = File.createTempFile("testDigitalExpression.", DIGITAL_EXPRESSION_EXTENSION);
        final File summaryFile = makeSummaryPathFromDgePath(outFile);
final File cellBarcodesFile = File.createTempFile("testDigitalExpression.", ".selectedCellBarcodes.txt");
        outFile.deleteOnExit();
        summaryFile.deleteOnExit();
cellBarcodesFile.deleteOnExit();
final ErrorCheckingPrintWriter writer = new ErrorCheckingPrintWriter(cellBarcodesFile);
for (final String cellBarcode : barcodes)
            writer.println(cellBarcode);
writer.close();
final DigitalExpression de = new DigitalExpression();
        de.INPUT = new File(basedir, IN_FILE.getPath());
        de.OUTPUT = outFile;
        de.SUMMARY = summaryFile;
de.CELL_BC_FILE = cellBarcodesFile;
de.OUTPUT_HEADER = outputHeader;
de.UNIQUE_EXPERIMENT_ID = "UIE" + new Random().nextInt();

Assert.assertEquals(de.doWork(), 0);
return outFile;
    }

    public static File makeSummaryPathFromDgePath(final File dge) {
        String dgeName = dge.getName();
        String basename = dgeName.substring(0, dgeName.length() - DIGITAL_EXPRESSION_EXTENSION.length());
        return new File(dge.getParentFile(),
                basename + DIGITAL_EXPRESSION_SUMMARY_EXTENSION);
    }
}
