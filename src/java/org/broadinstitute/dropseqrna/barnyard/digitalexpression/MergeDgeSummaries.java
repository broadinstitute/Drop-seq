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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.DigitalExpression;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

@CommandLineProgramProperties(
        summary = "Merges multiple Digital Gene Expression summary files into a single file.",
        oneLineSummary = "Merge Digital Gene Expression summary files.",
        programGroup = DropSeq.class
)
public class MergeDgeSummaries extends CommandLineProgram {

	private final Log log = Log.getInstance(MergeDgeSummaries.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The DGE files to be merged.", minElements = 1)
	public List<File> INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The merged DGE file.")
	public File OUTPUT;

    @Argument (doc="If true, sum metrics for same CBC.  Use when merging gene and metagene DGE summaries.  If false, " +
            "fail if a CBC appears more than once.")
    public boolean ACCUMULATE_CELL_BARCODE_METRICS = false;

	@Override
	protected int doWork() {
        IOUtil.assertFileIsWritable(OUTPUT);
        INPUT = FileListParsingUtils.expandFileList(INPUT);
        log.info("Merging [", INPUT.size()+"] DGE summary files");
        MetricsFile<DigitalExpression.DESummary, Integer> mergedMetricsFile = getMetricsFile();

        Map<String, DigitalExpression.DESummary> dgeSummaryMap = new HashMap<>();
        for (File summaryFile : INPUT) {
            final MetricsFile<DigitalExpression.DESummary, Integer> metricsFile = new MetricsFile<>();
            metricsFile.read(IOUtil.openFileForBufferedReading(summaryFile));
            List<DigitalExpression.DESummary> summaryList = metricsFile.getMetrics();
            for (DigitalExpression.DESummary dgeSummary : summaryList) {
                final DigitalExpression.DESummary existingSummary = dgeSummaryMap.get(dgeSummary.CELL_BARCODE);
                if (existingSummary == null) {
                    dgeSummaryMap.put(dgeSummary.CELL_BARCODE, dgeSummary);
                } else if (ACCUMULATE_CELL_BARCODE_METRICS) {
                    existingSummary.accumulate(dgeSummary);
                } else {
                    throw new RuntimeException("Cell barcode " + dgeSummary.CELL_BARCODE + " encountered multiple times");
                }
            }

            for (Header header : metricsFile.getHeaders()) {
                mergedMetricsFile.addHeader(header);
            }
        }

        DigitalExpression.writeSummary(dgeSummaryMap.values(), mergedMetricsFile, OUTPUT);

        return 0;
	}
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new MergeDgeSummaries().instanceMain(args));
	}
}
