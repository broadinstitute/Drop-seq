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
package org.broadinstitute.dropseqrna.metagene;

import java.io.File;
import java.util.List;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.dropseqrna.metagene.DiscoverMetaGenes;
import org.broadinstitute.dropseqrna.metagene.UMIMetaGeneAggregation;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(summary = "Merges metagene reports.",
        oneLineSummary = "Merges metagene reports.",
        programGroup = DropSeq.class)
public class MergeMetaGeneReports extends CommandLineProgram {

	private static final Log log = Log.getInstance(MergeMetaGeneReports.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input metagene reports to be merged.", minElements = 1)
	public List<File> INPUT;

    @Argument(doc="Should single gene UMI counts be written out as separate lines into the OUTPUT")
    public boolean WRITE_SINGLE_GENES = false;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The merged metagene report.")
	public File OUTPUT;

	@Override
	protected int doWork() {
        IOUtil.assertFileIsWritable(OUTPUT);

        UMIMetaGeneAggregation mergedAggregator = null;
        for (File reportFile : INPUT) {
            UMIMetaGeneAggregation aggregator = DiscoverMetaGenes.readReport(reportFile);
            if (mergedAggregator == null) {
                mergedAggregator = aggregator;
            } else {
                mergedAggregator.add(aggregator);
            }
        }

        DiscoverMetaGenes.writeReport(mergedAggregator, OUTPUT, WRITE_SINGLE_GENES);

        return 0;
	}
}
