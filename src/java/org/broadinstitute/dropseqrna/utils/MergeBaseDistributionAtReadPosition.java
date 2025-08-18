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
package org.broadinstitute.dropseqrna.utils;


import java.io.File;
import java.util.List;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import htsjdk.samtools.util.IOUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

/**
 * @author skashin
 *
 */
@CommandLineProgramProperties(summary = "Merges base composition per-position matrices",
        oneLineSummary = "Merges base composition per-position matrices",
        programGroup = DropSeq.class)
public class MergeBaseDistributionAtReadPosition extends CommandLineProgram {

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input histogram reports", minElements = 1)
	public List<File> INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output merged histogram report file.")
	public File OUTPUT;

	@Override
	protected int doWork() {
		IOUtil.assertFileIsWritable(OUTPUT);

        BaseDistributionMetricCollection mergedMetricCollection = null;
        for (final File reportFile : INPUT) {
            BaseDistributionMetricCollection metricCollection = BaseDistributionMetricCollection.readBaseDistribution(reportFile);
            if (mergedMetricCollection == null) {
                mergedMetricCollection = metricCollection;
            } else {
                mergedMetricCollection.mergeMetricCollections(metricCollection);
            }
        }
        mergedMetricCollection.writeOutput(OUTPUT);

		return 0;
	}
}
