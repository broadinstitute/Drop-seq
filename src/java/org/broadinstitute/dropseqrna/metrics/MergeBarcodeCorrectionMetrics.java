/*
 * MIT License
 *
 * Copyright 2025 Broad Institute
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
package org.broadinstitute.dropseqrna.metrics;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.beadsynthesis.BarcodeCorrectionMetrics;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.List;

@CommandLineProgramProperties(summary = "Merges BarcodeCorrectionMetrics files, produced by CorrectAndSplitScrnaReadPairs.",
        oneLineSummary = "Merges BarcodeCorrectionMetrics files.",
        programGroup = DropSeq.class)
public class MergeBarcodeCorrectionMetrics
extends CommandLineProgram {
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input read quality metrics files to be merged.", minElements = 1)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The merged read quality metrics file.")
    public File OUTPUT;

    @Argument(shortName = "D", doc="Delete input files after merging.")
    public Boolean DELETE_INPUTS;


    @Override
    protected int doWork() {
        new MergeMetricsHelper<BarcodeCorrectionMetrics, Integer>().
                mergeSingleMetric(OUTPUT, INPUT, DELETE_INPUTS, BarcodeCorrectionMetrics::new, BarcodeCorrectionMetrics::merge,
                        this.getMetricsFile());
        return 0;
    }
}
