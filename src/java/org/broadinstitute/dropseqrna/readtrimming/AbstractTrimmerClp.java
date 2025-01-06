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
package org.broadinstitute.dropseqrna.readtrimming;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Histogram;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

public abstract class AbstractTrimmerClp extends CommandLineProgram {
    static final int UNPAIRED = 0;
    static final int FIRST_OF_PAIR = 1;
    static final int SECOND_OF_PAIR = 2;
    protected static final Set<Integer> VALID_WHICH_READ = CollectionUtil.makeSet(UNPAIRED, FIRST_OF_PAIR, SECOND_OF_PAIR);
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.")
    public File INPUT;
    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM file")
    public File OUTPUT;
    @Argument(doc = "The output summary statistics", optional = true)
    public File OUTPUT_SUMMARY;
    @Argument(doc = "Which reads to trim.  0: unpaired reads; 1: first of pair; 2: second of pair")
    public List<Integer> WHICH_READ = new ArrayList<>(Arrays.asList(0));
    protected int readsTrimmed = 0;
    protected int numReadsTotal = 0;
    protected final Histogram<Integer> numBasesTrimmed = new Histogram<>();

    protected boolean shouldTrim(final SAMRecord r) {
        if (WHICH_READ.contains(UNPAIRED) && r.getReadPairedFlag() == false) {
            return true;
        }
        if (WHICH_READ.contains(FIRST_OF_PAIR) && r.getFirstOfPairFlag()) {
            return true;
        }
        if (WHICH_READ.contains(SECOND_OF_PAIR) && r.getSecondOfPairFlag()) {
            return true;
        }
        return false;
    }

    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> list = new ArrayList<>(1);
        if (!VALID_WHICH_READ.containsAll(WHICH_READ)) {
            list.add("WHICH_READ must be one of " + VALID_WHICH_READ);
        }
        if (WHICH_READ.isEmpty()) {
            list.add("WHICH_READ must be specified");
        }
        return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
    }
}
