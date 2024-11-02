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
package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.util.IntervalList;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class CreateIntervalsFilesTest {

    private static final File METADATA_DIR = new File("testdata/org/broadinstitute/transcriptome/annotation");
    private static final String REFERENCE_NAME = "mm10";
    private static final List<String> MT_SEQUENCE = List.of("MT");
    private static final List<String> NON_AUTOSOME_SEQUENCE = Arrays.asList("X", "Y", "MT");

    @Test
    public void testCreateIntervalsFiles() throws IOException {
        final CreateIntervalsFiles clp = new CreateIntervalsFiles();
        clp.SEQUENCE_DICTIONARY = new File(METADATA_DIR, REFERENCE_NAME + ".dict");
        clp.REDUCED_GTF = new File(METADATA_DIR, REFERENCE_NAME + ".reduced.gtf.gz");
        clp.OUTPUT = TestUtils.createTempDirectory("CreateIntervalsFilesTest");
        clp.PREFIX = REFERENCE_NAME;
        clp.MT_SEQUENCE = MT_SEQUENCE;
        clp.NON_AUTOSOME_SEQUENCE = NON_AUTOSOME_SEQUENCE;

        try {
            Assert.assertEquals(clp.doWork(), 0);

            final IntervalList mtIntervals = IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".mt.intervals"));
            Assert.assertEquals(mtIntervals.size(), 37);

            final IntervalList nonAutosomeSequences = IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".non_autosomes.intervals"));
            Assert.assertEquals(nonAutosomeSequences.size(), 0);
            Assert.assertEquals(nonAutosomeSequences.getHeader().getSequenceDictionary().size(), 3);

            final IntervalList rRnaIntervals = IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".rRNA.intervals"));
            Assert.assertEquals(rRnaIntervals.size(), 355);

            // This is hard to validate, so just confirm its readability.
            IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".intergenic.intervals"));

            final IntervalList genesIntervals = IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".genes.intervals"));
            Assert.assertEquals(genesIntervals.size(), 32976);

            final IntervalList exonsIntervals = IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".exons.intervals"));
            Assert.assertEquals(exonsIntervals.size(), 615275);

            final IntervalList consensusIntronsIntervals = IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".consensus_introns.intervals"));
            Assert.assertEquals(consensusIntronsIntervals.size(), 396038);

        } finally {
            FileUtils.deleteDirectory(clp.OUTPUT);
        }
    }

    @Test
    public void testCreateIntervalsFilesNoMT() throws IOException {
        final CreateIntervalsFiles clp = new CreateIntervalsFiles();
        clp.SEQUENCE_DICTIONARY = new File(METADATA_DIR, REFERENCE_NAME + ".dict");
        clp.REDUCED_GTF = new File(METADATA_DIR, REFERENCE_NAME + ".reduced.gtf.gz");
        clp.OUTPUT = TestUtils.createTempDirectory("CreateIntervalsFilesTest");
        clp.PREFIX = REFERENCE_NAME;
        clp.NON_AUTOSOME_SEQUENCE = NON_AUTOSOME_SEQUENCE;

        try {
            Assert.assertEquals(clp.doWork(), 0);

            Assert.assertFalse(new File(clp.OUTPUT, REFERENCE_NAME + ".mt.intervals").exists());

            final IntervalList nonAutosomeSequences =
                    IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".non_autosomes.intervals"));
            Assert.assertEquals(nonAutosomeSequences.size(), 0);
            Assert.assertEquals(nonAutosomeSequences.getHeader().getSequenceDictionary().size(), 3);

            final IntervalList rRnaIntervals =
                    IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".rRNA.intervals"));
            Assert.assertEquals(rRnaIntervals.size(), 355);

            // This is hard to validate, so just confirm its readability.
            IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".intergenic.intervals"));

            final IntervalList genesIntervals =
                    IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".genes.intervals"));
            Assert.assertEquals(genesIntervals.size(), 32976);

            final IntervalList exonsIntervals =
                    IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".exons.intervals"));
            Assert.assertEquals(exonsIntervals.size(), 615275);

            final IntervalList consensusIntronsIntervals =
                    IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".consensus_introns.intervals"));
            Assert.assertEquals(consensusIntronsIntervals.size(), 396038);

        } finally {
            FileUtils.deleteDirectory(clp.OUTPUT);
        }
    }

    @Test
    public void testCreateIntervalsFilesNoNonAutosome() throws IOException {
        final CreateIntervalsFiles clp = new CreateIntervalsFiles();
        clp.SEQUENCE_DICTIONARY = new File(METADATA_DIR, REFERENCE_NAME + ".dict");
        clp.REDUCED_GTF = new File(METADATA_DIR, REFERENCE_NAME + ".reduced.gtf.gz");
        clp.OUTPUT = TestUtils.createTempDirectory("CreateIntervalsFilesTest");
        clp.PREFIX = REFERENCE_NAME;
        clp.MT_SEQUENCE = MT_SEQUENCE;

        try {
            Assert.assertEquals(clp.doWork(), 0);

            final IntervalList mtIntervals =
                    IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".mt.intervals"));
            Assert.assertEquals(mtIntervals.size(), 37);

            Assert.assertFalse(new File(clp.OUTPUT, REFERENCE_NAME + ".non_autosomes.intervals").exists());

            final IntervalList rRnaIntervals =
                    IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".rRNA.intervals"));
            Assert.assertEquals(rRnaIntervals.size(), 355);

            // This is hard to validate, so just confirm its readability.
            IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".intergenic.intervals"));

            final IntervalList genesIntervals =
                    IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".genes.intervals"));
            Assert.assertEquals(genesIntervals.size(), 32976);

            final IntervalList exonsIntervals =
                    IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".exons.intervals"));
            Assert.assertEquals(exonsIntervals.size(), 615275);

            final IntervalList consensusIntronsIntervals =
                    IntervalList.fromFile(new File(clp.OUTPUT, REFERENCE_NAME + ".consensus_introns.intervals"));
            Assert.assertEquals(consensusIntronsIntervals.size(), 396038);

        } finally {
            FileUtils.deleteDirectory(clp.OUTPUT);
        }
    }
}
