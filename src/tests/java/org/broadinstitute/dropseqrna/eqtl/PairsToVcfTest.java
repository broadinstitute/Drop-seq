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
package org.broadinstitute.dropseqrna.eqtl;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class PairsToVcfTest {

    private static final File TEST_DATA_DIR =
            new File("testdata/org/broadinstitute/dropseq/eqtl/PairsToVcfTest");

    private static final File INPUT_PAIRS =
            new File(TEST_DATA_DIR, "input_pairs.tsv");

    private static final File MAP_COLUMNS =
            new File(TEST_DATA_DIR, "map_columns.tsv");

    private static final File EXPECTED_SORTED_PAIRS =
            new File(TEST_DATA_DIR, "expected_sorted_pairs.vcf");

    private static final File EXPECTED_UNSORTED_PAIRS =
            new File(TEST_DATA_DIR, "expected_unsorted_pairs.vcf");

    private static final File EXPECTED_UNSWAPPED_PAIRS =
            new File(TEST_DATA_DIR, "expected_unswapped_pairs.vcf");

    private static final File EXPECTED_MAPPED_PAIRS =
            new File(TEST_DATA_DIR, "expected_mapped_pairs.vcf");

    private static final File VCF =
            new File(TEST_DATA_DIR, "genotypes.vcf.gz");

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testStrictFailure() throws Exception {
        final PairsToVcf pairsToVcf = new PairsToVcf();
        final File outputFile = File.createTempFile("PairsToVcfTest.", ".vcf.gz");
        outputFile.deleteOnExit();
        final String[] argv = {
                "INPUT=" + INPUT_PAIRS.getAbsolutePath(),
                "VCF=" + VCF.getAbsolutePath(),
                "VARIANT_COLUMN=variant_id",
                "BETA_COLUMN=slope",
                "VALIDATION_STRINGENCY=STRICT",
                "OUTPUT=" + outputFile.getAbsolutePath(),
        };
        pairsToVcf.instanceMain(argv);
    }

    @Test
    public void testNoBeta() throws Exception {
        final PairsToVcf pairsToVcf = new PairsToVcf();
        final File outputFile = File.createTempFile("PairsToVcfTest.", ".vcf.gz");
        outputFile.deleteOnExit();
        final String[] argv = {
                "INPUT=" + INPUT_PAIRS.getAbsolutePath(),
                "VCF=" + VCF.getAbsolutePath(),
                "VARIANT_COLUMN=variant_id",
                "OUTPUT=" + outputFile.getAbsolutePath(),
        };
        Assert.assertEquals(pairsToVcf.instanceMain(argv), 0);
        testVcfs(outputFile, EXPECTED_UNSWAPPED_PAIRS);
    }

    @Test
    public void testSingleThread() throws Exception {
        final PairsToVcf pairsToVcf = new PairsToVcf();
        final File outputFile = File.createTempFile("PairsToVcfTest.", ".vcf.gz");
        outputFile.deleteOnExit();
        final String[] argv = {
                "INPUT=" + INPUT_PAIRS.getAbsolutePath(),
                "VCF=" + VCF.getAbsolutePath(),
                "VARIANT_COLUMN=variant_id",
                "BETA_COLUMN=slope",
                "VALIDATION_STRINGENCY=SILENT",
                "OUTPUT=" + outputFile.getAbsolutePath(),
        };
        Assert.assertEquals(pairsToVcf.instanceMain(argv), 0);
        testVcfs(outputFile, EXPECTED_SORTED_PAIRS);
    }

    @Test
    public void testNoCommandLine() throws Exception {
        final PairsToVcf pairsToVcf = new PairsToVcf();
        final File outputFile = File.createTempFile("PairsToVcfTest.", ".vcf.gz");
        outputFile.deleteOnExit();
        final String[] argv = {
                "INPUT=" + INPUT_PAIRS.getAbsolutePath(),
                "VCF=" + VCF.getAbsolutePath(),
                "VARIANT_COLUMN=variant_id",
                "BETA_COLUMN=slope",
                "VALIDATION_STRINGENCY=SILENT",
                "OUTPUT_COMMANDLINE=false",
                "OUTPUT=" + outputFile.getAbsolutePath(),
        };
        Assert.assertEquals(pairsToVcf.instanceMain(argv), 0);
        testVcfs(outputFile, EXPECTED_SORTED_PAIRS);
        try (final VCFReader reader = new VCFFileReader(outputFile)) {
            final VCFHeader header = reader.getHeader();
            final VCFHeaderLine commandLine = header.getMetaDataLine("PairsToVcfCommandLine");
            Assert.assertNull(commandLine);
        }
    }

    @Test
    public void testCommandLine() throws Exception {
        final PairsToVcf pairsToVcf = new PairsToVcf();
        final File outputFile = File.createTempFile("PairsToVcfTest.", ".vcf.gz");
        outputFile.deleteOnExit();
        final String[] argv = {
                "INPUT=" + INPUT_PAIRS.getAbsolutePath(),
                "VCF=" + VCF.getAbsolutePath(),
                "VARIANT_COLUMN=variant_id",
                "BETA_COLUMN=slope",
                "VALIDATION_STRINGENCY=SILENT",
                "OUTPUT_COMMANDLINE=true",
                "OUTPUT=" + outputFile.getAbsolutePath(),
        };
        Assert.assertEquals(pairsToVcf.instanceMain(argv), 0);
        testVcfs(outputFile, EXPECTED_SORTED_PAIRS);
        try (final VCFReader reader = new VCFFileReader(outputFile)) {
            final VCFHeader header = reader.getHeader();
            final VCFHeaderLine commandLine = header.getMetaDataLine("PairsToVcfCommandLine");
            Assert.assertNotNull(commandLine);
        }
    }

    @Test
    public void testParallel() throws Exception {
        final PairsToVcf pairsToVcf = new PairsToVcf();
        final File outputFile = File.createTempFile("PairsToVcfTest.", ".vcf.gz");
        outputFile.deleteOnExit();
        final String[] argv = {
                "INPUT=" + INPUT_PAIRS.getAbsolutePath(),
                "VCF=" + VCF.getAbsolutePath(),
                "VARIANT_COLUMN=variant_id",
                "BETA_COLUMN=slope",
                "NUM_THREADS=2",
                "VALIDATION_STRINGENCY=SILENT",
                "OUTPUT=" + outputFile.getAbsolutePath(),
        };
        Assert.assertEquals(pairsToVcf.instanceMain(argv), 0);
        testVcfs(outputFile, EXPECTED_SORTED_PAIRS);
    }

    @Test
    public void testUnsorted() throws Exception {
        final PairsToVcf pairsToVcf = new PairsToVcf();
        final File outputFile = File.createTempFile("PairsToVcfTest.", ".vcf.gz");
        outputFile.deleteOnExit();
        final String[] argv = {
                "INPUT=" + INPUT_PAIRS.getAbsolutePath(),
                "VCF=" + VCF.getAbsolutePath(),
                "VARIANT_COLUMN=variant_id",
                "BETA_COLUMN=slope",
                "NUM_THREADS=2",
                "SORT_OUTPUT=false",
                "VALIDATION_STRINGENCY=SILENT",
                "CREATE_INDEX=false",
                "OUTPUT=" + outputFile.getAbsolutePath(),
        };
        Assert.assertEquals(pairsToVcf.instanceMain(argv), 0);
        testVcfs(outputFile, EXPECTED_UNSORTED_PAIRS);
    }

    @Test
    public void testMapped() throws Exception {
        final PairsToVcf pairsToVcf = new PairsToVcf();
        final File outputFile = File.createTempFile("PairsToVcfTest.", ".vcf.gz");
        outputFile.deleteOnExit();
        final String[] argv = {
                "INPUT=" + INPUT_PAIRS.getAbsolutePath(),
                "VCF=" + VCF.getAbsolutePath(),
                "VARIANT_COLUMN=variant_id",
                "BETA_COLUMN=slope",
                "NUM_THREADS=2",
                "VALIDATION_STRINGENCY=SILENT",
                "MAP_COLUMNS=" + MAP_COLUMNS.getAbsolutePath(),
                "OUTPUT=" + outputFile.getAbsolutePath(),
        };
        Assert.assertEquals(pairsToVcf.instanceMain(argv), 0);
        testVcfs(outputFile, EXPECTED_MAPPED_PAIRS);
    }

    private static void testVcfs(final File actualVcf, final File expectedVcf) throws Exception {
        try (
                final VCFReader actualReader = new VCFFileReader(actualVcf, false);
                final VCFReader expectedReader = new VCFFileReader(expectedVcf, false)
        ) {
            final VCFHeader actualHeader = actualReader.getHeader();
            final VCFHeader expectedHeader = expectedReader.getHeader();
            Assert.assertEquals(actualHeader.getNGenotypeSamples(), 0);
            Assert.assertEquals(actualHeader.getContigLines(), expectedHeader.getContigLines());
            Assert.assertEquals(actualHeader.getIDHeaderLines(), expectedHeader.getIDHeaderLines());
            Assert.assertEquals(actualHeader.getFilterLines(), expectedHeader.getFilterLines());

            final CloseableIterator<VariantContext> actualIterator = actualReader.iterator();
            final CloseableIterator<VariantContext> expectedIterator = expectedReader.iterator();
            while (actualIterator.hasNext() && expectedIterator.hasNext()) {
                final VariantContext actualVariant = actualIterator.next();
                final VariantContext expectedVariant = expectedIterator.next();
                Assert.assertEquals(actualVariant.getContig(), expectedVariant.getContig());
                Assert.assertEquals(actualVariant.getStart(), expectedVariant.getStart());
                Assert.assertEquals(actualVariant.getEnd(), expectedVariant.getEnd());
                Assert.assertEquals(actualVariant.getReference(), expectedVariant.getReference());
                Assert.assertEquals(actualVariant.getAlternateAlleles(), expectedVariant.getAlternateAlleles());
                Assert.assertEquals(actualVariant.getAttributes(), expectedVariant.getAttributes());
            }
            Assert.assertFalse(actualIterator.hasNext());
            Assert.assertFalse(expectedIterator.hasNext());
        }
    }
}
