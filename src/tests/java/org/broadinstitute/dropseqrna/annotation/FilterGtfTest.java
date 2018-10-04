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

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class FilterGtfTest {
    private static final File INPUT_GTF = new File("testdata/org/broadinstitute/transcriptome/annotation/FilterGtfInput.gtf");
    private static final File SEQUENCE_DICTIONARY = new File("testdata/org/broadinstitute/transcriptome/annotation/FilterGtfInput.dict");

    @Test(dataProvider = "testGeneBiotypeFilterDataProvider")
    public void testGeneBiotypeFilter(final List<String> geneBiotypesToFilter, final int expectedOutputLines) throws IOException {
        final FilterGtf filterGtf = new FilterGtf();
        filterGtf.GTF = INPUT_GTF;
        filterGtf.OUTPUT = File.createTempFile("testGeneBiotypeFilter.", ".gtf");;
        filterGtf.OUTPUT.deleteOnExit();
        filterGtf.UNDESIRED_GENE_BIOTYPE = geneBiotypesToFilter;
        runClp(expectedOutputLines, filterGtf);
    }

    private void runClp(int expectedOutputLines, FilterGtf filterGtf) throws IOException {
        Assert.assertEquals(filterGtf.doWork(), 0);
        final BufferedReader reader = IOUtil.openFileForBufferedReading(filterGtf.OUTPUT);
        int lines = 0;
        while (reader.readLine() != null) {
            ++lines;
        }
        CloserUtil.close(reader);
        Assert.assertEquals(lines, expectedOutputLines);
    }

    @DataProvider(name="testGeneBiotypeFilterDataProvider")
    public Object[][] testGeneBiotypeFilterDataProvider() {
        return new Object[][]{
                {Arrays.asList("test-biotype", "test-biotype2"), 2},
                {Arrays.asList("test-biotype"), 3},
                {Arrays.asList("test-biotypex"), 4},
                {Collections.emptyList(), 4},
        };
    }

    @Test
    public void testSequenceDictionaryFilter() throws IOException {
        final FilterGtf filterGtf = new FilterGtf();
        filterGtf.GTF = INPUT_GTF;
        filterGtf.SEQUENCE_DICTIONARY = SEQUENCE_DICTIONARY;
        filterGtf.OUTPUT = File.createTempFile("testGeneBiotypeFilter.", ".gtf");;
        filterGtf.OUTPUT.deleteOnExit();
        filterGtf.UNDESIRED_GENE_BIOTYPE = Arrays.asList("test-biotype");
        // 1 line filtered by gene_biotype, another filtered by sequence name
        runClp(2, filterGtf);
    }
}
