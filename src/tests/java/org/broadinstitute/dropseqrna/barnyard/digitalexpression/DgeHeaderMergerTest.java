/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
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

import htsjdk.samtools.util.RuntimeIOException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class DgeHeaderMergerTest {

    private static final File HG19 = new File("/broad/mccarroll/software/metadata/individual_reference/hg19/hg19.fasta");
    private static final File HG19_MM10_TRANSGENES = new File("/broad/mccarroll/software/metadata/merged_reference/hg19_mm10_transgenes/hg19_mm10_transgenes.fasta");
    private static final String DGE_FILE_EXTENSION = ".digital_expression.txt.gz";

    private Random random = new Random();

    @Test(dataProvider = "stringencyDataProvider")
    public void testExpressionFormatMismatch(final DgeHeaderMerger.Stringency stringency) {
        final DgeHeader h1 = makeHeader(1);
        final DgeHeader h2 = makeHeader(1);
        h1.setExpressionFormat(DgeHeader.ExpressionFormat.raw);
        h2.setExpressionFormat(DgeHeader.ExpressionFormat.log10_normalized);
        final DgeHeader mergedHeader = mergeHeaders(h1, h2, stringency);
        switch (stringency) {
            case NONE:
                Assert.assertNotNull(mergedHeader); break;
            case LENIENT:
            case STRICT:
                Assert.assertNull(mergedHeader);  break;
        }
    }

    @Test(dataProvider = "stringencyDataProvider")
    public void testExpressionFormatFirstUnknown(final DgeHeaderMerger.Stringency stringency) {
        final DgeHeader h1 = makeHeader(1);
        final DgeHeader h2 = makeHeader(1);
        h1.setExpressionFormat(DgeHeader.ExpressionFormat.unknown);
        h2.setExpressionFormat(DgeHeader.ExpressionFormat.log10_normalized);
        final DgeHeader mergedHeader = mergeHeaders(h1, h2, stringency);
        switch (stringency) {
            case NONE:
            case LENIENT:
                Assert.assertNotNull(mergedHeader); break;
            case STRICT:
                Assert.assertNull(mergedHeader);  break;
        }
    }

    @Test(dataProvider = "stringencyDataProvider")
    public void testExpressionFormatSecondUnknown(final DgeHeaderMerger.Stringency stringency) {
        final DgeHeader h1 = makeHeader(1);
        final DgeHeader h2 = makeHeader(1);
        h1.setExpressionFormat(DgeHeader.ExpressionFormat.raw);
        h2.setExpressionFormat(DgeHeader.ExpressionFormat.unknown);
        final DgeHeader mergedHeader = mergeHeaders(h1, h2, stringency);
        switch (stringency) {
            case NONE:
            case LENIENT:
                Assert.assertNotNull(mergedHeader); break;
            case STRICT:
                Assert.assertNull(mergedHeader);  break;
        }
    }

    @Test(dataProvider = "stringencyDataProvider")
    public void testNoLibraries(final DgeHeaderMerger.Stringency stringency) {
        final DgeHeader h1 = makeHeader(1);
        final DgeHeader h2 = makeHeader(0);
        final DgeHeader mergedHeader = mergeHeaders(h1, h2, stringency);
        switch (stringency) {
            case NONE:
            case LENIENT:
                Assert.assertNotNull(mergedHeader); break;
            case STRICT:
                Assert.assertNull(mergedHeader);  break;
        }
    }

    @Test(dataProvider = "stringencyDataProvider")
    public void testNoPrefix(final DgeHeaderMerger.Stringency stringency) {
        final DgeHeader h1 = makeHeader(1);
        final DgeHeader h2 = makeHeader(1);
        h2.getLibrary(0).setPrefix(null);

        final DgeHeader mergedHeader = mergeHeaders(h1, h2, stringency);
        switch (stringency) {
            case NONE:
            case LENIENT:
                Assert.assertNotNull(mergedHeader); break;
            case STRICT:
                Assert.assertNull(mergedHeader);  break;
        }
    }

    @Test(dataProvider = "stringencyDataProvider")
    public void testReferenceMismatch(final DgeHeaderMerger.Stringency stringency) {
        final DgeHeader h1 = makeHeader(1);
        final DgeHeader h2 = makeHeader(1);
        h2.getLibrary(0).setReference(HG19_MM10_TRANSGENES);

        final DgeHeader mergedHeader = mergeHeaders(h1, h2, stringency);
        switch (stringency) {
            case NONE:
                Assert.assertNotNull(mergedHeader); break;
            case LENIENT:
            case STRICT:
                Assert.assertNull(mergedHeader);  break;
        }
    }

    @Test(dataProvider = "stringencyDataProvider")
    public void testOneReferenceNull(final DgeHeaderMerger.Stringency stringency) {
        final DgeHeader h1 = makeHeader(1);
        final DgeHeader h2 = makeHeader(1);
        h2.getLibrary(0).setReference(null);

        final DgeHeader mergedHeader = mergeHeaders(h1, h2, stringency);
        switch (stringency) {
            case NONE:
            case LENIENT:
                Assert.assertNotNull(mergedHeader); break;
            case STRICT:
                Assert.assertNull(mergedHeader);  break;
        }
    }

    @Test(dataProvider = "stringencyDataProvider")
    public void testBothReferenceNull(final DgeHeaderMerger.Stringency stringency) {
        final DgeHeader h1 = makeHeader(1);
        h1.getLibrary(0).setReference(null);
        final DgeHeader h2 = makeHeader(1);
        h2.getLibrary(0).setReference(null);

        final DgeHeader mergedHeader = mergeHeaders(h1, h2, stringency);
        switch (stringency) {
            case NONE:
            case LENIENT:
                Assert.assertNotNull(mergedHeader); break;
            case STRICT:
                Assert.assertNull(mergedHeader);  break;
        }
    }

    @Test(dataProvider = "stringencyDataProvider")
    public void testUeiCollision(final DgeHeaderMerger.Stringency stringency) {
        final DgeHeader h1 = makeHeader(1);
        final DgeHeaderLibrary lib1 = h1.getLibrary(0);
        final DgeHeader h2 = makeHeader(0);
        final DgeHeaderLibrary lib2 = new DgeHeaderLibrary(lib1.getUei());
        h2.addLibrary(lib2);

        final DgeHeader mergedHeader = mergeHeaders(h1, h2, stringency);
        switch (stringency) {
            case NONE:
                Assert.assertNotNull(mergedHeader); break;
            case LENIENT:
                Assert.assertNotNull(mergedHeader); break;
            case STRICT:
                Assert.assertNull(mergedHeader);  break;
        }
    }

    @Test(dataProvider = "stringencyDataProvider")
    public void testPrefixCollision(final DgeHeaderMerger.Stringency stringency) {
        final DgeHeader h1 = makeHeader(1);
        final DgeHeaderLibrary lib1 = h1.getLibrary(0);
        final DgeHeader h2 = makeHeader(1);
        final DgeHeaderLibrary lib2 = h2.removeLibrary(0);
        lib2.setPrefix(lib1.getPrefix());
        h2.addLibrary(lib2);

        final DgeHeader mergedHeader = mergeHeaders(h1, h2, stringency);
        switch (stringency) {
            case NONE:
            case LENIENT:
            case STRICT:
                Assert.assertNull(mergedHeader);  break;
        }
    }

    @Test
    public void testMergePrefixFromInputs() {
        final DgeHeader h1 = makeHeader(1);
        // It's fine to merge an input DGE with multiple libraries, as long as the input DGE has the prefixes.
        final DgeHeader h2 = makeHeader(2);
        Assert.assertNotNull(mergeHeaders(h1, h2, Collections.<String>emptyList(), DgeHeaderMerger.Stringency.STRICT));
    }

    @Test
    public void testMergePrefixFromCommandLine() {
        final DgeHeader h1 = makeHeader(1);
        final DgeHeader h2 = makeHeader(1);
        final DgeHeaderLibrary l1 = h1.getLibrary(0);
        final DgeHeaderLibrary l2 = h2.getLibrary(0);
        final List<String> prefix = Arrays.asList(l1.getPrefix(), l2.getPrefix());
        l1.setPrefix(null);
        l2.setPrefix(null);
        Assert.assertNotNull(mergeHeaders(h1, h2, prefix, DgeHeaderMerger.Stringency.STRICT));
    }

    @Test(dataProvider = "stringencyDataProvider")
    public void testMergeOverwritePrefixFromCommandLine(final DgeHeaderMerger.Stringency stringency) {
        final DgeHeader h1 = makeHeader(1);
        final DgeHeader h2 = makeHeader(1);
        final List<String> prefix = Arrays.asList(generateUniqueString(), generateUniqueString());
        h2.getLibrary(0).setPrefix(null);
        final DgeHeader mergedHeader = mergeHeaders(h1, h2, prefix, stringency);
        switch (stringency) {
            case NONE:
            case LENIENT:
                Assert.assertNotNull(mergedHeader); break;
            case STRICT:
                Assert.assertNull(mergedHeader);  break;
        }
    }

    @Test(dataProvider = "stringencyDataProvider")
    public void testMergeMultiplePrefixFromCommandLine(final DgeHeaderMerger.Stringency stringency) {
        final DgeHeader h1 = makeHeader(1);
        final DgeHeader h2 = makeHeader(2);
        final List<String> prefix = Arrays.asList(generateUniqueString(), generateUniqueString());
        final DgeHeader mergedHeader = mergeHeaders(h1, h2, prefix, stringency);
        switch (stringency) {
            case NONE:
            case LENIENT:
            case STRICT:
                Assert.assertNull(mergedHeader);  break;
        }
    }

    /**
     * @return the merged header, or null if there was an exception
     */
    private DgeHeader mergeHeaders(final DgeHeader h1,
                                   final DgeHeader h2,
                                   final DgeHeaderMerger.Stringency stringency) {
        try {
            final DgeHeaderMerger merger = new DgeHeaderMerger(stringency);
            merger.mergeDgeHeader(h1);
            merger.mergeDgeHeader(h2);
            return merger.getMergedHeader();
        } catch (DgeHeaderMerger.DgeMergerException e) {
            System.err.println(getClass().getSimpleName() + "." + getMethodName(1) + ":" + stringency + ": " + e.getMessage());
            return null;
        }
    }

    /**
     * Command-line merge functionality.
     * @return the merged header, or null if there was an exception
     */
    private DgeHeader mergeHeaders(final DgeHeader h1,
                                   final DgeHeader h2,
                                   final List<String> prefix,
                                   final DgeHeaderMerger.Stringency stringency) {
        final File f1 = writeHeaderToFile(h1);
        final File f2 = writeHeaderToFile(h2);
        final List<File> input = Arrays.asList(f1, f2);
        try {
            return DgeHeaderMerger.mergeDgeHeaders(input, prefix, stringency);
        } catch (DgeHeaderMerger.DgeMergerException e) {
            System.err.println(getClass().getSimpleName() + "." + getMethodName(1) + ":" + stringency + ": " + e.getMessage());
            return null;
        }
    }

    /**
     * @return Temporary file, which will deleteOnExit
     */
    public File writeHeaderToFile(final DgeHeader header) {
        try {
            final File f = File.createTempFile("DgeHeaderMergerTest.", DGE_FILE_EXTENSION);
            f.deleteOnExit();
            final DgeHeaderCodec codec = new DgeHeaderCodec();
            codec.encode(f, header);
            return f;
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    @DataProvider(name="stringencyDataProvider")
    public Object[][] stringencyDataProvider() {
        return new Object[][] {
                {DgeHeaderMerger.Stringency.STRICT},
                {DgeHeaderMerger.Stringency.LENIENT},
                {DgeHeaderMerger.Stringency.NONE},
        };
    }

    public DgeHeader makeHeader(final int numLibraries) {
        final DgeHeader ret = new DgeHeader();
        ret.setExpressionFormat(DgeHeader.ExpressionFormat.raw);
        for (int i = 0; i < numLibraries; ++i) {
            final DgeHeaderLibrary lib = new DgeHeaderLibrary(generateUniqueString());
            lib.setPrefix(generateUniqueString());
            lib.setReference(HG19);
            ret.addLibrary(lib);
        }
        return ret;
    }

    private String generateUniqueString() {
        return Integer.toString(random.nextInt());
    }

    private static String getMethodName(final int depth) {
        final StackTraceElement[] stack = Thread.currentThread().getStackTrace();
        // [0] == getStackTrace(), [1] == getMethodName, so nearest caller is [2]
        Assert.assertTrue(stack.length >= depth + 3);
        return stack[depth + 2].getMethodName();
    }
}
