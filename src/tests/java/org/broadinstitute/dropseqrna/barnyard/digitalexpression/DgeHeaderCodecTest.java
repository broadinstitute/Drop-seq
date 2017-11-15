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

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;

public class DgeHeaderCodecTest {

    @DataProvider(name = "numberOfLibrariesDataProvider")
    public Object[][] numberOfLibrariesDataProvider() {
        return new Object[][]{
                {Integer.valueOf(0)},
                {Integer.valueOf(1)},
                {Integer.valueOf(2)},
        };
    }


    @Test(dataProvider = "numberOfLibrariesDataProvider")
    public void testBasic(final Integer numberOfLibraries) {
        final DgeHeader header = makeDgeHeader(numberOfLibraries);

        final StringWriter writer = new StringWriter();
        final DgeHeaderCodec codec = new DgeHeaderCodec();
        codec.encode(writer, header);
        System.out.print(writer.toString());

        final BufferedReader reader = new BufferedReader(new StringReader(writer.toString()));
        final DgeHeader decodedHeader = codec.decode(reader, "test input");
        Assert.assertEquals(decodedHeader, header);
    }

    /**
     * Confirm that position of reader is immediately after header, after header is decoded.
     */
    @Test(dataProvider = "numberOfLibrariesDataProvider")
    public void testHeaderWithTrailingContents(final Integer numberOfLibraries) throws IOException {
        final DgeHeader header = makeDgeHeader(numberOfLibraries);

        final StringWriter writer = new StringWriter();
        final DgeHeaderCodec codec = new DgeHeaderCodec();
        codec.encode(writer, header);
        final String contentAfterHeader = "Goodbye, cruel world!";
        writer.append(contentAfterHeader);
        System.out.print(writer.toString());
        final BufferedReader reader = new BufferedReader(new StringReader(writer.toString()));
        final DgeHeader decodedHeader = codec.decode(reader, "test input");
        Assert.assertEquals(decodedHeader, header);
        final String actualContentAfterHeader = reader.readLine();
        Assert.assertEquals(actualContentAfterHeader, contentAfterHeader);
    }

    /**
     * Confirm that position of reader is correct after looking for a header and finding none.
     */
    @Test
    public void testNoHeaderWithTrailingContents() throws IOException {
        final DgeHeader header = new DgeHeader();
        // Clear default values so the header looks like a codec not finding a header
        header.setVersion(null);
        header.setExpressionFormat(DgeHeader.ExpressionFormat.unknown);

        final StringWriter writer = new StringWriter();
        final String contentAfterHeader = "Goodbye, cruel world!";
        writer.append(contentAfterHeader);
        final BufferedReader reader = new BufferedReader(new StringReader(writer.toString()));
        final DgeHeader decodedHeader = new DgeHeaderCodec().decode(reader, "test input");
        Assert.assertEquals(decodedHeader, header);
        final String actualContentAfterHeader = reader.readLine();
        Assert.assertEquals(actualContentAfterHeader, contentAfterHeader);
    }


    @Test(dataProvider = "numberOfLibrariesDataProvider")
    public void testBasicInputStream(final Integer numberOfLibraries) {
        final DgeHeader header = makeDgeHeader(numberOfLibraries);

        final StringWriter writer = new StringWriter();
        final DgeHeaderCodec codec = new DgeHeaderCodec();
        codec.encode(writer, header);
        System.out.print(writer.toString());

        final BufferedInputStream inputStream = new BufferedInputStream(new ByteArrayInputStream(writer.toString().getBytes()));
        final DgeHeader decodedHeader = codec.decode(inputStream, "test input");
        Assert.assertEquals(decodedHeader, header);
    }

    /**
     * Confirm that position of reader is immediately after header, after header is decoded.
     */
    @Test(dataProvider = "numberOfLibrariesDataProvider")
    public void testHeaderWithTrailingContentsInputStream(final Integer numberOfLibraries) throws IOException {
        final DgeHeader header = makeDgeHeader(numberOfLibraries);

        final StringWriter writer = new StringWriter();
        final DgeHeaderCodec codec = new DgeHeaderCodec();
        codec.encode(writer, header);
        final String contentAfterHeader = "Goodbye, cruel world!";
        writer.append(contentAfterHeader);
        System.out.print(writer.toString());
        final BufferedInputStream inputStream = new BufferedInputStream(new ByteArrayInputStream(writer.toString().getBytes()));
        final DgeHeader decodedHeader = codec.decode(inputStream, "test input");
        Assert.assertEquals(decodedHeader, header);
        final String actualContentAfterHeader = new BufferedReader(new InputStreamReader(inputStream)).readLine();
        Assert.assertEquals(actualContentAfterHeader, contentAfterHeader);
    }

    /**
     * Confirm that position of reader is correct after looking for a header and finding none.
     */
    @Test
    public void testNoHeaderWithTrailingContentsInputStream() throws IOException {
        final DgeHeader header = new DgeHeader();
        // Clear default values so the header looks like a codec not finding a header
        header.setVersion(null);
        header.setExpressionFormat(DgeHeader.ExpressionFormat.unknown);

        final StringWriter writer = new StringWriter();
        final String contentAfterHeader = "Goodbye, cruel world!";
        writer.append(contentAfterHeader);
        final BufferedInputStream inputStream = new BufferedInputStream(new ByteArrayInputStream(writer.toString().getBytes()));
        final DgeHeader decodedHeader = new DgeHeaderCodec().decode(inputStream, "test input");
        Assert.assertEquals(decodedHeader, header);
        final String actualContentAfterHeader = new BufferedReader(new InputStreamReader(inputStream)).readLine();
        Assert.assertEquals(actualContentAfterHeader, contentAfterHeader);
    }


    private DgeHeader makeDgeHeader(final int numLibraries) {
        final DgeHeader header = new DgeHeader();
        final String testVersion = "1234";
        header.setVersion(testVersion);
        header.setExpressionFormat(DgeHeader.ExpressionFormat.log10_normalized);
        if (numLibraries > 0) {
            final DgeHeaderLibrary lib1 = new DgeHeaderLibrary("uei1");
            lib1.setInput(new File("/this/is/a/dummy/file"));
            lib1.setReference(new File("/dummy/reference"));
            lib1.setTag("test_key", "test_value");
            lib1.setTag("test_key2", "test_value2");
            header.addLibrary(lib1);
            if (numLibraries > 1) {
                final DgeHeaderLibrary lib2 = new DgeHeaderLibrary("uei2");
                lib2.setInputDge(new File("/dummy/dge"));
                lib2.setPrefix("a prefix");
                header.addLibrary(lib2);
            }
        }
        return header;
    }
}
