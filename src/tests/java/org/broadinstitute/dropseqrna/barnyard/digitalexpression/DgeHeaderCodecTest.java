/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2016 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
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
        header.setExpressionFormat(null);

        final StringWriter writer = new StringWriter();
        final String contentAfterHeader = "Goodbye, cruel world!";
        writer.append(contentAfterHeader);
        final BufferedReader reader = new BufferedReader(new StringReader(writer.toString()));
        final DgeHeader decodedHeader = new DgeHeaderCodec().decode(reader, "test input");
        Assert.assertEquals(decodedHeader, header);
        final String actualContentAfterHeader = reader.readLine();
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
