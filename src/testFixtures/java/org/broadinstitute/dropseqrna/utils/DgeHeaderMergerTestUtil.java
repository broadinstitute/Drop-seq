/*
 * MIT License
 *
 * Copyright 2023 Broad Institute
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

import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeader;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderCodec;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderLibrary;

import java.io.File;
import java.io.IOException;
import java.util.Random;

public class DgeHeaderMergerTestUtil {
    private static final File HG19 = new File("/broad/mccarroll/software/metadata/individual_reference/hg19/hg19.fasta");
    private static final String DGE_FILE_EXTENSION = ".digital_expression.txt.gz";
    private static final Random random = new Random();

    /**
     * @return Temporary file, which will deleteOnExit
     */
    public static File writeHeaderToFile(final DgeHeader header) {
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

    public static DgeHeader makeHeader(final int numLibraries) {
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

    public static String generateUniqueString() {
        return Integer.toString(random.nextInt());
    }
}
