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

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashSet;

public class ReportFileUtil {

    public static final String COMMENT_PREFIX = "#";

    public static Collection<String> readComments(File reportFile) {
        Collection<String> comments = new ArrayList<>();

        BufferedReader reader;
        try {
            IOUtil.assertFileIsReadable(reportFile);
            reader = IOUtil.openFileForBufferedReading(reportFile);

            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(COMMENT_PREFIX))
                    comments.add(line);
                else if (line.trim().length() > 0)
                    break;
            }
        } catch (IOException ex) {
            throw new RuntimeIOException("Exception reading " + reportFile.getAbsolutePath());
        }
        CloserUtil.close(reader);

        return comments;
    }

    public static Collection<String> readComments(Collection<File> reportFiles) {
        Collection<String> comments = new LinkedHashSet<>();
        for (File reportFile : reportFiles) {
            comments.addAll(readComments(reportFile));
        }

        return comments;
    }
}
