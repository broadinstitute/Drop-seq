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
package org.broadinstitute.dropseqrna.utils.io;

import java.io.*;

public class ErrorCheckingPrintWriter extends PrintWriter {

    private String path = "unknown";

    public ErrorCheckingPrintWriter(Writer out) {
        super(out);
    }

    public ErrorCheckingPrintWriter(Writer out, boolean autoFlush) {
        super(out, autoFlush);
    }

    public ErrorCheckingPrintWriter(OutputStream out) {
        super(out);
    }

    public ErrorCheckingPrintWriter(OutputStream out, boolean autoFlush) {
        super(out, autoFlush);
    }

    public ErrorCheckingPrintWriter(String fileName) throws FileNotFoundException {
        super(fileName);
        path = fileName;
    }

    public ErrorCheckingPrintWriter(String fileName, String csn) throws FileNotFoundException, UnsupportedEncodingException {
        super(fileName, csn);
        path = fileName;
    }

    public ErrorCheckingPrintWriter(File file) throws FileNotFoundException {
        super(file);
        path = file.getAbsolutePath();
    }

    public ErrorCheckingPrintWriter(File file, String csn) throws FileNotFoundException, UnsupportedEncodingException {
        super(file, csn);
        path = file.getAbsolutePath();
    }

    private final ThreadLocal<Boolean> nestedFlush = new ThreadLocal<Boolean>() {
        @Override
        protected Boolean initialValue() {
            return false;
        }
    };

    @Override
    public void flush() {
        // checkError calls flush() so if inside flush, don't recurse.
        if (nestedFlush.get()) return;
        try {
            nestedFlush.set(true);
            if (checkError()) {
                throw new RuntimeException("Exception writing file " + path);
            }
            super.flush();
        } finally {
            nestedFlush.set(false);
        }
    }

    @Override
    public void close() {
        if (checkError()) {
            throw new RuntimeException("Exception writing file " + path);
        }
        super.close();
    }

    /** For better error reporting */
    public void setPath(String path) {
        this.path = path;
    }
}
