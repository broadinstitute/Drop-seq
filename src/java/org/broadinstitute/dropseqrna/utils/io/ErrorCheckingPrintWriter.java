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
