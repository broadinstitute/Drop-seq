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

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableAdapter;

import java.io.BufferedReader;
import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class DgeHeaderMerger {
    public enum Stringency {
        STRICT, // Every library must have a reference, and they must all match; and expression format must match.
        LENIENT,// Reference must match or be undefined; likewise for expression format.
        NONE}   // Anything goes (but UEI collisions are always an error).

    private final Stringency stringency;
    private final Set<DgeHeader.ExpressionFormat> formatSet = new HashSet<>();
    private  File reference = null;
    private final Set<String> ueiSet = new HashSet<>();
    private final Set<String> prefixSet = new HashSet<>();
    private DgeHeader mergedHeader = null;

    public DgeHeaderMerger(Stringency stringency) {
        this.stringency = stringency;
    }

    public void mergeDgeHeader(final DgeHeader header) {
        if (mergedHeader == null) {
            mergedHeader = new DgeHeader();
            mergedHeader.setExpressionFormat(header.getExpressionFormat());
        } else if (mergedHeader.getExpressionFormat() != header.getExpressionFormat()) {
            if (stringency == Stringency.STRICT) {
                throw new DgeMergerException("Expression formats do not agree when merging DGE");
            } else if (mergedHeader.getExpressionFormat() == DgeHeader.ExpressionFormat.unknown) {
                mergedHeader.setExpressionFormat(header.getExpressionFormat());
            } else if (header.getExpressionFormat() != DgeHeader.ExpressionFormat.unknown && stringency != Stringency.NONE) {
                throw new DgeMergerException("Expression formats do not agree when merging DGE");
            }
        }
        if (header.getNumLibraries() == 0) {
            if (stringency == Stringency.STRICT) {
                throw new DgeMergerException("DGE Header does have any libraries");
            }
        } else {
            for (final String prefix: new IterableAdapter<>(header.iterateLibraries())) {
                addLibrary(header.getLibrary(prefix));
            }
        }
    }

    private void addLibrary(final DgeHeaderLibrary lib) {
        if (ueiSet.contains(lib.getUei())) {
            throw new DgeMergerException("UEI " + lib.getUei() + " appears more than once in DGE headers to be merged");
        }
        final String prefix = lib.getPrefix();
        if (prefix == null) {
            if (stringency == Stringency.STRICT) {
                throw new DgeMergerException("Prefix is required when merging DGE headers.");
            }
        } else if (prefixSet.contains(prefix)) {
            throw new DgeMergerException("Prefix collision during DGE header merging: " + prefix);
        }
        if (lib.getReference() == null) {
            if (stringency == Stringency.STRICT) {
                throw new DgeMergerException("Reference cannot be null when doing strict DGE header merging");
            }
        } else if (reference == null) {
            reference = lib.getReference();
        } else if (!reference.equals(lib.getReference()) && stringency != Stringency.NONE) {
            throw new DgeMergerException("References must be null or must agree when doing STRICT or LENIENT DGE header merging");
        }
        ueiSet.add(lib.getUei());
        if (prefix != null) {
            prefixSet.add(prefix);
        }
        mergedHeader.addLibrary(lib);
    }

    public DgeHeader getMergedHeader() {
        return mergedHeader;
    }

    public static DgeHeader mergeDgeHeaders(final List<File> input, final List<String> prefix, final Stringency stringency) {
        final DgeHeaderMerger headerMerger = new DgeHeaderMerger(stringency);
        final DgeHeaderCodec codec = new DgeHeaderCodec();
        for (int i = 0; i < input.size(); ++i) {
            final File file = input.get(i);
            final BufferedReader reader = IOUtil.openFileForBufferedReading(file);
            final DgeHeader dgeHeader = codec.decode(reader, file.getAbsolutePath());
            CloserUtil.close(reader);
            if (!prefix.isEmpty()) {
                if (dgeHeader.getNumLibraries() > 1) {
                    throw new DgeMergerException("Cannot set PREFIX when input DGE has more than one LIBRARY");
                } else if (dgeHeader.getNumLibraries() == 1){
                    final DgeHeaderLibrary lib = dgeHeader.getLibrary(0);
                    if (stringency == DgeHeaderMerger.Stringency.STRICT && lib.getPrefix() != null) {
                        throw new DgeMergerException("Overwriting prefix in DGE header for " + file.getAbsolutePath());
                    }
                    lib.setPrefix(prefix.get(i));
                }
            }
            headerMerger.mergeDgeHeader(dgeHeader);
        }
        return headerMerger.getMergedHeader();

    }

    public static class DgeMergerException extends RuntimeException {
        public DgeMergerException() {
        }

        public DgeMergerException(String message) {
            super(message);
        }

        public DgeMergerException(String message, Throwable cause) {
            super(message, cause);
        }

        public DgeMergerException(Throwable cause) {
            super(cause);
        }

        public DgeMergerException(String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace) {
            super(message, cause, enableSuppression, writableStackTrace);
        }
    }
}
