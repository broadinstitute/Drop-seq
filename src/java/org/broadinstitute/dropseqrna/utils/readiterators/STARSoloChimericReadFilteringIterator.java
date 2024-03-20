/*
 * MIT License
 *
 * Copyright 2024 Broad Institute
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
package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import java.util.Iterator;

/**
 * Filter out reads that STARsolo has marked as chimeric.
 */
public class STARSoloChimericReadFilteringIterator
        extends FilteredIterator<SAMRecord> {

    private final String molecularBarcodeTag;

    public STARSoloChimericReadFilteringIterator(Iterator<SAMRecord> underlyingIterator, String molecularBarcodeTag) {
        super(underlyingIterator);
        this.molecularBarcodeTag = molecularBarcodeTag;
    }

    @Override
    public boolean filterOut(SAMRecord rec) {
        final String umi = rec.getStringAttribute(molecularBarcodeTag);
        if (umi == null) {
            throw new IllegalArgumentException("Read " + rec.getSAMString() + " does not have a UMI tag " + molecularBarcodeTag);
        }
        return umi.equals(CHIMERIC_UMI);
    }

    /**
     * The UMI value that STARsolo uses to indicate a chimeric read.
     */
    public static final String CHIMERIC_UMI = "-";
}
