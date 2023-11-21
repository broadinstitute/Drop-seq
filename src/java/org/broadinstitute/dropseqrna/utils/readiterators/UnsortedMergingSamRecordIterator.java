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
package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.dropseqrna.utils.IteratorOfIterators;

import java.util.Iterator;
import java.util.List;

/**
 * Iterate over the SAMRecords in order of readers, and for each reader, in the order they appear in the filw.
 */
public class UnsortedMergingSamRecordIterator
implements CloseableIterator<SAMRecord> {
    private final SAMFileHeader outputHeader;
    private final List<SamReader> readers;
    private final IteratorOfIterators<SAMRecord> it;

    /**
     *
     * @param outputHeader If non-null, each SAMRecord has its header set to this.
     */
    public UnsortedMergingSamRecordIterator(SAMFileHeader outputHeader, List<SamReader> readers) {
        this.outputHeader = outputHeader;
        this.readers = readers;
        this.it = new IteratorOfIterators<>(this.readers.stream().map(reader -> (Iterator<SAMRecord>)(reader.iterator())).iterator());
    }

    @Override
    public void close() {
        CloserUtil.close(readers);
    }

    @Override
    public boolean hasNext() {
        return it.hasNext();
    }

    @Override
    public SAMRecord next() {
        final SAMRecord next = it.next();
        if (outputHeader != null) {
            next.setHeader(outputHeader);
        }
        return next;
    }
}
