/*
 * MIT License
 *
 * Copyright 2026 Broad Institute
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
package org.broadinstitute.dropseqrna.utils.readiterators.umisort;

import htsjdk.samtools.util.SortingCollection;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import static org.broadinstitute.dropseqrna.utils.readiterators.umisort.UMISortCodecIO.*;

/**
 * Spills {@link UMIRead}s to disk for sorting by read name. A flag byte says which of the nullable fields
 * are present, followed by the strings.
 */
public class UMIReadCodec implements SortingCollection.Codec<UMIRead> {
    private static final int NEGATIVE_STRAND = 1;
    private static final int HAS_NAME = 2;
    private static final int HAS_STRAND = 4;
    private static final int HAS_FUNCTION = 8;

    private DataOutputStream out;
    private DataInputStream in;

    @Override
    public void setOutputStream(final OutputStream os) {
        this.out = new DataOutputStream(os);
    }

    @Override
    public void setInputStream(final InputStream is) {
        this.in = new DataInputStream(is);
    }

    @Override
    public void encode(final UMIRead r) {
        try {
            int flags = 0;
            if (r.negativeStrand) flags |= NEGATIVE_STRAND;
            if (r.readName != null) flags |= HAS_NAME;
            if (r.geneStrand != null) flags |= HAS_STRAND;
            if (r.geneFunction != null) flags |= HAS_FUNCTION;
            out.writeByte(flags);
            writeString(out, r.cell);
            writeString(out, r.gene);
            writeString(out, r.umi);
            if (r.readName != null) writeString(out, r.readName);
            if (r.geneStrand != null) writeString(out, r.geneStrand);
            if (r.geneFunction != null) writeString(out, r.geneFunction);
        } catch (IOException e) {
            throw wrap(e);
        }
    }

    @Override
    public UMIRead decode() {
        try {
            final int flags = in.read();
            if (flags == -1) {
                return null;
            }
            final String cell = readString(in);
            final String gene = readString(in);
            final String umi = readString(in);
            final String name = (flags & HAS_NAME) != 0 ? readString(in) : null;
            final String strand = (flags & HAS_STRAND) != 0 ? readString(in) : null;
            final String function = (flags & HAS_FUNCTION) != 0 ? readString(in) : null;
            return new UMIRead(name, (flags & NEGATIVE_STRAND) != 0, cell, gene, strand, function, umi);
        } catch (IOException e) {
            throw wrap(e);
        }
    }

    @Override
    public SortingCollection.Codec<UMIRead> clone() {
        return new UMIReadCodec();
    }
}
