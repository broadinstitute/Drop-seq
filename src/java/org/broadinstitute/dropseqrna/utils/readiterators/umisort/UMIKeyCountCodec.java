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

/** Spills {@link UMIKeyCount}s to disk: cell, gene, UMI strings, and a varint count. */
public class UMIKeyCountCodec implements SortingCollection.Codec<UMIKeyCount> {
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
    public void encode(final UMIKeyCount k) {
        try {
            writeString(out, k.cell);
            writeString(out, k.gene);
            writeString(out, k.umi);
            writeVarInt(out, k.getCount());
        } catch (IOException e) {
            throw wrap(e);
        }
    }

    @Override
    public UMIKeyCount decode() {
        try {
            // The first varint of a record is the cell length, so -1 here means a clean EOF.
            final int cellLength = readVarInt(in);
            if (cellLength < 0) {
                return null;
            }
            final byte[] cellBytes = new byte[cellLength];
            in.readFully(cellBytes);
            final String cell = new String(cellBytes, java.nio.charset.StandardCharsets.UTF_8);
            final String gene = readString(in);
            final String umi = readString(in);
            final int count = readVarInt(in);
            return new UMIKeyCount(cell, gene, umi, count);
        } catch (IOException e) {
            throw wrap(e);
        }
    }

    @Override
    public SortingCollection.Codec<UMIKeyCount> clone() {
        return new UMIKeyCountCodec();
    }
}
