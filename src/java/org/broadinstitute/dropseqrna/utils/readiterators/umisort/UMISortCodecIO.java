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

import htsjdk.samtools.util.RuntimeIOException;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.charset.StandardCharsets;

/**
 * Length-prefixed UTF-8 string and varint helpers shared by the UMI sorting codecs.
 * Unlike DataOutputStream.writeUTF there is no 64K limit.
 */
final class UMISortCodecIO {
    private UMISortCodecIO() {}

    static void writeVarInt(final DataOutputStream out, int value) throws IOException {
        while ((value & ~0x7F) != 0) {
            out.writeByte((value & 0x7F) | 0x80);
            value >>>= 7;
        }
        out.writeByte(value);
    }

    /** @return the value, or -1 if the stream is at EOF before the first byte. Values must be non-negative. */
    static int readVarInt(final DataInputStream in) throws IOException {
        int b = in.read();
        if (b == -1) {
            return -1;
        }
        int value = 0;
        int shift = 0;
        while ((b & 0x80) != 0) {
            value |= (b & 0x7F) << shift;
            shift += 7;
            b = in.read();
            if (b == -1) {
                throw new EOFException("Truncated varint");
            }
        }
        return value | (b << shift);
    }

    static void writeString(final DataOutputStream out, final String s) throws IOException {
        final byte[] bytes = s.getBytes(StandardCharsets.UTF_8);
        writeVarInt(out, bytes.length);
        out.write(bytes);
    }

    static String readString(final DataInputStream in) throws IOException {
        final int length = readVarInt(in);
        if (length < 0) {
            throw new EOFException("Truncated string");
        }
        final byte[] bytes = new byte[length];
        in.readFully(bytes);
        return new String(bytes, StandardCharsets.UTF_8);
    }

    static RuntimeIOException wrap(final IOException e) {
        return new RuntimeIOException(e);
    }
}
