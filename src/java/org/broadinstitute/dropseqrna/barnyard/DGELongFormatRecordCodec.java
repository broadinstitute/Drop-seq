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
package org.broadinstitute.dropseqrna.barnyard;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.SortingCollection.Codec;
import picard.PicardException;

public class DGELongFormatRecordCodec implements SortingCollection.Codec<DGELongFormatRecord> {

	private DataOutputStream outputStream = null;
	private DataInputStream inputReader = null;

	@Override
	public void setOutputStream(final OutputStream stream) {
		// this.outputStream = new ErrorCheckingPrintStream(stream);
		this.outputStream = new DataOutputStream(stream);
	}

	@Override
	public void setInputStream(final InputStream stream) {
		this.inputReader = new DataInputStream(stream);
	}

	@Override
	public void encode(final DGELongFormatRecord val) {
		try {
			this.outputStream.writeUTF(val.getCell());
			this.outputStream.writeUTF(val.getGene());
			this.outputStream.writeInt(val.getCount());
		} catch (final IOException ioe) {
			throw new RuntimeIOException("Could not encode DGELongFormat record for a sorting collection: " + ioe.getMessage(), ioe);
		}
	}

	@Override
	public DGELongFormatRecord decode() {
		try {
            return new DGELongFormatRecord (this.inputReader.readUTF(), this.inputReader.readUTF(), this.inputReader.readInt());
        } catch (EOFException e) {
            return null;
        } catch (IOException e) {
            throw new PicardException("Exception reading DGELongFormat from temporary file.", e);
        }
	}



	@Override
	public Codec<DGELongFormatRecord> clone() {
		return new DGELongFormatRecordCodec();
	}

}