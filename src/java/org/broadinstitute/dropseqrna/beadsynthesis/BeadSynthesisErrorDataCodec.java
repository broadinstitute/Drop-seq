/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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

package org.broadinstitute.dropseqrna.beadsynthesis;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;

import org.broadinstitute.dropseqrna.TranscriptomeException;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.SortingCollection.Codec;
import picard.PicardException;

public class BeadSynthesisErrorDataCodec implements SortingCollection.Codec<BeadSynthesisErrorData> {

	private ObjectOutputStream outputStream = null;
	private ObjectInputStream inputReader = null;

	@Override
	public Codec<BeadSynthesisErrorData> clone() {
		return new BeadSynthesisErrorDataCodec();
	}

	@Override
	public BeadSynthesisErrorData decode() {
		try {
			Object o = this.inputReader.readObject();
			return (BeadSynthesisErrorData) o;
        } catch (EOFException e) {
            return null;
        } catch (IOException e) {
            throw new PicardException("Exception reading BeadSynthesisErrorData from temporary file.", e);
        } catch (ClassNotFoundException cnfe) {
        	throw new PicardException("Could not find class BeadSynthesisErrorData");
        }

	}

	@Override
	public void encode(final BeadSynthesisErrorData val) {
		try {
			this.outputStream.writeObject(val);
		} catch (final IOException ioe) {

			throw new RuntimeIOException("Could not encode BeadSynthesisErrorData record for a sorting collection: " + ioe.getMessage(), ioe);
		}

	}

	@Override
	public void setInputStream(final InputStream stream) {
		try {
			this.inputReader = new ObjectInputStream(stream);
		} catch (IOException e) {
			throw new TranscriptomeException(e.getMessage());
		}
	}

	@Override
	public void setOutputStream(final OutputStream stream) {
		try {
			this.outputStream = new ObjectOutputStream(stream);
		} catch (IOException e) {
			throw new TranscriptomeException(e.getMessage());
		}

	}


}
