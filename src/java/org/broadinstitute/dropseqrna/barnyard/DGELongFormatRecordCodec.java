/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
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