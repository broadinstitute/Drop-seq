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
package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Iterator;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import htsjdk.samtools.SAMRecord;

public class BamTagCountingIterator extends FilteredIterator<SAMRecord> {

	private String tag;
    private ObjectCounter<String> counter;

	public BamTagCountingIterator(final Iterator<SAMRecord> underlyingIterator, final String tag) {
		super(underlyingIterator);
		if (tag!=null) this.counter = new ObjectCounter<>();
		this.tag=tag;
	}

	@Override
	public boolean filterOut(final SAMRecord rec) {
		if (this.tag==null) return false;
		String value = rec.getStringAttribute(tag);
		if (value!=null) counter.increment(value);
		return false;
	}

	public ObjectCounter<String> getCounts () {
		return this.counter;
	}


}
