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

import htsjdk.samtools.SAMRecord;

public class EditDistanceFilteringIterator extends FilteredIterator<SAMRecord> {

	private Integer maxEditDistance;
	private static final String EDIT_DISTANCE_TAG="NM";

	public EditDistanceFilteringIterator (final Iterator<SAMRecord> underlyingIterator, final Integer maxEditDistance) {
		super(underlyingIterator);
		this.maxEditDistance=maxEditDistance;
	}

	@Override
	public boolean filterOut(final SAMRecord r) {
		if (maxEditDistance==null) return false;
		Object o = r.getAttribute(EDIT_DISTANCE_TAG);
		if (o==null) return false;
		int readEditDistance = ((Integer) o).intValue();
    	return readEditDistance > this.maxEditDistance;
    }
}
