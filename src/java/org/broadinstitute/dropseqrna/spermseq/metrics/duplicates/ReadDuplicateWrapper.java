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
 package org.broadinstitute.dropseqrna.spermseq.metrics.duplicates;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.dropseqrna.utils.CountChangingIteratorWrapper;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;

import java.util.Iterator;

/**
 * Makes a new tag out of the chromosome and position and adds it to the read.
 * This uses the same methodology as Picard's MarkDuplicates, where fragments
 * are tagged from the unclipped 5' end of the genome.
 *
 * @author nemesh
 *
 */
public class ReadDuplicateWrapper extends
		CountChangingIteratorWrapper<SAMRecord> {

	private final String tag;

	protected ReadDuplicateWrapper(final Iterator<SAMRecord> underlyingIterator, final String tag) {
		super(underlyingIterator);
		this.tag = tag;
	}

	@Override
	protected void processRecord(final SAMRecord rec) {
		// String recName = rec.getReadName();
		int coordinate = getCoordinate(rec);
		//Interval i = new Interval(rec.getReferenceName(), coordinate, coordinate);
		//String value = IntervalTagComparator.toString(i);
		// encode the string manually as an interval, don't use the full encoding as it's not needed!
		String value = rec.getReferenceName() + IntervalTagComparator.ENCODE_DELIMITER + coordinate;
		rec.setAttribute(this.tag, value);
		queueRecordForOutput(rec);
	}

	static int getCoordinate (final SAMRecord rec) {
		int coordinate=-1;
		if (rec.getReadNegativeStrandFlag())
			coordinate=rec.getUnclippedEnd();
		else
			coordinate=rec.getUnclippedStart();
		if (coordinate <1) coordinate=1;  // thanks for the negative coordinates, makes sense.
		return coordinate;
	}

}
