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
import htsjdk.samtools.SAMTagUtil;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;

/**
 * Takes in reads that have been sorted in cell/umi/position order and groups them into clusters that are <threshold> size from start to end.
 * Each cluster is for a single cell/umi.
 * @author nemesh
 *
 */
public class ReadClusterIterator implements CloseableIterator<ReadCluster> {

	private final GroupingIterator<SAMRecord> iterator;
	private final int windowSize;
	private PeekableIterator<SAMRecord> reads; // hold the group of reads here.  consume this before getting the next group.

	private final short cellBarcodeTag;
	private final short umiBarcodeTag;

	public ReadClusterIterator (final GroupingIterator<SAMRecord> iterator, final int windowSize, final String cellBarcodeTag, final String umiBarcodeTag) {
		this.iterator=iterator;
		this.windowSize=windowSize;
		this.cellBarcodeTag = SAMTagUtil.getSingleton().makeBinaryTag(cellBarcodeTag);
		this.umiBarcodeTag = SAMTagUtil.getSingleton().makeBinaryTag(umiBarcodeTag);
	}

	@Override
	/**
	 * Get the next cluster of reads that occur within the set threshold distance.
	 */
	public ReadCluster next() {
		// populate reads from the next group
		if (reads==null || !reads.hasNext())
			reads=new PeekableIterator<SAMRecord>(iterator.next().iterator());

		String cellBC = (String) reads.peek().getAttribute(this.cellBarcodeTag);
		String umiBC = (String)reads.peek().getAttribute(this.umiBarcodeTag);

		ReadCluster c = new ReadCluster(windowSize, cellBC, umiBC);
		while (reads.hasNext()) {
			SAMRecord r = reads.peek();
			boolean added = c.addReadNearestNeighbor(r);
			if (added)
				reads.next(); // consume the read instead of peeking.
			else // read failed to add, this cluster is done.
				return c;
		} // you consumed all the reads in the group, the cluster is done.
		return c;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}

	@Override
	/**
	 * Are there more reads available to cluster?
	 *
	 */
	public boolean hasNext() {
		// if there are more reads, or more groups of reads.
		if (this.iterator.hasNext()) return true;
		if (this.reads!=null)
			return this.reads.hasNext();
		return false;
	}

	@Override
	public void close() {
		CloserUtil.close(iterator);
	}

}
