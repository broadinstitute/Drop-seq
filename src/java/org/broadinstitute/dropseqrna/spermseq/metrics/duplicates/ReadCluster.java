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
import htsjdk.samtools.util.Interval;

import java.util.ArrayList;
import java.util.List;

public class ReadCluster {

	private List<SAMRecord> reads = new ArrayList<SAMRecord>();
	private final int windowSize;
	private Integer clusterStart=null;
	private Integer clusterEnd=null;
	private String contig=null;
	private Boolean unmapped=null;
	private final String cell;
	private final String umi;

	/**
	 * Generate an empty ReadCluster with a set window size.  All reads must fit within this size.
	 *
	 * @param windowSize The size in base pairs that all reads will be contained in.
	 */
	public ReadCluster(final int windowSize, final String cellBarcode, final String umiBarcode) {
		this.cell=cellBarcode;
		this.umi=umiBarcode;
		this.windowSize=windowSize;
	}

	/**
	 * True if the cluster contains mapped reads, false if the cluster contains unmapped reads, or null if no reads have yet been added.
	 * @return
	 */
	public Boolean isUnmapped() {
		return this.unmapped;
	}

	public String getCellBarcode () {
		return this.cell;
	}

	public String getUMIBarcode () {
		return this.umi;
	}

	public List<SAMRecord> getReads () {
		return this.reads;
	}

	public Interval getClusterInterval () {
		if (this.isUnmapped()) return null;
		return new Interval(contig, clusterStart, clusterEnd);
	}

	public int getNumReadsPositiveStrand () {
		int count =0;
		for (SAMRecord r: this.reads)
			if (!r.getReadNegativeStrandFlag())
				count++;
		return count;
	}

	public int getNumReadsNegativeStrand () {
		int count =0;
		for (SAMRecord r: this.reads)
			if (r.getReadNegativeStrandFlag())
				count++;
		return count;
	}


	/**
	 * Attempts to add this read to the ReadCluster.
	 *
	 * @param r The read to add
	 * @return True if the read fits inside the bounds of this cluster as defined by the windowSize and is successfully added, otherwise false.
	 * A read on a different chromosome than the current collection will also be rejected.
	 */
	public boolean addReadInWindow (final SAMRecord r) {
		boolean readUnmapped=r.getReadUnmappedFlag();

		// short circuit for an unmapped read cluster.
		// if there are no reads, the next read determines if the read is mapped or not.
		if (this.reads.isEmpty())
			this.unmapped=readUnmapped;

		// if the cluster is unmapped
		if (this.unmapped)
			if (readUnmapped) { // if unmapped cluster & read, accept the read and return true.
				this.reads.add(r);
				return true;
			} else // if unmapped cluster, but mapped read, reject the read.
				return false;

		int s = r.getAlignmentStart();
		int e = r.getAlignmentEnd();

		// if the start and end are unset, set them.
		if (clusterStart==null || clusterEnd ==null) {
			this.contig=r.getContig();
			clusterStart=s;
			clusterEnd=e;
		}

		// test the bounds.
		if (!this.contig.equals(r.getContig()))
			return false; // the contig is different from reads already in the cluster.
		if (e>=clusterEnd) {
			int len=e-clusterStart+1;
			if (len>windowSize) return false; // this end position is too far.
		}
		if (s<=clusterStart) {
			int len=clusterEnd-s+1;
			if (len>windowSize) return false; // this start position is too far.
		}
		// you're within the bounds, add the read and modify the bounds.
		if (s<clusterStart)
			clusterStart=s;
		if (e>clusterEnd)
			clusterEnd=e;

		// add the read.
		this.reads.add(r);
		return true;
	}

	/**
	 * Attempts to add this read to the ReadCluster.
	 * This greedily adds reads to a cluster until the next read is more than some distance away, and then closes the cluster and starts a new one with the next read.
	 * In this add method, the window size is used determine the maximum distance a read can be to the next read before the cluster is closed.
	 * @param r The read to add
	 * @return True if the read fits inside the bounds of this cluster as defined by the windowSize and is successfully added, otherwise false.
	 * A read on a different chromosome than the current collection will also be rejected.
	 */
	//TODO: I think this may not be safe for split read data.  Might need to switch to alignment blocks if we ever want to use this class in that context.
	public boolean addReadNearestNeighbor (final SAMRecord r) {
		boolean readUnmapped=r.getReadUnmappedFlag();

		// short circuit for an unmapped read cluster.
		// if there are no reads, the next read determines if the read is mapped or not.
		if (this.reads.isEmpty())
			this.unmapped=readUnmapped;

		// if the cluster is unmapped
		if (this.unmapped)
			if (readUnmapped) { // if unmapped cluster & read, accept the read and return true.
				this.reads.add(r);
				return true;
			} else // if unmapped cluster, but mapped read, reject the read.
				return false;

		int s = r.getAlignmentStart();
		int e = r.getAlignmentEnd();

		// if the start and end are unset, set them.
		if (clusterStart==null || clusterEnd ==null) {
			this.contig=r.getContig();
			clusterStart=s;
			clusterEnd=e;
		}

		// test the bounds.
		if (!this.contig.equals(r.getContig()))
			return false; // the contig is different from reads already in the cluster.

		// length from the end of the cluster to this read's start.
		int lenEnd=Math.min(Math.abs(s-this.clusterEnd), Math.abs(e-this.clusterEnd));
		// length from this read to the start of the cluster, covers the 3' side.
		int lenStart=Math.min(Math.abs(s-this.clusterStart), Math.abs(e-this.clusterStart));
		int maxLen=Math.min(lenEnd, lenStart);

		// if too far away, return false.
		if (maxLen>this.windowSize)
			return false;

		// you're within the bounds, add the read and modify the bounds.
		if (s<clusterStart)
			clusterStart=s;
		if (e>clusterEnd)
			clusterEnd=e;

		// add the read.
		this.reads.add(r);
		return true;
	}



	@Override
	public String toString () {
		StringBuilder b = new StringBuilder();
		Interval i = this.getClusterInterval();
		b.append("cell [" + this.cell + "] ");
		b.append("umi [" + this.umi +"] ");
		b.append("Window Size [" + this.windowSize +"] ");
		if (!this.isUnmapped()) {
			b.append("Interval [" + i+"] ");
			b.append("Cluster Length [" + i.length()+"] ");
		}
		b.append("Reads [" + this.reads.size()+"]");
		return b.toString();
	}
}
