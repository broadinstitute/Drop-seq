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

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * Writes a variety of reports from a streamed in set of ReadClusters.
 *
 * IntervalReport: for each cell/molecular barcode, write out the cell / umi / interval / #reads forward / #reads reverse
 * BedReport: the chromosome/start/end/name/strand.  The name is the cell barcode + molecular barcode colon separated.
 * The strand is the strand of the majority of reads.
 * IntervalMinimumDistanceReport for each cell/molecular barcode, output the minimum distance between the intervals of the read islands.
 * @author nemesh
 *
 */
public class ReadClusterReportWriter {

	private static final Log log = Log.getInstance(ReadClusterReportWriter.class);

	private final PrintStream intervalReport;
	private final PrintStream bedReport;
	private final PrintStream minDistanceReport;
	private IntervalAggregator ia;
	private final Set<String> cellBarcodes;

	/**
	 * Initialize the report writer.  All reports are optional.
	 * @param intervalReportFile for each cell/molecular barcode, write out the cell / umi / interval / #reads forward / #reads reverse
	 * @param bedReportFile the chromosome/start/end/name/strand.  The name is the cell barcode + molecular barcode colon separated.
	 * The strand is the strand of the majority of reads.
	 * @param cellBarcodes If this parameter is set, only report cells in this set.
	 * @param intervalMinimumDistanceReportFile for each cell/molecular barcode, output the minimum distance between the intervals of the read islands.
	 */
	public ReadClusterReportWriter (final File intervalReportFile, final File bedReportFile, final File intervalMinimumDistanceReportFile, final Collection<String> cellBarcodes) {
		this.intervalReport = getPrintStream(intervalReportFile);
		this.bedReport = getPrintStream(bedReportFile);
		this.minDistanceReport = getPrintStream(intervalMinimumDistanceReportFile);
		writeIntervalReportHeader(this.intervalReport);
		writeBedReportHeader(this.bedReport);
		writeMinDistanceReportHeader(this.minDistanceReport);
		this.ia=new IntervalAggregator();
		if (cellBarcodes==null)
			this.cellBarcodes=new HashSet<String>(0);
		else
			this.cellBarcodes=new HashSet<String>(cellBarcodes);
	}

	private static PrintStream getPrintStream (final File f) {
		if (f==null) return null;

		return new ErrorCheckingPrintStream(IOUtil.openFileForWriting(f));
	}

	public void addReadCluster (final ReadCluster c) {
		// short circuit - if there's a list of cell barcodes and this isn't in it, don't report the cell barcode.
		if (!hasCellBarcode(c.getCellBarcode()))
			return;
		writeIntervalReport(c);
		writeBedReport(c);
		// calculate minimum distance between intervals.
		this.ia.add(c);
	}


	private static void writeIntervalReportHeader (final PrintStream p) {
		if (p==null) return;
		List<String> h = new ArrayList<String>();
		h.add("cell");
		h.add("umi");
		h.add("chr");
		h.add("start");
		h.add("end");
		h.add("length");
		h.add("reads_forward");
		h.add("reads_reverse");
		String line = StringUtils.join(h, "\t");
		p.println(line);
	}

	void writeIntervalReport (final ReadCluster c) {
		if (c.isUnmapped()) return;
		List<String> h = new ArrayList<String>();
		h.add(c.getCellBarcode());
		h.add(c.getUMIBarcode());
		Interval i = c.getClusterInterval();
		h.add(i.getContig());
		h.add(Integer.toString(i.getStart()));
		h.add(Integer.toString(i.getEnd()));
		h.add(Integer.toString(i.getEnd()-i.getStart()+1));
		h.add(Integer.toString(c.getNumReadsPositiveStrand()));
		h.add(Integer.toString(c.getNumReadsNegativeStrand()));
		String line = StringUtils.join(h, "\t");
		this.intervalReport.println(line);
	}

	private static void writeBedReportHeader(final PrintStream p) {
		if (p==null) return;
		String header = "track name=\"SpermSeq ReadClusters\"";
		p.println(header);

	}
	/**
	 * Bed file format: chrom, start (base 0), end, name, score, strand
	 * @param c
	 */
	void writeBedReport (final ReadCluster c) {
		if (c.isUnmapped()) return;
		Interval i = c.getClusterInterval();
		String name = c.getCellBarcode()+":"+c.getUMIBarcode();
		int pos=c.getNumReadsPositiveStrand();
		int neg = c.getNumReadsNegativeStrand();
		int majority=pos;
		String strand="+";
		if (neg>pos) {
			majority=neg;
			strand="-";
		}

		String[] l = {i.getContig(), Integer.toString(i.getStart()-1), Integer.toString(i.getEnd()), name, Integer.toString(majority), strand};
		String line = StringUtils.join(l, "\t");
		this.bedReport.println(line);
	}

	private static void writeMinDistanceReportHeader(final PrintStream p) {
		if (p==null) return;
		List<String> h = new ArrayList<String>();
		h.add("cell");
		h.add("umi");
		h.add("minimum_cluster_distance");
		String line = StringUtils.join(h, "\t");
		p.println(line);
	}

	/**
	 * Only call this report after all clusters have been added.
	 * @param c
	 */
	public void writeMinDistanceReport () {
		Map<String, Set<String>> barcodes = this.ia.getCellUMIs();
		for (String cellBarcode: barcodes.keySet()) {
			Set<String> umiBarcodes = barcodes.get(cellBarcode);
			for (String umiBarcode: umiBarcodes) {
				int distance = ia.getMinimumDistance(cellBarcode, umiBarcode);
				if (distance==Integer.MAX_VALUE) {
					String [] l = {cellBarcode, umiBarcode, "NA"};
					String line = StringUtils.join(l, "\t");
					this.minDistanceReport.println(line);
				} else {
					String [] l = {cellBarcode, umiBarcode, Integer.toString(distance)};
					String line = StringUtils.join(l, "\t");
					this.minDistanceReport.println(line);
				}
			}
		}
	}

	/**
	 * Keeps track of the last seen interval and distance between intervals for a cell/umi
	 *
	 */
	class IntervalAggregator {
		// key: cell barcode : umi barcode
		// value = minimum distance
		private Map<String, Integer> distanceMap = new HashMap<String,Integer>();
		private Map<String, Interval> intervalMap = new HashMap<String,Interval>();

		public void add (final ReadCluster c) {
			if (c.isUnmapped()) return;
			String key = getKey(c.getCellBarcode(), c.getUMIBarcode());
			Interval previousInterval = this.intervalMap.get(key);
			Interval currentInterval = c.getClusterInterval();
			Integer dist = Integer.MAX_VALUE;

			// if there's no previous interval, set up the current distance as the maximum.
			if (previousInterval==null)
				distanceMap.put(key, dist);

			if (previousInterval!=null) {
				// if the distance is on a separate chromosome, the distance should be something large like Integer.MAX

				if (currentInterval.getContig().equals(previousInterval.getContig()))
					dist = currentInterval.getStart() - (previousInterval.getEnd()+1);

				if (dist<0)
					log.warn("Overlapping ReadCluster detected: "+ c.toString() +" previous interval " + previousInterval.toString());

				// update distance compared to old distance if smaller.
				int oldDist = distanceMap.get(key);
				if (dist< oldDist)
					distanceMap.put(key, dist);
			}
			intervalMap.put(key, currentInterval);
		}

		public int getMinimumDistance (final String cellBarcode, final String umiBarcode) {
			String key = getKey(cellBarcode, umiBarcode);
			return distanceMap.get(key);
		}

		private String getKey (final String cellBarcode, final String umiBarcode) {
			return cellBarcode+":"+umiBarcode;
		}

		public Map<String,Set<String>> getCellUMIs () {
			Map<String, Set<String>> result = new HashMap<String, Set<String>>();
			for (String key : this.intervalMap.keySet()) {
				String [] vals = key.split(":");
				Set<String> umiSet = result.get(vals[0]);
				if (umiSet==null)
					umiSet=new HashSet<String>();
				umiSet.add(vals[1]);
				result.put(vals[0], umiSet);
			}
			return (result);
		}
	}

	private boolean hasCellBarcode (final String cellBC) {
		if (this.cellBarcodes.isEmpty()) return true;
		return this.cellBarcodes.contains(cellBC);
	}

}
