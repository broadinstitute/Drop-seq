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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.utils.CountChangingIteratorWrapper;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;

import java.util.*;

public class SNPUMICellReadIteratorWrapper extends CountChangingIteratorWrapper<SAMRecord> {

	private String cellBarcodeTag;
	private Set<String> cellBarcodeList;
	private String geneTag;
	private final String snpTag;
	private final OverlapDetector<Interval> snpIntervals;


	/**
	 * Filters and copies reads for generating SNPCellUMIBasePileups.
	 * Filters out reads in which the read cell barcode does not match a barcode in the list (if list is not null)
	 * Reads that are marked as secondary or supplementary are filtered out
	 * Filters reads based on read map quality, removing reads below that quality
	 * Optionally filters reads where the annotated gene and the strand of the read don't match, or can clone a read and return it multiple times
	 * if the read maps to more than one gene and <assignReadsToAllGenes> is true.
	 * @param underlyingIterator Source of reads
	 * @param cellBarcodeTag The cell barcode BAM tag
	 * @param cellBarcodeList A list of cell barcodes, or null to ignore.  If populated, reads where the cellBarcodeTag matches one of these Strings will be retained
	 * @param geneTag The gene/exon tag.
	 * @param strandTag The strand tag
	 * @param readMQ The minimum map quality of a read to be retained.
	 */
	public SNPUMICellReadIteratorWrapper(final Iterator<SAMRecord> underlyingIterator,
                                         final IntervalList snpIntervals,
                                         final String cellBarcodeTag,
                                         final Collection<String> cellBarcodeList,
                                         final String geneTag,
                                         final String snpTag,
                                         final int readMQ) {
        super(underlyingIterator);
		this.cellBarcodeTag = cellBarcodeTag;
		this.cellBarcodeList = new HashSet<String>(cellBarcodeList);
		this.geneTag=geneTag;
		this.snpTag=snpTag;


		// construct OverlapDetector
		OverlapDetector<Interval> od = new OverlapDetector<>(0, 0);
		od.addAll(snpIntervals.getIntervals(), snpIntervals.getIntervals());
		this.snpIntervals=od;
	}

    @Override
    protected void processRecord(final SAMRecord r) {
        String cellBC=r.getStringAttribute(cellBarcodeTag);
        String geneList = r.getStringAttribute(this.geneTag);
        List<SAMTagAndValue> allAttributes= r.getAttributes();

        // if there are cell barcodes to filter on, and this read's cell barcode isn't one of them, then move on to the next read;
        if (this.cellBarcodeList!=null && !cellBarcodeList.contains(cellBC))
			return;
        // if the read have any genes.
        if (geneList==null)
			return;
        processGene(r);
    }

	/**
	 * Check if a read overlaps any SNPs in the OverlapDetector.  Tag reads with SNPs.
	 * If more than 1 SNP tags a read, make a read for each SNP.
	 * Simplified since data goes through GeneFunctionIteratorWrapper to take care of how reads/genes interact.
	 */
	private void processSNP (final SAMRecord r) {
		List<AlignmentBlock> blocks = r.getAlignmentBlocks();

		Collection<String> snps = new HashSet<>();

		for (AlignmentBlock b: blocks) {
			int start = b.getReferenceStart();
			int end = start + b.getLength() -1;

			Interval i = new Interval(r.getReferenceName(), start, end);
			Collection<Interval> overlaps = this.snpIntervals.getOverlaps(i);
			for (Interval o: overlaps)
				snps.add(IntervalTagComparator.toString(o));
		}

		// 1 read per SNP.
		for (String snp:snps) {
			SAMRecord rr = Utils.getClone(r);
			rr.setAttribute(this.snpTag, snp);
			queueRecordForOutput(rr);
		}
	}


	/**
	 * For a read, check and see if the read maps to more than 1 gene, and if the gene matches the strand (optional, controlled by useStrandInfo flag).
	 * If useStrandInfo is true and the read strand does not match the gene strand, discard the read.
	 * If useStrandInfo is true and the read strand matches the gene strand, keep the read.
	 * If useStrandInfo is false, keep the read
	 * If there is more than 1 gene tagging the read, make a copy of the read for each gene.  This also uses the strand option.
	 */
	private void processGene (final SAMRecord r) {
		String geneList = r.getStringAttribute(this.geneTag);
		String [] genes = geneList.split(",");


		// if there's just 1 gene (the common case) avoid cloning the read and use strand info if required to accept/reject the read.
		if (genes.length==1) {
			processSNP(r);
			return;
		}

		// for more than 1 case, clone the subsequent reads, update the gene to the new
		for (String g : genes) {
			SAMRecord rr = Utils.getClone(r);
			rr.setAttribute(geneTag, g);
			processSNP(rr);
		}
	}

}
