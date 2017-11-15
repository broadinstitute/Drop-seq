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

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.utils.CountChangingIteratorWrapper;

import htsjdk.samtools.SAMRecord;

//TODO: Better name?
public class UmiIteratorWrapper extends CountChangingIteratorWrapper<SAMRecord> {

	private String cellBarcodeTag;
	private Set<String> cellBarcodeList;
	private String geneExonTag;
	private String strandTag;
	private int readMQ;
	private boolean assignReadsToAllGenes;
	private boolean useStrandInfo;
	private Random rand = new Random();

	/**
	 * Filters and copies reads for generating UMICollections.
	 * Only accepts reads that where the read cell barcode matches a barcode in the list (if not null)
	 * Reads that are marked as secondary or supplementary are rejected
	 * Filters reads based on read map quality, removing reads below that quality
	 * Optionally filters reads where the annotated gene and the strand of the read don't match, or can clone a read and return it multiple times
	 * if the read maps to more than one gene and <assignReadsToAllGenes> is true.
	 * @param cellBarcodeTag The cell barcode BAM tag
	 * @param cellBarcodeList A list of cell barcodes, or null to ignore.  If populated, reads where the cellBarcodeTag matches one of these Strings will be retained
	 * @param geneExonTag The gene/exon tag.
	 * @param strandTag The strand tag
	 * @param readMQ The minimum map quality of a read to be retained.
	 * @param assignReadsToAllGenes Clone SAMRecord for each gene
	 * @param useStrandInfo Only queue SAMRecord if gene strand agrees with read strand.
	 */
	public UmiIteratorWrapper(final Iterator<SAMRecord> underlyingIterator,
                              final String cellBarcodeTag,
                              final Collection<String> cellBarcodeList,
                              final String geneExonTag,
                              final String strandTag,
                              final int readMQ,
                              final boolean assignReadsToAllGenes,
                              final boolean useStrandInfo) {
        super(underlyingIterator);
		this.cellBarcodeTag = cellBarcodeTag;
		this.cellBarcodeList = new HashSet<>(cellBarcodeList);
		this.geneExonTag=geneExonTag;
		this.strandTag= strandTag;
		this.readMQ = readMQ;
		this.assignReadsToAllGenes = assignReadsToAllGenes;
		this.useStrandInfo = useStrandInfo;
	}

    @Override
    protected void processRecord(final SAMRecord r) {
		String cellBC=r.getStringAttribute(cellBarcodeTag);
		String geneList = r.getStringAttribute(this.geneExonTag);

		// if there are cell barcodes to filter on, and this read's cell barcode isn't one of them, then move on to the next read;
		if (this.cellBarcodeList!=null && !cellBarcodeList.contains(cellBC))
			return;
		// if the read doesn't pass map quality, etc, more on.
		if (r.isSecondaryOrSupplementary() || r.getMappingQuality()<this.readMQ || geneList==null)
			return;

		// there's at least one good copy of the read.  Does the read match on strand/gene, or is it assigned to multiple genes?
		String [] genes = geneList.split(",");
		String [] strands = null;

		if (this.useStrandInfo) {
			strands = r.getStringAttribute(strandTag).split(",");
			if (strands==null)
				throw new IllegalStateException("For read [" + r.getReadName()+"] gene tags found [" + geneList +" but no strand info set when strand use was requested.");
			if (genes.length!=strands.length)
				throw new IllegalStateException("For read [" + r.getReadName()+"] gene tags found [" + geneList +"] but a different number of strand tags were set [" + r.getStringAttribute(strandTag) +"]");
		}

		if (this.assignReadsToAllGenes)
			for (int i=0; i<genes.length; i++) {
				String g = genes[i];
				SAMRecord rr = Utils.getClone(r);
				rr.setAttribute(geneExonTag, g);

				// if you use strand info, then the gene has to match the read strand for the read to be added.
				if (useStrandInfo) {
					String geneStrand = strands[i];
					String readStrandString = Utils.strandToString(!r.getReadNegativeStrandFlag());
					if (geneStrand.equals(readStrandString))
						queueRecordForOutput(rr);
					else {
						// rejected read
					}
				} else
					queueRecordForOutput(rr);
			}
		else {
			// pick a random read.
			int randomNum = rand.nextInt((genes.length-1) + 1);
			String g = genes[randomNum];
			r.setAttribute(geneExonTag, g);
            queueRecordForOutput(r);
		}
	}
}
