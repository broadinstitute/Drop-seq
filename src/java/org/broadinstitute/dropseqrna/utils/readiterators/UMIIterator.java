package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.PeekableIterator;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;

public class UMIIterator implements CloseableIterator<UMICollection>  {
	
	private final AggregatedTagOrderIterator atoi;
	private final String geneExonTag;
	private final String cellBarcodeTag;
	private final String molecularBarcodeTag;
	
	/**
	 * Construct an object that generates UMI objects from a BAM file
	 * @param bamFile The BAM file to extract UMIs from
	 * @param geneExonTag The geneExon tag on BAM records
	 * @param cellBarcodeTag The cell barcode tag on BAM records
	 * @param molecularBarcodeTag The molecular barcode tag on BAM records
	 * @param strandTag The strand tag on BAM records
	 * @param maxRecordsInRAM the maximum number of records that can be stored in RAM.
	 * @param readMQ The minimum map quality of the reads 
	 * @param assignReadsToAllGenes Should records tagged with multiple genes be double counted, once for each gene?
	 * @param useStrandInfo should the gene and read strand match for the read to be accepted
	 * @param cellBarcodes The list of cell barcode tag values that match the <cellBarcodeTag> tag on the BAM records.  Only reads with these values will be used.  If set to null, all cell barcodes are used. 
	 */
	public UMIIterator(File bamFile, String geneExonTag, String cellBarcodeTag, String molecularBarcodeTag, String strandTag, int readMQ, boolean assignReadsToAllGenes, boolean useStrandInfo, List<String> cellBarcodes, int maxRecordsInRAM) {
		this.geneExonTag=geneExonTag;
		this.cellBarcodeTag=cellBarcodeTag;
		this.molecularBarcodeTag=molecularBarcodeTag;
		
		// assign the tags in the order you want data sorted.
		List<String> sortingTags = new ArrayList<String>();
		sortingTags = new ArrayList<String>();
		sortingTags.add(geneExonTag);
		sortingTags.add(cellBarcodeTag);
		
		UMIReadProcessor f = new UMIReadProcessor(cellBarcodeTag, cellBarcodes, geneExonTag, strandTag, readMQ, assignReadsToAllGenes, useStrandInfo);
		
		TagOrderIterator toi = new TagOrderIterator(bamFile, sortingTags, f, true);
		this.atoi = new AggregatedTagOrderIterator(toi);
		
	}
	
	/**
	 * Gets the next UMI Collection - all the molecular barcodes and reads on those for a particular gene/cell.
	 * @return Null if there are no reads left in the iterator.  Otherwise, returns a UMICollection.
	 */
	@Override
	public UMICollection next () {
		if (!this.atoi.hasNext()) {
			return null;
		}
		
		Collection<SAMRecord> records = this.atoi.next();
		PeekableIterator<SAMRecord> recordCollectionIter = new PeekableIterator<SAMRecord>(records.iterator());
		
		// the next record is the first of the "batch"
		SAMRecord r = recordCollectionIter.peek();
		String currentGene = r.getStringAttribute(geneExonTag);
		String currentCell = Utils.getCellBC(r, cellBarcodeTag);
		UMICollection umi = new UMICollection(currentCell, currentGene);		
		
		while (recordCollectionIter.hasNext()) {
			// if there's a next read, set up the current gene/cell variables.
			r=recordCollectionIter.next();
			String molecularBarcode = r.getStringAttribute(this.molecularBarcodeTag);
			umi.incrementMolecularBarcodeCount(molecularBarcode);
						
		}
		// I ran out of reads
		recordCollectionIter.close();
		return (umi);
	}
	
	@Override
	public void remove() {
		this.atoi.remove();
	}

	@Override
	public void close() {
		this.atoi.close();
	}

	@Override
	public boolean hasNext() {
		return this.atoi.hasNext();
	}
	
	
}
