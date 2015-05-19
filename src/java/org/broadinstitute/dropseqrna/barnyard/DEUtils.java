package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.BAMTagComparator;

public class DEUtils {

	private static final Log log = Log.getInstance(DEUtils.class);
	private ProgressLogger progress = new ProgressLogger(log, 100000);
	
	private Integer MAX_RECORDS_IN_RAM;
	private String DEFAULT_CELL_BARCODE="DEFAULT";
	
	private Random rand = new Random();
	
	public DEUtils () {
		MAX_RECORDS_IN_RAM = SAMFileWriterImpl.getDefaultMaxRecordsInRam();
	}
	
	public DEUtils (int maxRecordsInRam) {
		MAX_RECORDS_IN_RAM = maxRecordsInRam;
	}
	
	/**
	 * Gets the next UMI Collection - all the molecular barcodes and reads on those for a particular gene/cell.
	 * @param iter A peekable iterator of SAM Records.
	 * @return Null if there are no reads left in the iterator, or a UMICollection.
	 */
	public UMICollection getUMICollection (PeekableIterator<SAMRecord> iter, String geneExonTag, String cellBarcodeTag, String molecularBCTag) {
		if (iter==null || iter.hasNext()==false) return null;
		
		// the next record is the first of the "batch"
		SAMRecord r = iter.peek();
		String currentGene = r.getStringAttribute(geneExonTag);
		String currentCell = getCellBC(r, cellBarcodeTag);
		UMICollection umi = new UMICollection(currentCell, currentGene);
		
		
		while (iter.hasNext()) {
			// if there's a next read, set up the current gene/cell variables.
			r=iter.peek();
			String nextGene = r.getStringAttribute(geneExonTag);
			String nextCell = getCellBC(r, cellBarcodeTag);
				
							
			// you're done with this cell/gene, or you're out of data, so return the UMI
			if ((!nextGene.equals(currentGene) || !nextCell.equals(currentCell))) {
				return (umi); 
			} else {  // you're in the same cell/gene as last read, add the record.
				r=iter.next();
				this.progress.record(r);
				String molecularBarcode = r.getStringAttribute(molecularBCTag);
				umi.incrementMolecularBarcodeCount(molecularBarcode);
			}			
		}
		
		// I ran out of reads
		return (umi);
	}
	
	
	/**
	 * Implements a custom iteration over the bamFile to sort the reads in cell and gene order for digital expression.
	 * This sorts by the GENE_EXON tag (so genes are in order) followed by the cell tag (so all results for a gene are clumped together, but cells can be processed one at a time.) 
	 * @param inFile
	 * @param assignReadsToAllGenes If true, reads that overlap multiple genes will be placed into the iterator once per gene.  If false, the read will be assigned to a random gene
	 * out of the available options.
	 * @cellBarcodes a list of cell barcodes to retain.  If set to null, this is ignored.
	 * @return An iterator of reads, sorted by the GENE_EXON_TAG then the CELL_BARCODE_TAG.
	 */
	CloseableIterator<SAMRecord> getReadsInTagOrder (File inFile, String geneExonTag, String strandTag, String cellBarcodeTag, int readMQ, boolean assignReadsToAllGenes, boolean useStrandInfo, List<String> cellBarcodes) {
		// set up a set for speed.
		Set<String> cellBCSet=null;
		if (cellBarcodes!=null) {
			cellBCSet=new HashSet<String>(cellBarcodes);
		}
		
		SamReader reader = SamReaderFactory.makeDefault().open(inFile);
		// the list of tags to sort in order.
		List<String> l = new ArrayList<String>();
		l.add(geneExonTag);
		l.add(cellBarcodeTag);
		
		SAMSequenceDictionary dict= reader.getFileHeader().getSequenceDictionary();
		List<SAMProgramRecord> programs =reader.getFileHeader().getProgramRecords();
		
		final SAMFileHeader writerHeader = new SAMFileHeader();
        // writerHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        writerHeader.setSequenceDictionary(dict);
        for (SAMProgramRecord spr : programs) {
        	writerHeader.addProgramRecord(spr);
        }
		SortingCollection<SAMRecord> alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
	            new BAMRecordCodec(writerHeader), new BAMTagComparator(l),
	                MAX_RECORDS_IN_RAM);
		
		log.info("Reading in records for TAG name sorting");
		
		ProgressLogger prog = new ProgressLogger(log);
		
		List<SAMRecord> tempList = new ArrayList<SAMRecord>(10);
		int numReadsAdded=0;
		for (SAMRecord r: reader) {
			prog.record(r);
			String cellBC=r.getStringAttribute(cellBarcodeTag);
			// if there are cell barcodes to filter on, and this read's cell barcode isn't one of them, then move on to the next read;
			if (cellBCSet!=null && !cellBCSet.contains(cellBC)) {
				continue;
			}
			
			tempList=getSAMRecords(r, tempList, geneExonTag, strandTag, readMQ, assignReadsToAllGenes, useStrandInfo);
			if (tempList.size()>0) numReadsAdded++;
			for (SAMRecord rr: tempList) {
				alignmentSorter.add(rr);
			}
			
		}
		log.info("Added " + numReadsAdded + " to iterator out of " +prog.getCount());
		
		if (numReadsAdded==0) log.warn("The number of reads added to the iterator was 0.  Either you have all low map quality data, or there's no geneExon or cellBarcode tags on the reads");
		CloserUtil.close(reader);
		CloseableIterator<SAMRecord> result = alignmentSorter.iterator();
		log.info("Sorting finished.");
		return (result);
	}
	
	
	
	
	public static String strandToString(boolean strand) {
		if (strand) return "+";
		return "-";
	}

	/**
	 * 
	 * @param r The read to possibly add to the iterator
	 * @param tempList The list to add reads to
	 * @param assignReadsToAllGenes If true, reads that overlap multiple genes will be placed into the iterator once per gene.  If false, the read will be assigned to a random gene
	 * out of the available options.
	 * @return the tempList, filled with 0 or more instances of the read.
	 */
	
	private List<SAMRecord> getSAMRecords (SAMRecord r, List<SAMRecord> tempList, String geneExonTag, String strandTag, int readMQ, boolean assignReadsToAllGenes, boolean useStrandInfo) {
		tempList.clear();
		String geneList = r.getStringAttribute(geneExonTag);
		
		if (r.isSecondaryOrSupplementary() || r.getMappingQuality()<readMQ || geneList==null) {
			return (tempList);
		}
		
		// there's at least one good copy of the read.
		String [] genes = geneList.split(",");
		String [] strands = null;
		
		if (useStrandInfo) {
			strands = r.getStringAttribute(strandTag).split(",");
		}
		
		if (assignReadsToAllGenes) {
			for (int i=0; i<genes.length; i++) {
				String g = genes[i];
				SAMRecord rr = getClone(r);
				rr.setAttribute(geneExonTag, g);
				
				// if you use strand info, then the gene has to match the read strand for the read to be added.
				if (useStrandInfo) {
					String geneStrand = strands[i];
					String readStrandString = DEUtils.strandToString(!r.getReadNegativeStrandFlag());
					if (geneStrand.equals(readStrandString)) {
						tempList.add(rr);
					} else {
						// rejected read
						// log.info("Read wrong strand");
					}
				} else { // if you don't use strand info, add the read to all genes.
					tempList.add(rr);	
				}
			}
		} 
		else {
			// pick a random read.
			int randomNum = rand.nextInt((genes.length-1) + 1);
			String g = genes[randomNum];
			r.setAttribute(geneExonTag, g);
		}
		return (tempList);
		
	}
	
	private SAMRecord getClone (SAMRecord r) {
		SAMRecord rr=null;
		try {
			rr = (SAMRecord) r.clone();
		} catch (CloneNotSupportedException e) {
			log.info("This should never happen.  SAMRecord can't be cloned?");
		}
		return (rr);
	}
	
	private String getCellBC (SAMRecord r, String cellBCTag) {
		String currentCell = r.getStringAttribute(cellBCTag);
		if (currentCell==null) return (DEFAULT_CELL_BARCODE);
		return (currentCell);
	}
		
}
