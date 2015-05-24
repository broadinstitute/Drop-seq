package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.editdistance.EDUtils;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        usage = "Measures the digital expression of a library.  " +
			"Method: 1) For each gene, find the molecular barcodes on the exons of that gene.  " +
			"2) Determine how many HQ mapped reads are assigned to each barcode.  " +
			"3) Collapse barcodes by edit distance.  " +
			"4) Throw away barcodes with less than threshold # of reads. " +
			"5) Count the number of remaining unique molecular barcodes for the gene." +
			"This program requires a tag for what gene a read is on, a molecular barcode tag, and a exon tag.  The exon and gene tags may not be present on every read." +
			"When filtering the data for a set of barcodes to use, the data is filtered by ONE of the following methods (and if multiple params are filled in, the top one takes precidence):" +
			"1) CELL_BC_FILE, to filter by the some fixed list of cell barcodes" + 
			"2) MIN_NUM_GENES_PER_CELL " + 
			"3) MIN_NUM_TRANSCRIPTS_PER_CELL " +
			"4) NUM_CORE_BARCODES " +
			"5) MIN_NUM_READS_PER_CELL",
        usageShort = "Calculate Digital Expression",
        programGroup = DropSeq.class
)

public class DigitalExpression extends DGECommandLineBase {

	private static final Log log = Log.getInstance(DigitalExpression.class);

	@Option(doc="A summary of the digital expression output, containing 3 columns - the cell barcode, the #genes, and the #transcripts.", optional=true)
	public File SUMMARY;
	
	@Option(doc="Output number of reads instead of number of unique molecular barcodes.", optional=true)
	public boolean OUTPUT_READS_INSTEAD=false;
	
	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file of DGE Matrix.  Genes are in rows, cells in columns.  The first column contains the gene name. This supports zipped formats like gz and bz2.")
	public File OUTPUT;
	
	private boolean OUTPUT_EXPRESSED_GENES_ONLY=false;
	
	@Override
	/**
	 * This is a revision of the original DGE code to implement a more complicated state machine in the main loop and in exchange get rid of the batch system.
	 * This change lets you store the counts of molecular barcodes for gene/cell instead of the reads for one gene/cell, which can save huge amounts of memory when
	 * cells have millions of reads on housekeeping genes. 
	 */
	protected int doWork() {
		
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		if (this.SUMMARY!=null) IOUtil.assertFileIsWritable(this.SUMMARY);
		
		// boolean check = new Utils().validateGetCellBarcodeListParams(this.INPUT, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, 
		//		this.GENE_TAG, this.EXON_TAG, this.CELL_BC_FILE, this.READ_MQ, this.MIN_NUM_TRANSCRIPTS_PER_CELL, 
		//		this.MIN_NUM_GENES_PER_CELL, this.MIN_NUM_READS_PER_CELL, this.NUM_CORE_BARCODES);
		
		List<String> cellBarcodes=new BarcodeListRetrieval().getCellBarcodes(this.INPUT, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, 
				this.GENE_EXON_TAG, this.STRAND_TAG, this.CELL_BC_FILE, this.READ_MQ, this.MIN_NUM_TRANSCRIPTS_PER_CELL, 
				this.MIN_NUM_GENES_PER_CELL, this.MIN_NUM_READS_PER_CELL, this.NUM_CORE_BARCODES, this.EDIT_DISTANCE, this.MIN_BC_READ_THRESHOLD, USE_STRAND_INFO);
		
		if (cellBarcodes.isEmpty()) {
			log.error("Running digital expression without somehow selecting a set of barcodes to process no longer supported.");
			return (1);
		} else {
			log.info("Calculating digital expression for [" + cellBarcodes.size()+ "] cells.");
			digitalExpression(cellBarcodes);
		}
			
		return 0;
	}
	

	
	private void digitalExpression(List<String> cellBarcodes) {
		PrintStream out = new PrintStream(IOUtil.openFileForWriting(OUTPUT));
		
		writeHeader(out, cellBarcodes);

		UMIIterator umiIterator = new UMIIterator(this.INPUT, this.GENE_EXON_TAG, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.STRAND_TAG, 
		this.READ_MQ, true, this.USE_STRAND_INFO, cellBarcodes, this.MAX_RECORDS_IN_RAM);
		
		String gene = null;
		Map<String, Integer> countMap = new HashMap<String, Integer>();
		Map<String, DESummary> summaryMap = initializeSummary(cellBarcodes);
		
		
		UMICollection batch;
		while ((batch=umiIterator.next())!=null) {
			if (batch==null || batch.isEmpty()){
				continue;	
			}
			String currentGene = batch.getGeneName();
			// if just starting the loop
			if (gene==null) gene=currentGene;
			// if the gene is the same, you're still gathering expression on that gene.
			if (gene.equals(currentGene)) {
				
				if (this.RARE_UMI_FILTER_THRESHOLD>0) batch.filterByUMIFrequency(this.RARE_UMI_FILTER_THRESHOLD);
				int molBCCount = batch.getDigitalExpression(this.MIN_BC_READ_THRESHOLD, this.EDIT_DISTANCE, this.OUTPUT_READS_INSTEAD);
				countMap.put(batch.getCellBarcode(), molBCCount);
			}
			// you've gathered all the data for the gene, write it out and start on the next.
			if (!gene.equals(currentGene)) {
				writeStats (gene, countMap, cellBarcodes, out);
				addToSummary(countMap, summaryMap);
				countMap.clear();
				// start the next gene
				if (this.RARE_UMI_FILTER_THRESHOLD>0) batch.filterByUMIFrequency(this.RARE_UMI_FILTER_THRESHOLD);
				int molBCCount = batch.getDigitalExpression(this.MIN_BC_READ_THRESHOLD, this.EDIT_DISTANCE, this.OUTPUT_READS_INSTEAD);
				countMap.put(batch.getCellBarcode(), molBCCount);
				
				gene=currentGene;
			}
		}
		// write out remainder
		if (countMap.isEmpty()==false) {
			writeStats (gene, countMap, cellBarcodes, out);
			addToSummary(countMap, summaryMap);
		}
		out.close();
		if (this.SUMMARY!=null) {
			writeSummary(summaryMap.values(), this.SUMMARY);
		}
		
		
	}
	 
	
	
	/**
	 * Collapses a bunch of strings by the edit distance.
	 * If edit distance computations indicate it's greater than threshold edit distance, then the threshold is returned.
	 * This is to avoid hard work on indel calculations when edit distance between two strings is high.
	 * You can safely set threshold to be 3 * edit distance.
	 * @param barcodes
	 * @param editDistance
	 * @param threshold
	 * @return
	 */
	public ObjectCounter <String> collapseByEditDistance (ObjectCounter<String> barcodes, int editDistance, int threshold) {
		// map the barcode to the object so I can look up counts
		
		
		ObjectCounter <String> result = new ObjectCounter<String>();
		List<String> barcodeList = barcodes.getKeysOrderedByCount(true);
		
		// short circuit for ED=0
		if (this.EDIT_DISTANCE==0) {
			for (String barcode: barcodeList) {
				int count=barcodes.getCountForKey(barcode);
				result.setCount(barcode, count);
			}
			return (result);
		}
		
		while (barcodeList.isEmpty()==false) {
			String b = barcodeList.get(0);
			barcodeList.remove(b);
			// this is still the "old" single core version.  Molecular barcode counts are small, so this may be ok.
			Set<String> closeBC = EDUtils.getInstance().getStringsWithinEditDistanceWithIndel(b,barcodeList, editDistance);
			barcodeList.removeAll(closeBC);
			// for counting.
			closeBC.add(b);
			int totalCount = 0;
			for (String bc: closeBC) {
				int count = barcodes.getCountForKey(bc);
				totalCount+=count;
			}
			result.setCount(b, totalCount);
		}		
		return (result);
	}
		
	
	private void writeStats (String gene, Map<String, Integer> countMap, List<String> cellBarcodes, PrintStream out) {
		
		int totalCount=0;
		List<String> line = new ArrayList<String>(cellBarcodes.size()+1);
		line.add(gene);
		for (String b: cellBarcodes) {
			Integer count = countMap.get(b);
			String v = "0";
			if (count!=null) {
				totalCount+=count;
				v = count.toString();
			}
			line.add(v);
		}
		if (OUTPUT_EXPRESSED_GENES_ONLY & totalCount==0) return;
		String h = StringUtils.join(line, "\t");
		out.println(h);
		//OutputWriterUtil.writeResult(h, out);	
	}
	
	
	private void writeHeader(PrintStream out, List<String> cellBarcodes) {
		List<String> header = new ArrayList<String>(cellBarcodes.size()+1);
		header.add("GENE");
		for (String c: cellBarcodes) {
			header.add(c);
		}
		String h = StringUtils.join(header, "\t");
		out.println(h);
	}
	
	public class DESummary extends MetricBase {
		
		public String CELL_BARCODE;
		public int NUM_GENES;
		public int NUM_TRANSCRIPTS;
		
		public DESummary (String cellBarcode) {
			this.CELL_BARCODE=cellBarcode;
			this.NUM_GENES=0;
			this.NUM_TRANSCRIPTS=0;
		}
		
	}
	
	static final Comparator<DESummary> TRANSCRIPT_ORDER_DESCENDING =  new Comparator<DESummary>() {
		public int compare(DESummary e1, DESummary e2) {
			return (e1.NUM_TRANSCRIPTS > e2.NUM_TRANSCRIPTS ? -1 : e1.NUM_TRANSCRIPTS == e2.NUM_TRANSCRIPTS ? 0 : 1);
		}
};
	
	public Map<String, DESummary> initializeSummary(Collection<String> cellBarcodes) {
		Map<String, DESummary> map = new HashMap<String, DESummary>();
		
		for (String s: cellBarcodes) {
			DESummary des = new DESummary(s);
			map.put(s, des);
		}
		return (map);
	}
	
	public Map<String, DESummary> addToSummary(Map<String, Integer> countMap, Map<String, DESummary> summaryMap) {
		for (String cellBC: countMap.keySet()) {
			DESummary sum = summaryMap.get(cellBC);
			// for genes, it doesn't matter what the count is as long as it's > 0.  Increment by 1.
			sum.NUM_GENES++;
			// for transcripts, increment by the count.
			sum.NUM_TRANSCRIPTS+=countMap.get(cellBC);
		}
		return (summaryMap);
	}
	
	public void writeSummary(Collection<DESummary> summaryCollection, File outFile) {
		MetricsFile<DESummary, Integer> out = new MetricsFile<DESummary, Integer>();
		List<DESummary> sc = new ArrayList<DESummary>(summaryCollection);
		Collections.sort(sc, DigitalExpression.TRANSCRIPT_ORDER_DESCENDING);
		for (DESummary z: sc) {
			out.addMetric(z);
		}
		out.write(outFile);
	}
	

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new DigitalExpression().instanceMain(args));
	}
	
	
}
