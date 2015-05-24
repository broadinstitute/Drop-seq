package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;

import java.io.BufferedWriter;
import java.io.File;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.OutputWriterUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        usage = "For each gene, count the number of times each molecular barcode is observed [The UMI]" +
			"Similar to digital expression, reads are filtered on map quality, and must overlap exons as well as genes. "+ 
			"This program requires a tag for what gene a read is on, a molecular barcode tag, and a exon tag.  The exon and gene tags may not be present on every read." +
			"When filtering the data for a set of barcodes to use, the data is filtered by ONE of the following methods (and if multiple params are filled in, the top one takes precidence):\n" +
			"1) CELL_BC_FILE, to filter by the some fixed list of cell barcodes" + 
			"2) MIN_NUM_GENES_PER_CELL " + 
			"3) MIN_NUM_TRANSCRIPTS_PER_CELL " +
			"4) NUM_CORE_BARCODES " +
			"5) MIN_NUM_READS_PER_CELL",
        usageShort = "Get the number of reads for each UMI",
        programGroup = DropSeq.class
)
public class GatherMolecularBarcodeDistributionByGene extends DGECommandLineBase {
	
	private static final Log log = Log.getInstance(GatherMolecularBarcodeDistributionByGene.class);
	
	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file of with 4 columns: CELL, GENE, MOLECULAR BC, #Observations. This supports zipped formats like gz and bz2.")
	public File OUTPUT;
	
	@Override
	protected int doWork() {
				
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT);
		
		writePerTranscriptHeader(out);
		
		List<String> barcodes=new BarcodeListRetrieval().getCellBarcodes(this.INPUT, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, 
				this.GENE_EXON_TAG, this.STRAND_TAG, this.CELL_BC_FILE, this.READ_MQ, this.MIN_NUM_TRANSCRIPTS_PER_CELL, 
				this.MIN_NUM_GENES_PER_CELL, this.MIN_NUM_READS_PER_CELL, this.NUM_CORE_BARCODES, this.EDIT_DISTANCE, this.MIN_BC_READ_THRESHOLD, USE_STRAND_INFO);
				
		UMIIterator umiIterator = new UMIIterator(this.INPUT, this.GENE_EXON_TAG, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.STRAND_TAG, 
				this.READ_MQ, true, this.USE_STRAND_INFO, barcodes, this.MAX_RECORDS_IN_RAM);
				
		
		UMICollection batch;
		
		while ((batch=umiIterator.next())!=null) {
			if (batch.isEmpty()==false) {
				String cellTag = batch.getCellBarcode();
				if (barcodes.contains(cellTag) || barcodes.isEmpty()) {
					writePerTranscriptStats (batch.getGeneName(), batch.getCellBarcode(), batch.getMolecularBarcodeCountsCollapsed(this.EDIT_DISTANCE), out);
				}
			}
		}
		
		
		CloserUtil.close(out);
		return 0;
	}
	
		
	private void writePerTranscriptStats (String gene, String cellBarcode, ObjectCounter<String> counts, BufferedWriter out) {
		for (String key: counts.getKeys()) {
			int value = counts.getCountForKey(key);
			String [] line ={cellBarcode, gene, key, value+""};
			String h = StringUtils.join(line, "\t");
			OutputWriterUtil.writeResult(h, out);
		}
	}
	
	
	private void writePerTranscriptHeader(BufferedWriter out) {
		String [] header = {"Cell Barcode", "Gene", "Molecular_Barcode", "Num_Obs"};
		String h = StringUtils.join(header, "\t");
		OutputWriterUtil.writeResult(h, out);
	}
	
	
	
	/**
	 * Get the number of transcripts for each cell barcode.
	 * @param bamFile The input BAM file
	 * @param cellTag the tag for the cell barcode
	 * @param molecularBarcodeTag the tag for the molecular barcode
	 * @param geneExonTag
	 * @param mapQuality The minimum map quality for each read to be considered.
	 * @param useStrandInfo
	 * @return
	 */
	public ObjectCounter<String> getNumTranscriptsPerCell (File bamFile, String cellTag, String molecularBarcodeTag, 
			String geneExonTag, String strandTag, Integer mapQuality, int editDistance, int minNumReadsMolBarcode, boolean useStrandInfo) {
		
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT);
		
		writePerTranscriptHeader(out);
		
		UMIIterator umiIterator = new UMIIterator(bamFile, geneExonTag, cellTag, molecularBarcodeTag, strandTag, 
				mapQuality, true, useStrandInfo, null, this.MAX_RECORDS_IN_RAM);
		
		
		ObjectCounter<String> transcriptsPerCell = new ObjectCounter<String>();
		
		UMICollection batch;
		while ((batch=umiIterator.next())!=null) {
			if (batch.isEmpty()==false) {
				int numTranscripts = batch.getMolecularBarcodeCounts().getSize();
				transcriptsPerCell.incrementByCount(batch.getCellBarcode(), numTranscripts);
			}
		}
		return (transcriptsPerCell);
		
	}
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new GatherMolecularBarcodeDistributionByGene().instanceMain(args));
	}

}
