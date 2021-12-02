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
package org.broadinstitute.dropseqrna.barnyard;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.OutputWriterUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.annotation.LocusFunction;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        summary = "For each gene, count the number of times each molecular barcode is observed [The UMI]" +
			"Similar to digital expression, reads are filtered on map quality, and must overlap exons as well as genes. "+
			"This program requires a tag for what gene a read is on, a molecular barcode tag, and a exon tag.  The exon and gene tags may not be present on every read." +
			"When filtering the data for a set of barcodes to use, the data is filtered by ONE of the following methods (and if multiple params are filled in, the top one takes precedence):\n" +
			"1) CELL_BC_FILE, to filter by the some fixed list of cell barcodes" +
			"2) MIN_NUM_GENES_PER_CELL " +
			"3) MIN_NUM_TRANSCRIPTS_PER_CELL " +
			"4) NUM_CORE_BARCODES " +
			"5) MIN_NUM_READS_PER_CELL",
        oneLineSummary = "Get the number of reads for each UMI",
        programGroup = DropSeq.class
)
public class GatherMolecularBarcodeDistributionByGene extends DGECommandLineBase {

	private static final Log log = Log.getInstance(GatherMolecularBarcodeDistributionByGene.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<File> INPUT;
	
	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file of with 4 columns: CELL, GENE, MOLECULAR BC, #Observations. This supports zipped formats like gz and bz2.")
	public File OUTPUT;

	@Override
	protected int doWork() {

		INPUT = FileListParsingUtils.expandFileList(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT);

		writePerTranscriptHeader(out);

		Set<String> cellBarcodes=new HashSet<>(new BarcodeListRetrieval().getCellBarcodes(this.INPUT, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG,
                this.GENE_NAME_TAG, this.GENE_STRAND_TAG, this.GENE_FUNCTION_TAG, this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST,
                this.CELL_BC_FILE, this.READ_MQ, this.MIN_NUM_TRANSCRIPTS_PER_CELL,
                this.MIN_NUM_GENES_PER_CELL, this.MIN_NUM_READS_PER_CELL, this.NUM_CORE_BARCODES, this.EDIT_DISTANCE, this.MIN_BC_READ_THRESHOLD));

		UMIIterator umiIterator = new UMIIterator(SamFileMergeUtil.mergeInputs(this.INPUT, false),
				GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
        		this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG,
        		this.READ_MQ, false, cellBarcodes, true, false);

		UMICollection batch;

		while ((batch=umiIterator.next())!=null)
			if (!batch.isEmpty()) {
				String cellTag = batch.getCellBarcode();
				if (cellBarcodes.contains(cellTag) || cellBarcodes.isEmpty())
					writePerTranscriptStats (batch.getGeneName(), batch.getCellBarcode(), batch.getMolecularBarcodeCountsCollapsed(this.EDIT_DISTANCE), out);
			}

		CloserUtil.close(umiIterator);

		try {
			out.close();
		} catch (IOException io) {
			throw new TranscriptomeException("Problem writing file", io);
		}

		return 0;
	}


	private void writePerTranscriptStats (final String gene, final String cellBarcode, final ObjectCounter<String> counts, final BufferedWriter out) {
		for (String key: counts.getKeys()) {
			int value = counts.getCountForKey(key);
			String [] line ={cellBarcode, gene, key, value+""};
			String h = StringUtils.join(line, "\t");
			OutputWriterUtil.writeResult(h, out);
		}
	}


	public static void writePerTranscriptHeader(final BufferedWriter out) {
		String [] header = {"Cell Barcode", "Gene", "Molecular_Barcode", "Num_Obs"};
		String h = StringUtils.join(header, "\t");
		OutputWriterUtil.writeResult(h, out);
	}



	/**
	 * I think this might still be horribly inefficient for some tasks - it calculates the number of UMIs for every cell, doing collapse along the way.
	 * This works on even cells that have a single molecular barcode, meaning that if you're using this to find cells with > 1000 UMis (for example), it could be pretty slow.
	 * Get the number of transcripts for each cell barcode.
	 * @param bamFile The input BAM file
	 * @param cellBarcodeTag the tag for the cell barcode
	 * @param molecularBarcodeTag the tag for the molecular barcode
	 * @param mapQuality The minimum map quality for each read to be considered.
	 * @return
	 */
	public ObjectCounter<String> getNumTranscriptsPerCell (final List<File> bamFile, final String cellBarcodeTag, final String molecularBarcodeTag,
			final String geneNameTag, final String strandTag, final String geneFunctionTag, final StrandStrategy strategy, final Collection<LocusFunction> locusFunctionList,
			final Integer mapQuality, final int editDistance, final int minNumReadsMolBarcode, final Collection<String> cellBarcodes) {

		SamReaderFactory factory= SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE);		
		SamHeaderAndIterator headerIterator= SamFileMergeUtil.mergeInputs(bamFile, false, factory);
		
		UMIIterator umiIterator = new UMIIterator(headerIterator,
				geneNameTag, strandTag, geneFunctionTag,
				strategy, locusFunctionList, cellBarcodeTag, molecularBarcodeTag,
				mapQuality, false, cellBarcodes);

		ObjectCounter<String> transcriptsPerCell = new ObjectCounter<>();

		UMICollection batch;
		while ((batch=umiIterator.next())!=null)
			if (!batch.isEmpty()) {
				int numTranscripts = batch.getMolecularBarcodeCountsCollapsed(editDistance).getSize();
				transcriptsPerCell.incrementByCount(batch.getCellBarcode(), numTranscripts);
			}
		umiIterator.close();
		return (transcriptsPerCell);

	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new GatherMolecularBarcodeDistributionByGene().instanceMain(args));
	}

}
