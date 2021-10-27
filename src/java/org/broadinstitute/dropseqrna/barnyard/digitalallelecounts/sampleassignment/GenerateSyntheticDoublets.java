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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.BarcodeListRetrieval;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(summary = "Takes a fraction of the cell barcodes and merges the barcode labels and reads of those cells into synthetic doublets.  For example if cell A and B were"
		+ "selected to be a doublet, then a new cell called A_B would be generated that contained all the reads from the A and B cell barcodes.", oneLineSummary = "", programGroup = DropSeq.class)
public class GenerateSyntheticDoublets extends CommandLineProgram {
	private static final Log log = Log
			.getInstance(GenerateSyntheticDoublets.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<File> INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM file with generated multi-cells.")
	public File OUTPUT;

	@Argument (doc="Should singleton cells be emitted?  If set to false, only the multicells that are generated are emitted.")
	public Boolean EMIT_SINGLETONS=false;

	@Argument(doc="Output file of the newly generated multicells.  This file has 1 column with the names of the cells that were combined, colon seperated. "
			+ "This is the same format as MULTI_CELL_FILE, and can be used as input in a different run. ", optional=true)
	public File SUMMARY;

	@Argument(doc="How many multicells should be generated from the data?  If the NUMBER_MULTICELL*NUM_CELLS_PER_MULTICELL > "
			+ "number of total cells, the program will report this as an error and exit.", mutex={"MULTI_CELL_FILE"})
	public Integer NUMBER_MULTICELL;

	@Argument(doc="How many cells should each multi-cell be made from?", mutex={"MULTI_CELL_FILE"})
	public int MULTIPLICITY=2;

	@Argument(doc="Instead of randomly selecting a fraction of cells as multicells, input a file defining the cell barcodes that will generate a doublet. "
			+ "The file has 1 column, for each cell barcode.  Each row defines a single multicell. The entry is a colon seperated list of cell barcodes that represent all the cells that should be merged. "
			+ "For example, if you wanted to merge cells A B and C, the entry would be A:B;C. The file is tab delimited with no column header.", mutex={"NUMBER_MULTICELL"})
	public File MULTI_CELL_FILE;

	@Argument(doc = "The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
	public String CELL_BARCODE_TAG = "XC";

	@Argument(doc = "Number of cells that you want to extract from the BAM. The program will picks the top <NUM_BARCODES> barcodes by read count.", mutex = { "CELL_BC_FILE" }, optional=true)
	public Integer NUM_BARCODES = null;

	@Argument(doc = "Override NUM_CORE_BARCODES and process reads that have the cell barcodes in this file instead.  The file has 1 column with no header.", mutex = { "NUM_BARCODES" }, optional=true)
	public File CELL_BC_FILE = null;

	@Argument(doc = "The map quality of the read to be included.")
	public Integer READ_MQ = 10;

	@Override
	protected int doWork() {
		this.INPUT = FileListParsingUtils.expandFileList(INPUT);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		if (this.SUMMARY!=null)
			IOUtil.assertFileIsWritable(this.SUMMARY);

		List<MultiCell> multiCells = getMultiCells(this.MULTI_CELL_FILE, this.NUMBER_MULTICELL, this.MULTIPLICITY);
		Map<String, List<String>> multiCellMap = MultiCell.getMap(multiCells);

		// now just iterate through reads and write out the BAM with modified cell barcode tags.
		SamHeaderAndIterator headerAndIter= SamFileMergeUtil.mergeInputs(this.INPUT, false, SamReaderFactory.makeDefault());
		
		SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(headerAndIter.header, true, OUTPUT);

		ProgressLogger pl = new ProgressLogger(this.log);

		for (Iterator<SAMRecord> i = headerAndIter.iterator; i.hasNext(); ) {
			SAMRecord r = i.next();		
			pl.record(r);
			String cellBarcode = r.getStringAttribute(this.CELL_BARCODE_TAG);
			List<String> newBarcodes = multiCellMap.get(cellBarcode);

			// cell is in MultiCell, change the barcode and write.
			if (newBarcodes!=null) {
				// write out the original read before writing out retagged reads.
				if (EMIT_SINGLETONS)
					writer.addAlignment(r);
				// for each multicell this cell barcode partipates in, write out the read with the new tag.
				for (String newBarcode: newBarcodes) {
					r.setAttribute(this.CELL_BARCODE_TAG, newBarcode);
					writer.addAlignment(r);
				}
			}
		}		
		writer.close();

		// write summary
		if (this.SUMMARY!=null)
			MultiCell.writeToFile(multiCells, this.SUMMARY);

		return 0;
	}



	/**
	 * Generate a list of doublets
	 * This is either parsed in from a doubletFile, or if that parameter is null
	 * is generated from a list of cell barcodes and a fraction of doublets.
	 * @param doubletFile A file containing a list of cells to be merged into doublets
	 * @param cellBarcodes If doubletFile is null, a list of cell barcodes to merge into doublets
	 * @param  numMultiCells The number of multi-cells created.
	 * @param multiplicity how many cells go into the multi-cell
	 * @return A list of doublets.
	 */
	private List<MultiCell> getMultiCells (final File multiCellFile, final Integer numMultiCells, final Integer multiplicity) {
		// if there's a doublet file
		if (multiCellFile!=null) {
			IOUtil.assertFileIsReadable(multiCellFile);
			return MultiCell.parseFile(multiCellFile);
		}
		// do it the hard way and generate multiCells
		List<String> cellBarcodes = getCellBarcodes();
		List<MultiCell> result = MultiCell.generateMultiCells(cellBarcodes, numMultiCells, multiplicity);
		return result;
	}


	// done before, pretty well tested elsewhere boilerplate.
	private List<String> getCellBarcodes() {
		if (this.CELL_BC_FILE != null) {
			IOUtil.assertFileIsReadable(this.CELL_BC_FILE);
			List<String> cellBarcodes = ParseBarcodeFile
					.readCellBarcodeFile(this.CELL_BC_FILE);
			log.info("Found " + cellBarcodes.size() + " cell barcodes in file");
			return (cellBarcodes);
		}
		log.info("Gathering barcodes for the top [" + this.NUM_BARCODES
				+ "] cells");
		return new BarcodeListRetrieval().getListCellBarcodesByReadCount(
				this.INPUT, this.CELL_BARCODE_TAG, this.READ_MQ, null,
				this.NUM_BARCODES);
	}


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new GenerateSyntheticDoublets().instanceMain(args));
	}
}
