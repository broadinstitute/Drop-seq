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
package org.broadinstitute.dropseqrna.beadsynthesis;

import java.io.File;
import java.util.List;

import org.broadinstitute.dropseqrna.barnyard.BarcodeListRetrieval;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

/**
 * Generate a collection of BeadSynthesisErrorData objects from a BAM file.
 * @author nemesh
 *
 */
public class BiasedBarcodeCollectionFactory {

	private static final Log log = Log.getInstance(BiasedBarcodeCollectionFactory.class);

	private SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE);

	/**
	 * Set up the UMI Iterator.
	 * @return
	 */
	public UMIIterator prepareUMIIterator(final List<File> inputFiles, final String geneExonTag, final String cellBarcodeTag, final String molBCTag, final String strandTag,
			final int readMQ, final List<String> cellBarcodes) {

		return new UMIIterator(SamFileMergeUtil.mergeInputs(inputFiles, false, samReaderFactory),
				geneExonTag, cellBarcodeTag, molBCTag, strandTag, readMQ, false, false, cellBarcodes, true);

	}

	/**
	 * Get a list of cell barcodes that can be queried in the BAM.
	 * @param inputFiles One or more BAM files
	 * @param cellBarcodeFile A file containing a list of cell barcodes, can be null.
	 * @param numBarcodes A count of the number of barcodes to extract.  Use if cellBarcodeFile is null.
	 * @param cellBarcodeTag The cell barcode tag.  Use if cellBarcodeFile is null.
	 * @param readMQ The map quality of reads to use if cellBarcodeFile is null.
	 * @return
	 */
	public List<String> getCellBarcodes (final List<File> inputFiles, final File cellBarcodeFile, final int numBarcodes, final String cellBarcodeTag, final int readMQ) {

		if (cellBarcodeFile!=null) {
			IOUtil.assertFileIsReadable(cellBarcodeFile);
			List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBarcodeFile);
			log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
			return (cellBarcodes);
		}
		log.info("Gathering barcodes for the top [" + numBarcodes +"] cells");
        return new BarcodeListRetrieval().getListCellBarcodesByReadCount(
                SamFileMergeUtil.mergeInputs(inputFiles, false, samReaderFactory).iterator,
                cellBarcodeTag, readMQ, null, numBarcodes);
	}




}
