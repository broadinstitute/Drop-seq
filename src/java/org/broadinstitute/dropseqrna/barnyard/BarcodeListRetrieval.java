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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.dropseqrna.metrics.BAMTagHistogram;
import org.broadinstitute.dropseqrna.metrics.BAMTagofTagCounts;
import org.broadinstitute.dropseqrna.metrics.TagOfTagResults;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

public class BarcodeListRetrieval {

	private static Log log = Log.getInstance(BarcodeListRetrieval.class);
	private ProgressLogger progress = new ProgressLogger(log, 1000000);


	/**
	 * Convenience method to select a set of cell barcodes to work with dependent on which parameters are filled in.
	 * @param bamFile
	 * @param cellBarcodeTag
	 * @param molecularBarcodeTag
	 * @param geneBarcodeTag
	 * @param exonTag
	 * @param cellBCFile
	 * @param readQuality
	 * @param minNumTranscriptsPerCell
	 * @param minNumNumGenesPerCell
	 * @param minNumReadsPerCell
	 * @param numCoreBarcodes
	 * @return
	 */
	public List<String> getCellBarcodes(final File bamFile, final String cellBarcodeTag, final String molecularBarcodeTag,
			final String geneExonBarcodeTag, final String strandTag, final File cellBCFile, final Integer readQuality, final Integer minNumTranscriptsPerCell,
			final Integer minNumNumGenesPerCell, final Integer minNumReadsPerCell, final Integer numCoreBarcodes, final Integer editDistance, final Integer minNumReadsMolBarcode, final boolean useStrandInfo) {
		List<String> cellBarcodes=new ArrayList<String>();

		if (cellBCFile!=null) {
			cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
			log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
			return (cellBarcodes);
		}

		if (minNumNumGenesPerCell!=null){
			log.info("Looking for cell barcodes that have at least " + minNumNumGenesPerCell + " genes");
			cellBarcodes = getListCellBarcodesByGeneCount(bamFile, cellBarcodeTag, geneExonBarcodeTag, readQuality, minNumNumGenesPerCell);
			return (cellBarcodes);
		}

		if (minNumTranscriptsPerCell!=null) {
			log.info("Looking for cell barcodes that have at least " + minNumTranscriptsPerCell + " transcripts");
			cellBarcodes = getListCellBarcodesByTranscriptCount(bamFile, cellBarcodeTag, molecularBarcodeTag, geneExonBarcodeTag, strandTag, readQuality, minNumTranscriptsPerCell, editDistance, minNumReadsMolBarcode, useStrandInfo);
		}

		if (minNumReadsPerCell!=null || numCoreBarcodes!=null){
			cellBarcodes = getListCellBarcodesByReadCount (bamFile, cellBarcodeTag, readQuality, minNumReadsPerCell, numCoreBarcodes);
			return (cellBarcodes);
		}

		return (cellBarcodes);
	}

	public boolean validateGetCellBarcodeListParams(final File bamFile, final String cellBarcodeTag, final String molecularBarcodeTag,
			final String geneBarcodeTag, final String exonTag, final File cellBCFile, final Integer readQuality, final Integer minNumTranscriptsPerCell,
			final Integer minNumNumGenesPerCell, final Integer minNumReadsPerCell, final Integer numCoreBarcodes) {

		IOUtil.assertFileIsReadable(bamFile);
		if (cellBarcodeTag==null) {
			log.error("cell barcode argument must be set");
			return false;
		}
		if (molecularBarcodeTag==null) {
			log.error("moleclar barcode argument must be set");
			return false;
		}
		if (cellBCFile!=null)  {
			IOUtil.assertFileIsReadable(bamFile);
			return true;
		}

		// everything below relies on read quality
		if (readQuality==null) {
			log.error("Read Quality must be set");
			return (false);
		}

		if (minNumNumGenesPerCell!=null)
			if (geneBarcodeTag==null) {
				log.error("Gene barcode tag argument must be set to filter by #genes");
				return false;
			}
		if (minNumTranscriptsPerCell!=null) {
			if (geneBarcodeTag==null) {
				log.error("Gene barcode tag argument must be set to filter by #transcripts");
				return false;
			}
			if (exonTag==null) {
				log.error("Exon barcode tag argument must be set to filter by #transcripts");
				return false;
			}

		}
		// if the minNumReadsPerCell is not null or the numCoreBarcodes then we have all the params set correctly and we're good.
		return true;
	}


    /**
     * Returns a list of cell barcodes with more at least MIN_NUM_READS_PER_CELL reads.
     * Or, if the numReadsCore is set, use that instead to get a list of cell barcodes.
     * If there are no cell barcodes, return an empty set.
     * @return
     */
    public List<String> getListCellBarcodesByReadCount (final File input, final String cellBarcodeTag, final int readQuality, final Integer minNumReads, final Integer numReadsCore) {
        return getListCellBarcodesByReadCount(
                SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(input).iterator(),
                cellBarcodeTag,
                readQuality,
                minNumReads,
                numReadsCore
        );
    }

    /**
     * Returns a list of cell barcodes with more at least MIN_NUM_READS_PER_CELL reads.
     * Or, if the numReadsCore is set, use that instead to get a list of cell barcodes.
     * If there are no cell barcodes, return an empty set.
     * Closes iterator before returning.
     */
    public List<String> getListCellBarcodesByReadCount(final CloseableIterator<SAMRecord> input, final String cellBarcodeTag, final int readQuality, final Integer minNumReads, final Integer numReadsCore) {

        BAMTagHistogram bth = new BAMTagHistogram();
        ObjectCounter<String> cellBarcodes = bth.getBamTagCounts (input, cellBarcodeTag, readQuality, false);

        List<String> result=null;

        if (minNumReads!=null)
			result = getCoreBarcodesByReadCount(cellBarcodes, minNumReads);
        if (numReadsCore!=null)
			result = getTopCoreBarcodesByReadCount (cellBarcodes, numReadsCore);
        CloserUtil.close(input);
        return (result);
    }

    public List<String> getListCellBarcodesByGeneCount(final File input, final String cellBarcodeTag, final String geneExonTag, final int readQuality, final int minNumGenes) {
		List<String> result = new ArrayList<String>();

		BAMTagofTagCounts tot = new BAMTagofTagCounts();
		TagOfTagResults<String,String> results = tot.getResults(input, cellBarcodeTag, geneExonTag, false, readQuality);
		for (String key: results.getKeys()) {
			int count = results.getCount(key);
			if (count>=minNumGenes)
				result.add(key);
		}
		return (result);
	}

	public List<String> getListCellBarcodesByTranscriptCount(final File input, final String cellBarcodeTag, final String molecularBarcodeTag, final String geneExonTag, final String strandTag, final int readQuality, final int editDistance, final int minNumReadsMolBarcode, final int minNumTranscripts, final boolean useStrandInfo) {
		List<String> result = new ArrayList<String>();
		GatherMolecularBarcodeDistributionByGene g = new GatherMolecularBarcodeDistributionByGene();
		ObjectCounter<String> numTranscripts = g.getNumTranscriptsPerCell (input,cellBarcodeTag, molecularBarcodeTag, geneExonTag, strandTag, readQuality, editDistance, minNumReadsMolBarcode, useStrandInfo);

		for (String key: numTranscripts.getKeys()) {
			int count = numTranscripts.getCountForKey(key);
			if (count>=minNumTranscripts)
				result.add(key);
		}
		return (result);
	}

	/*
	public List<String> getListCellBarcodesByGeneCount(File input, String cellBarcodeTag, int readQuality, int minNumGenes) {
		BAMTagofTagCounts bt = new BAMTagofTagCounts();
		// bt.
	}
	*/
	public List<String> getCoreBarcodesByReadCount(final ObjectCounter<String> barcodes, final Integer numReadsCore) {
		List<String> result = new ArrayList<String>();
		log.info("Looking for cell barcodes that have at least " + numReadsCore + " reads");
		for (String bc: barcodes.getKeysOrderedByCount(true)) {
			int count=barcodes.getCountForKey(bc);
			if (count>=numReadsCore)
				result.add(bc);
			else
				// sorted in order, so you're done.
				break;
		}
		log.info("Selected " + result.size() + " core barcodes");
		return (result);
	}

	public List<String> getTopCoreBarcodesByReadCount(final ObjectCounter<String> barcodes, final Integer numCells) {
		List<String> result = new ArrayList<String>();
		log.info("Looking for the top " + numCells +" cell barcodes");
		int numCellsAdded=0;
		for (String bc: barcodes.getKeysOrderedByCount(true))
			if (numCellsAdded<numCells) {
				result.add(bc);
				numCellsAdded++;
			} else
				break;
		log.info("Selected " + result.size() + " core barcodes");
		return (result);
	}


}
