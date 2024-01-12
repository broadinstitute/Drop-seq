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
import htsjdk.samtools.util.*;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorStrategyEnum;
import org.broadinstitute.dropseqrna.metrics.BamTagHistogram;
import org.broadinstitute.dropseqrna.metrics.BamTagOfTagCounts;
import org.broadinstitute.dropseqrna.metrics.TagOfTagResults;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import picard.annotation.LocusFunction;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class BarcodeListRetrieval {

	private static Log log = Log.getInstance(BarcodeListRetrieval.class);
	private ProgressLogger progress = new ProgressLogger(log, 1000000);


	/**
	 * Convenience method to select a set of cell barcodes to work with dependent on which parameters are filled in.
	 * @param bamFile
	 * @param cellBarcodeTag
	 * @param molecularBarcodeTag
	 * @param cellBCFile
	 * @param readQuality
	 * @param minNumTranscriptsPerCell
	 * @param minNumNumGenesPerCell
	 * @param minNumReadsPerCell
	 * @param numCoreBarcodes
	 * @return
	 */
	public List<String> getCellBarcodes(final List<File> bamFile, final String cellBarcodeTag, final String molecularBarcodeTag,
			final String geneNameTag, final String strandTag, final String geneFunctionTag, final StrandStrategy strategy, final Collection<LocusFunction> locusFunctionList, final FunctionalDataProcessorStrategyEnum functionStrategy,
										final File cellBCFile, final Integer readQuality, final Integer minNumTranscriptsPerCell,
			final Integer minNumNumGenesPerCell, final Integer minNumReadsPerCell, final Integer numCoreBarcodes, final Integer editDistance, final Integer minNumReadsMolBarcode) {
		List<String> cellBarcodes=new ArrayList<String>();

		if (cellBCFile!=null) {
			cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
			log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
			return (cellBarcodes);
		}

		if (minNumNumGenesPerCell!=null){
			log.info("Looking for cell barcodes that have at least " + minNumNumGenesPerCell + " genes");
			cellBarcodes = getListCellBarcodesByGeneCount(bamFile, cellBarcodeTag, geneNameTag, readQuality, minNumNumGenesPerCell);
			return (cellBarcodes);
		}

		if (minNumTranscriptsPerCell!=null) {
			log.info("Looking for cell barcodes that have at least " + minNumTranscriptsPerCell + " transcripts");

			/*
			 * (final File input, final String cellBarcodeTag, final String molecularBarcodeTag, final String geneNameTag,
			final String strandTag, final String geneFunctionTag, final StrandStrategy strategy, final Collection<LocusFunction> locusFunctionList,
			final int readQuality, final int editDistance, final int minNumReadsMolBarcode, final int minNumTranscripts) {
			 */
			cellBarcodes = getListCellBarcodesByTranscriptCount(bamFile, cellBarcodeTag, molecularBarcodeTag, geneNameTag, strandTag, geneFunctionTag,
					 strategy, locusFunctionList, functionStrategy, readQuality, editDistance, minNumReadsMolBarcode, minNumTranscriptsPerCell);
			return (cellBarcodes);
		}

		if (minNumReadsPerCell!=null || numCoreBarcodes!=null){
			cellBarcodes = getListCellBarcodesByReadCount (bamFile, cellBarcodeTag, readQuality, minNumReadsPerCell, numCoreBarcodes);
			return (cellBarcodes);
		}

		return (cellBarcodes);
	}

    /**
     * Returns a list of cell barcodes with more at least MIN_NUM_READS_PER_CELL reads.
     * Or, if the numReadsCore is set, use that instead to get a list of cell barcodes.
     * If there are no cell barcodes, return an empty set.
     * @return
     */
    public List<String> getListCellBarcodesByReadCount (final List<File> input, final String cellBarcodeTag, final int readQuality, final Integer minNumReads, final Integer numReadsCore) {
    	SamHeaderAndIterator headerAndIter= SamFileMergeUtil.mergeInputs(input, false, SamReaderFactory.makeDefault());
        return getListCellBarcodesByReadCount(
        		headerAndIter.iterator,
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

        BamTagHistogram bth = new BamTagHistogram();
        ObjectCounter<String> cellBarcodes = bth.getBamTagCounts (input, cellBarcodeTag, readQuality, false);

        List<String> result=null;

        if (minNumReads!=null)
			result = getCoreBarcodesByReadCount(cellBarcodes, minNumReads);
        if (numReadsCore!=null)
			result = getTopCoreBarcodesByReadCount (cellBarcodes, numReadsCore);
        CloserUtil.close(input);
        return (result);
    }

    public List<String> getListCellBarcodesByGeneCount(final List<File> input, final String cellBarcodeTag, final String geneExonTag, final int readQuality, final int minNumGenes) {
		List<String> result = new ArrayList<String>();

		BamTagOfTagCounts tot = new BamTagOfTagCounts();
		TagOfTagResults<String,String> results = tot.getResults(input, cellBarcodeTag, geneExonTag, false, readQuality);
		for (String key: results.getKeys()) {
			int count = results.getCount(key);
			if (count>=minNumGenes)
				result.add(key);
		}
		return (result);
	}

	public List<String> getListCellBarcodesByTranscriptCount(final List<File> input, final String cellBarcodeTag, final String molecularBarcodeTag, final String geneNameTag,
			final String strandTag, final String geneFunctionTag, final StrandStrategy strategy, final Collection<LocusFunction> locusFunctionList, final FunctionalDataProcessorStrategyEnum functionStrategy,
			final int readQuality, final int editDistance, final int minNumReadsMolBarcode, final int minNumTranscripts) {
		
		List<String> result = new ArrayList<String>();
		List<String> initialCells = getListCellBarcodesByReadCount(input, cellBarcodeTag, readQuality, minNumTranscripts, null);

		GatherMolecularBarcodeDistributionByGene g = new GatherMolecularBarcodeDistributionByGene();
		/**
		 * (final File bamFile, final String cellBarcodeTag, final String molecularBarcodeTag,
			final String geneNameTag, final String strandTag, final String geneFunctionTag, final StrandStrategy strategy, final Collection<LocusFunction> locusFunctionList,
			final Integer mapQuality, final int editDistance, final int minNumReadsMolBarcode, final boolean useStrandInfo, final List<String> cellBarcodes)
		 */
		ObjectCounter<String> numTranscripts = g.getNumTranscriptsPerCell (input,cellBarcodeTag, molecularBarcodeTag,
				geneNameTag, strandTag, geneFunctionTag, strategy, locusFunctionList, functionStrategy,
				readQuality, editDistance, minNumReadsMolBarcode, initialCells);

		for (String key: numTranscripts.getKeys()) {
			int count = numTranscripts.getCountForKey(key);
			if (count>=minNumTranscripts)
				result.add(key);
		}
		return (result);
	}

	/*
	public List<String> getListCellBarcodesByGeneCount(File input, String cellBarcodeTag, int readQuality, int minNumGenes) {
		BamTagOfTagCounts bt = new BamTagOfTagCounts();
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
