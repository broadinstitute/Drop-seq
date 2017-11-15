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
package org.broadinstitute.dropseqrna.utils.editdistance;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.broadinstitute.dropseqrna.barnyard.BarcodeListRetrieval;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.metrics.BAMTagHistogram;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

/**
 * Fold down barcodes, possibly in the context of another barcode (that has been folded down already.)
 * @author nemesh
 *
 */
@CommandLineProgramProperties(usage = "Fold down barcodes, possibly in the context of another barcode (that has been folded down already.)",
        usageShort = "Fold down barcodes, possibly in the context of another barcode (that has been folded down already.)",
        programGroup = DropSeq.class)
public class CollapseBarcodesInPlace extends CommandLineProgram {


	private final Log log = Log.getInstance(CollapseBarcodesInPlace.class);
	private ProgressLogger pl = new ProgressLogger(this.log);

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.  Must be coordinate sorted. ",
            minElements = 1)
	public List<File> INPUT;

	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BAM file with an extra tag.")
	public File OUTPUT;

	@Option(doc="Barcode to collapse")
	public String PRIMARY_BARCODE;

	@Option(doc="The edit distance to collapse barcodes")
	public Integer EDIT_DISTANCE=1;

	@Option(doc = "Should indels be considered in edit distance calculations?  Doing this correctly is far slower than a simple edit distance test, but gives a more complete result.")
	public boolean FIND_INDELS=true;

	@Option(doc="The output barcode tag for the newly collapsed barcodes")
	public String OUT_BARCODE;

	@Option(doc="Read quality filter.  Filters all reads lower than this mapping quality.  Defaults to 10.  Set to 0 to not filter reads by map quality.")
	public Integer READ_QUALITY=10;

	@Option(doc="Number of reads a barcode would need to have in order to have other barcodes get merged into it.  All barcodes are candidates to be merged into another barcode." +
			"For cell barcodes you probably want to set this to a relatively high number like 100, since we expect cells to have thousands or more reads, and this signficantly speeds up analysis.  " +
			"For molecular barcodes, you probably want to set this to 1, as you want to include all molecular barcodes, unless you have very high sequencing depth.", optional=true)
	public Integer MIN_NUM_READS_CORE=null;

	@Option(doc="Number of cells that you think are in the library.  This accomplishes the same goals as the MIN_NUM_READS_CORE argument, but instead of defining barcodes as important based on the number of reads, it picks the top <X> barcodes as core.", optional=true)
	public Integer NUM_CORE_BARCODES=null;

	@Option(doc="The number of reads a non-core barcode must have to be merged with a core barcode.", optional=true)
	public Integer MIN_NUM_READS_NONCORE=1;

	@Option(doc="Filter PCR Duplicates.  Defaults to false")
	public boolean FILTER_PCR_DUPLICATES=false;

	@Option(doc="Number of threads to use.  Defaults to 1.")
	public int NUM_THREADS=1;


	private int REPORT_PROGRESS_INTERVAL=100;
	private CollapseBarcodeThreaded cbt=null;
	private int threadedBlockSize=20000;

	@Override
	protected int doWork() {
		log.info("Number of cores selected [" + Integer.toString(this.NUM_THREADS) + "]");

		if (this.NUM_THREADS>1) cbt= new CollapseBarcodeThreaded(this.threadedBlockSize, this.NUM_THREADS);

		IOUtil.assertFileIsWritable(OUTPUT);
		for (final File inputFile: INPUT)
			IOUtil.assertFileIsReadable(inputFile);

        processOnlyPrimary();

		return 0;
	}


	public void processOnlyPrimary () {
        final SamHeaderAndIterator inputs = openInputs();
		CloseableIterator<SAMRecord> inputSam = inputs.iterator;
		SAMFileHeader header = inputs.header;
		header.addComment("Edit distance collapsed tag " +  this.PRIMARY_BARCODE + " to new tag " + this.OUT_BARCODE+ " with edit distance "+ this.EDIT_DISTANCE);
        SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, this.OUTPUT);

		// gather up the barcodes that exist in the BAM
        final SamHeaderAndIterator inputs2 = openInputs();
		ObjectCounter<String> barcodes = new BAMTagHistogram().getBamTagCounts(inputs2.iterator, this.PRIMARY_BARCODE,this.READ_QUALITY, this.FILTER_PCR_DUPLICATES);
        CloserUtil.close(inputs2.iterator);

		// filter barcodes by #reds in each barcode.
		barcodes=filterBarcodesByNumReads(barcodes, this.MIN_NUM_READS_NONCORE);

		// collapse them
		Map<String, String> childParentBarcodes=collapseBarcodes(this.MIN_NUM_READS_CORE, this.NUM_CORE_BARCODES, barcodes, this.FIND_INDELS, this.EDIT_DISTANCE);
		// iterate through the reads and retag with the proper reads.
		// log.info("STUFF");
		retagReads(inputSam, writer, childParentBarcodes, this.PRIMARY_BARCODE, this.OUT_BARCODE);
		// collapsed.size();

		CloserUtil.close(inputSam);
		writer.close();
	}

    private SamHeaderAndIterator openInputs() {
        final SamHeaderAndIterator ret = SamFileMergeUtil.mergeInputs(INPUT, true);
        if (SAMFileHeader.SortOrder.coordinate != ret.header.getSortOrder())
			throw new PicardException("Input files are not coordinate sorted");
        return ret;
    }

	private ObjectCounter<String> filterBarcodesByNumReads (final ObjectCounter<String> barcodes, final int minNumReads) {

		ObjectCounter<String> result = new ObjectCounter<>();
		for (String k: barcodes.getKeys()) {
			int count = barcodes.getCountForKey(k);
			if (count>=minNumReads)
				result.setCount(k, count);
		}
		log.info("Filtering barcodes by min non-core reads.  Started with [" + barcodes.getSize()+ "] ended with ["+ result.getSize()+"]");
		return (result);
	}

	private void retagReads (final Iterator<SAMRecord> inputSam, final SAMFileWriter writer, final Map<String, String> collapsed, final String tag, final String outTag) {
		for (final SAMRecord r : new IterableAdapter<>(inputSam)) {
			pl.record(r);
			String s1 = r.getStringAttribute(tag);
			String newTag=collapsed.get(s1);
			if (newTag==null) newTag=s1;
			r.setAttribute(outTag, newTag);
			writer.addAlignment(r);
		}
	}

	/**
	 * Collapse barcodes.  Does this for all barcodes.
	 * @param barcodes
	 * @return A map of each child barcode to it's parent.  Many keys will point to the same value.
	 */

	private Map<String, String> collapseBarcodes(final ObjectCounter<String> barcodes, final boolean findIndels, final int editDistance) {
		List<String> barcodeList = barcodes.getKeysOrderedByCount(true);
		Map<String, String> result = collapseBarcodes(barcodeList, barcodes, findIndels, editDistance);
		return (result);
	}


	/**
	 * Convenience method
	 * @param numReadsCore
	 * @param barcodes
	 * @param findIndels
	 * @param editDistance
	 * @return
	 */
	private Map<String, String> collapseBarcodes(final Integer numReadsCore, final Integer numCells, final ObjectCounter<String> barcodes, final boolean findIndels, final int editDistance) {
		if (numReadsCore==null && numCells==null) return (collapseBarcodes(barcodes, findIndels, editDistance));
		// otherwise, select core barcodes and run.
		List<String> core=null;
		BarcodeListRetrieval u = new BarcodeListRetrieval();

		if (numReadsCore!=null)
			core = u.getCoreBarcodesByReadCount(barcodes, numReadsCore);
		else if (numCells!=null)
			core = u.getTopCoreBarcodesByReadCount (barcodes, numCells);

		return (collapseBarcodes(core, barcodes, findIndels, editDistance));
	}


	public Map<String, String> collapseBarcodes(final List<String> coreBarcodes, final ObjectCounter<String> barcodes, final boolean findIndels, final int editDistance) {
		Map<String, String> result = new HashMap<>();

		MapBarcodesByEditDistance med = new MapBarcodesByEditDistance(true, this.NUM_THREADS, 10000);
		Map<String, List<String>> r = med.collapseBarcodes(coreBarcodes, barcodes, findIndels, editDistance);
		for (String key: r.keySet())
			for (String value: r.get(key))
				result.put(value, key);
		return (result);

	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CollapseBarcodesInPlace().instanceMain(args));
	}

}
