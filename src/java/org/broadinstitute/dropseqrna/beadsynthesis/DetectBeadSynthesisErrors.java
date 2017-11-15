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
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.BaseDistributionMetric;
import org.broadinstitute.dropseqrna.utils.BaseDistributionMetricCollection;
import org.broadinstitute.dropseqrna.utils.Bases;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;

/**
 *
 * @author nemesh
 *
 */

@CommandLineProgramProperties(
        usage = "For each cell, gather up all the UMIs.  An error in synthesis will result in the last base of the synthesis being fixed in >90% of the UMIs for that cell, across all genes." +
			"This fixed base is T.  For cell barcodes where this occurs, output the cell barcode in a file, as well as (optionally) pad the cell barcodes with N for the error bases.",
        usageShort = "Detect barcode synthesis errors where the final base of a UMI is fixed across all UMIs of a cell.",
        programGroup = DropSeq.class
)

public class DetectBeadSynthesisErrors extends AbstractDetectBeadSynthesisErrors {

	private static final Log log = Log.getInstance(DetectBeadSynthesisErrors.class);

	@Option(doc="Output of detailed information on each cell barcode analyzed.  Each row is a single cell barcode.  "
			+ "The data has multiple columns: the cell barcode, the number of UMIs, then one column per UMI base position containing the count of the reads, with a | "
			+ "delimiter between bases.  Bases are ordered A,C,G,T for these columns.  An example output with a single base UMI would be:"
			+ "AAAAAA	20		5|4|6|5.")
	public File OUTPUT_STATS;

	private Double EXTREME_BASE_RATIO=0.8;


	private DetectPrimerInUMI detectPrimerTool=null;

	@Override
	protected int doWork() {
		// primer detection if requested.
		if (this.PRIMER_SEQUENCE!=null)
			detectPrimerTool = new DetectPrimerInUMI(this.PRIMER_SEQUENCE);

		UMIIterator iterator = prepareUMIIterator();

		BiasedBarcodeCollection biasedBarcodeCollection = findBiasedBarcodes(iterator);
		Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions = biasedBarcodeCollection.getBiasedBarcodes();
		int numCellsFilteredLowUMIs = biasedBarcodeCollection.getNumBarcodesFilteredLowUMIs();

        // initialize output writers.
        PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT_STATS));
		writeFile (errorBarcodesWithPositions.values(), out);
		writeSummary(errorBarcodesWithPositions.values(), numCellsFilteredLowUMIs, SUMMARY);

		// clean up the BAM if desired.
		if (this.OUTPUT!=null)
			cleanBAM(errorBarcodesWithPositions);
		return 0;
	}



	/**
	 * For each problematic cell, replace cell barcodes positions with N.
	 * Take the replaced bases and prepend them to the UMI, and trim the last <X> bases off the end of the UMI.
	 */
	private void cleanBAM (final Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions) {
		log.info("Cleaning BAM");
        final SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(INPUT, true);
		SamHeaderUtil.addPgRecord(headerAndIterator.header, this);

		SAMFileWriter writer= new SAMFileWriterFactory().setCreateIndex(CREATE_INDEX).makeSAMOrBAMWriter(headerAndIterator.header, true, OUTPUT);
		ProgressLogger pl = new ProgressLogger(log);
		for (SAMRecord r: new IterableAdapter<>(headerAndIterator.iterator)) {
			pl.record(r);
			r=padCellBarcodeFix(r, errorBarcodesWithPositions, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.EXTREME_BASE_RATIO);
			if (r!=null)
				writer.addAlignment(r);
		}
		CloserUtil.close(headerAndIterator.iterator);
		writer.close();
	}


	/**
	 * @return null if the read should not be included in the output BAM.
	 */
	SAMRecord padCellBarcodeFix (final SAMRecord r, final Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions, final String cellBarcodeTag, final String molecularBarcodeTag, final double extremeBaseRatio) {
		String cellBC=r.getStringAttribute(cellBarcodeTag);

		BeadSynthesisErrorData bsed = errorBarcodesWithPositions.get(cellBC);
		if (bsed==null) return (r); // no correction data, no fix.

		// we're only going to fix cells where there's one or more synthesis errors
		BeadSynthesisErrorTypes bset = getEnhancedErrorType(bsed, extremeBaseRatio, this.detectPrimerTool);
		if (bset==BeadSynthesisErrorTypes.NO_ERROR) return (r); // no error, return.
		// has an error, not a synthesis error...
		if (bset!=BeadSynthesisErrorTypes.SYNTH_MISSING_BASE)
			return (null);

		// has a synthesis error
		int polyTErrorPosition = bsed.getPolyTErrorPosition(this.EXTREME_BASE_RATIO);
		int umiLength = bsed.getBaseLength();
		int numErrors= umiLength-polyTErrorPosition+1;
		// if there are too many errors, or the errors aren't all polyT, return null.
		if (numErrors > MAX_NUM_ERRORS)
			return null;

		// apply the fix and return the fixed read.
		String umi = r.getStringAttribute(molecularBarcodeTag);
		String cellBCFixed = padCellBarcode(cellBC, polyTErrorPosition, umiLength);
		String umiFixed = fixUMI(cellBC, umi, polyTErrorPosition);
		r.setAttribute(cellBarcodeTag, cellBCFixed);
		r.setAttribute(molecularBarcodeTag, umiFixed);
		return r;
	}





	private void writeSummary (final Collection <BeadSynthesisErrorData> data, final int numCellsFilteredLowUMIs, final File out) {
		// skip if no data.
		if (data.size()==0)
			return;

		// gather up error types.
		BeadSynthesisErrorsSummaryMetric m = new BeadSynthesisErrorsSummaryMetric();
		m.LOW_UMI_COUNT=numCellsFilteredLowUMIs;

		for (BeadSynthesisErrorData bsde: data) {
			BeadSynthesisErrorTypes t = getEnhancedErrorType(bsde, EXTREME_BASE_RATIO, this.detectPrimerTool);
			m.NUM_BEADS++;
			switch (t) {
				case SYNTH_MISSING_BASE: m.SYNTHESIS_MISSING_BASE++; m.incrementSynthesisMissingBase(bsde.getPolyTErrorPosition(this.EXTREME_BASE_RATIO)); break;
				case PRIMER: m.PRIMER_MATCH++; break;
				case SINGLE_UMI: m.SINGLE_UMI_ERROR++; break;
				case FIXED_FIRST_BASE: m.FIXED_FIRST_BASE++; break;
				case OTHER_ERROR: m.OTHER_ERROR_COUNT++; break;
				default: m.NO_ERROR++; break;
			}

		}

		MetricsFile<BeadSynthesisErrorsSummaryMetric, Integer> outFile = new MetricsFile<>();
		outFile.addMetric(m);
		outFile.addHistogram(m.getHistogram());
		outFile.write(out);

	}

	private void writeFile (final Collection <BeadSynthesisErrorData> data, final PrintStream out) {
		if (data.size()==0) {
			out.close();
			return;
		}

		// if there are records, write out the file.
        final List<BeadSynthesisErrorData> dataArray = new ArrayList<>(data);
        Collections.sort(dataArray, new Comparator<BeadSynthesisErrorData>() {
            @Override
            public int compare(final BeadSynthesisErrorData o1, final BeadSynthesisErrorData o2) {
                // Note that this is backwards so record with largest UMICount comes first
                int cmp = Integer.compare(o2.getUMICount(), o1.getUMICount());
                if (cmp != 0)
					return cmp;
                return o1.getCellBarcode().compareTo(o2.getCellBarcode());
            }
        });

		BeadSynthesisErrorData first = dataArray.get(0);
		int umiLength = first.getBaseLength();
		writeBadBarcodeStatisticsFileHeader(umiLength, out);
		for (BeadSynthesisErrorData bsde: dataArray)
			writeBadBarcodeStatisticsFileEntry(bsde, out);
		out.close();
	}

	/**
	 * Write the header.
	 */
	private void writeBadBarcodeStatisticsFileHeader (final int umiLength, final PrintStream out) {
		List<String> header = new ArrayList<>();
		header.add("CELL_BARCODE");
		header.add("NUM_UMI");
		header.add("FIRST_BIASED_BASE");
		header.add(BeadSynthesisErrorTypes.SYNTH_MISSING_BASE.toString());
		header.add("ERROR_TYPE");
		for (int i=0; i<umiLength; i++)
			header.add("BASE_"+ Integer.toString(i+1));
		String h = StringUtils.join(header, "\t");
		out.println(h);
	}

	private void writeBadBarcodeStatisticsFileEntry (final BeadSynthesisErrorData data, final PrintStream out) {
		List<String> line = new ArrayList<>();
		line.add(data.getCellBarcode());
		line.add(Integer.toString(data.getUMICount()));
		int base = data.getErrorBase(EXTREME_BASE_RATIO);
		line.add(Integer.toString(base));
		int polyTErrorBase = data.getPolyTErrorPosition(EXTREME_BASE_RATIO);
		line.add(Integer.toString(polyTErrorBase));
		line.add(getEnhancedErrorType(data, EXTREME_BASE_RATIO, this.detectPrimerTool).toString());


		BaseDistributionMetricCollection bases = data.getBaseCounts();
		List<Integer> pos =bases.getPositions();
		for (Integer i: pos) {
			BaseDistributionMetric bdm = bases.getDistributionAtPosition(i);
			String formattedResult = format(bdm);
			line.add(formattedResult);
		}
		String outLine = StringUtils.join(line, "\t");
		out.println(outLine);
	}



	private String format (final BaseDistributionMetric bdm) {
		List<String> d = new ArrayList<>();

		for (Bases b: Bases.values()) {
			char bb = b.getBase();
			int count = bdm.getCount(bb);
			d.add(Integer.toString(count));
		}
		return StringUtils.join(d, "|");
	}

    @Override
    protected String[] customCommandLineValidation() {
        for (final File input : INPUT)
			IOUtil.assertFileIsReadable(input);
        IOUtil.assertFileIsWritable(this.OUTPUT_STATS);
        IOUtil.assertFileIsWritable(this.SUMMARY);

        if (OUTPUT!=null) IOUtil.assertFileIsWritable(this.OUTPUT);
        return super.customCommandLineValidation();
    }



	public class BeadSynthesisErrorsSummaryMetric extends MetricBase {
		public int NUM_BEADS;
		public int NO_ERROR;
		public int LOW_UMI_COUNT;
		public int SYNTHESIS_MISSING_BASE;
		public int SINGLE_UMI_ERROR;
		public int PRIMER_MATCH;
		public int FIXED_FIRST_BASE;
		public int OTHER_ERROR_COUNT;


		/** The distribution of  SYNTHESIS_MISSING_BASE error positions */
		private Histogram <Integer> histogram = null;

		public BeadSynthesisErrorsSummaryMetric () {
			this.NUM_BEADS=0;
			this.NO_ERROR=0;
			this.SYNTHESIS_MISSING_BASE=0;
			this.SINGLE_UMI_ERROR=0;
			this.PRIMER_MATCH=0;
			this.OTHER_ERROR_COUNT=0;
			this.LOW_UMI_COUNT=0;
			histogram = new Histogram<>("SYNTHESIS_ERROR_BASE", "num cells");
		}

		public void incrementSynthesisMissingBase (final int position) {
			histogram.increment(position);
		}

		public Histogram<Integer> getHistogram() {
			return histogram;
		}

	}


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new DetectBeadSynthesisErrors().instanceMain(args));
	}
}
