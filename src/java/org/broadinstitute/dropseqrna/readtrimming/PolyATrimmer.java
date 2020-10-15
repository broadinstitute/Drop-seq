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
package org.broadinstitute.dropseqrna.readtrimming;

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.ClippingUtility;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;

@CommandLineProgramProperties(summary = "", oneLineSummary = "", programGroup = DropSeq.class)
public class PolyATrimmer extends CommandLineProgram {

	private final Log log = Log.getInstance(PolyATrimmer.class);

	// In debug mode, print a message if there is this much adapter match but no
	// poly A found.
	private static final int NO_POLY_A_ADAPTER_DEBUG_THRESHOLD = 6;

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.")
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM file")
	public File OUTPUT;

	@Argument(doc = "The output summary statistics", optional = true)
	public File OUTPUT_SUMMARY;

	@Argument(shortName = "NEW", doc="The old polyA trimmer looks for the longest run of As that is at least NUM_BASES" +
			" long and has no more than MISMATCHES bases that are not A.  For the new polyA trimmer, it is assumed that" +
			"the polyA run either extends to the end of the read, or to the beginning of adapter sequence at the end of " +
			"the read.  If the entire read is adapter sequence, the entire read is trimmed by setting all base qualities low. " +
			"If adapter sequence was found, the polyA run must be at least MIN_POLY_A_LENGTH long.  If no adapter sequence" +
			" was found, the polyA run must be at least MIN_POLY_A_LENGTH_NO_ADAPTER_MATCH long, on the assumption " +
			"that the polyA tail may extend beyond the end of the read.  The fraction of allowed non-A bases in a polyA " +
			"run must be <= MAX_POLY_A_ERROR_RATE.")
	public boolean USE_NEW_TRIMMER = false;

	@Argument(doc = "How many mismatches are acceptable in the sequence (old trim algo).")
	public Integer MISMATCHES = 0;

	@Argument(doc = "How many bases of polyA qualifies as a run of A's (old trim algo).")
	public Integer NUM_BASES = 6;

	@Argument(doc = "The tag to set for trimmed reads.  This tags the first base to exclude in the read.  37 would mean to retain the first 36 bases.")
	public String TRIM_TAG = "ZP";

	@Argument(doc = "Symbolic & literal specification of adapter sequence.  This is a combination of fixed bases to match, "
			+ " and references to SAMRecord tag values.  "
			+ "E.g. '~XM^XCACGT' means 'RCed value of XM tag' + 'value of XC tag' + 'ACGT'. "
			+ "Ideally this is at least as long as the read (new trim algo)")
	public AdapterDescriptor ADAPTER = new AdapterDescriptor(AdapterDescriptor.DEFAULT_ADAPTER);

	@Argument(doc = "Fraction of bases that can mismatch when looking for adapter match  (new trim algo)")
	public double MAX_ADAPTER_ERROR_RATE = ClippingUtility.MAX_ERROR_RATE;

	@Argument(doc = "Minimum number of bases for adapter match (new trim algo)")
	public int MIN_ADAPTER_MATCH = 4;

	@Argument(doc = "Minimum length of a poly A run, except when start of end of read intervenes (new trim algo)")
	public int MIN_POLY_A_LENGTH = 20;

	@Argument(doc = "Minimum length of poly A run at end of read, if there is no adapter match (new trim algo)")
	public int MIN_POLY_A_LENGTH_NO_ADAPTER_MATCH = 6;

	@Argument(doc = "If adapter match is at end of read, with fewer than this many bases matching the read, and not enough "
			+ "poly A is found preceding it, then ignore the adapter match and try again from the end of the read (new trim algo)")
	public int DUBIOUS_ADAPTER_MATCH_LENGTH = 6;

	@Argument(doc = "When looking for poly A, allow this fraction of bases not to be A (new trim algo)")
	public double MAX_POLY_A_ERROR_RATE = 0.1;

	private Integer readsTrimmed = 0;
	private int readsCompletelyTrimmed = 0;
	final private Histogram<Integer> numBasesTrimmed = new Histogram<>();

	// The following are only used if LOG_LEVEL==DEBUG
	int numDiffs = 0;
	int numOldDidntClip = 0;
	int numNewDidntClip = 0;

	@Override
	protected int doWork() {

		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		final ProgressLogger progress = new ProgressLogger(log);

		final SamReader bamReader = SamReaderFactory.makeDefault().open(INPUT);
		final SAMFileHeader header = bamReader.getFileHeader();
		SamHeaderUtil.addPgRecord(header, this);
		final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);
		PolyAFinder simplePolyAFinder = new SimplePolyAFinder(this.NUM_BASES, this.MISMATCHES);
		PolyAFinder polyAWithAdapterFinder = new PolyAWithAdapterFinder(ADAPTER, MIN_ADAPTER_MATCH,
				MAX_ADAPTER_ERROR_RATE, MIN_POLY_A_LENGTH, MIN_POLY_A_LENGTH_NO_ADAPTER_MATCH, MAX_POLY_A_ERROR_RATE,
				DUBIOUS_ADAPTER_MATCH_LENGTH);
		final PolyAFinder polyAFinder;
		if (USE_NEW_TRIMMER)
			polyAFinder = polyAWithAdapterFinder;
		else
			polyAFinder = simplePolyAFinder;
		for (SAMRecord r : bamReader) {
			final SimplePolyAFinder.PolyARun polyARun = polyAFinder.getPolyAStart(r);
			final int polyAStart = polyARun.startPos;

			if (log.isEnabled(Log.LogLevel.DEBUG)) {
				final PolyAFinder.PolyARun simple;
				final PolyAFinder.PolyARun withAdapter;
				if (USE_NEW_TRIMMER) {
					withAdapter = polyARun;
					simple = simplePolyAFinder.getPolyAStart(r);
				} else {
					simple = polyARun;
					withAdapter = polyAWithAdapterFinder.getPolyAStart(r);
				}
				logTrimDifference(simple, withAdapter, r);
			}

			hardClipPolyAFromRecord(r, polyAStart);
			writer.addAlignment(r);
			progress.record(r);
		}
		CloserUtil.close(bamReader);
		writer.close();
		log.info("Total " + progress.getCount() + " reads processed.");
		log.info("Number of reads trimmed: ", this.readsTrimmed);
		log.info("Number of reads completely trimmed: ", this.readsCompletelyTrimmed);
		log.debug(String.format("differences: %d; old didn't clip: %d; new didn't clip: %d", numDiffs, numOldDidntClip,
				numNewDidntClip));
		if (this.OUTPUT_SUMMARY != null)
			writeSummary(this.numBasesTrimmed);

		return 0;
	}

	private void logTrimDifference(final PolyAFinder.PolyARun simpleRun, final PolyAFinder.PolyARun withAdapterRun,
			final SAMRecord r) {
		final String newAdapter = ADAPTER.getAdapterSequence(r);
		if (simpleRun.startPos != withAdapterRun.startPos) {
			++numDiffs;
			if (withAdapterRun.isNoMatch())
				++numNewDidntClip;
			else if (simpleRun.isNoMatch())
				++numOldDidntClip;
			System.out.println("\nREADNAME: " + r.getReadName());
			final String readString = r.getReadString();
			System.out.println("READ:" + readString);
			System.out.print("OLD: ");
			if (simpleRun.isNoMatch())
				System.out.println("NOCLIP");
			else {
				indent(System.out, simpleRun.startPos);
				System.out.print("^");
				indent(System.out, simpleRun.length - 1);
				System.out.println(readString.substring(simpleRun.endPos() + 1));
			}
			if (withAdapterRun.isNoMatch()) {
				System.out.print("XEW: ");
				indent(System.out, simpleRun.startPos + simpleRun.length);
				System.out.println(newAdapter);
				if (withAdapterRun.adapterStartPos != SimplePolyAFinder.NO_MATCH
						&& readString.length() - withAdapterRun.adapterStartPos >= NO_POLY_A_ADAPTER_DEBUG_THRESHOLD) {
					System.out.print("NPA: ");
					indent(System.out, withAdapterRun.adapterStartPos);
					System.out.println(newAdapter);
				}
			} else {
				System.out.print("NEW:");
				// If no poly A (adapter starts at beginning of read, put caret
				// before the read starts
				// so that adapter aligns properly.
				if (withAdapterRun.length > 0)
					System.out.print(" ");
				indent(System.out, withAdapterRun.startPos);
				System.out.print("^");
				indent(System.out, withAdapterRun.length - 1);
				System.out.println(newAdapter);
				int positionAfterNewAdapter = withAdapterRun.endPos() + 1 + newAdapter.length();
				if (positionAfterNewAdapter < r.getReadLength()) {
					System.out.print("ANA: ");
					indent(System.out, positionAfterNewAdapter);
					System.out.println(readString.substring(positionAfterNewAdapter));
				}
			}
			/*
			 * Code below prints out adapter sequence and stuff after adapter,
			 * in case where old and new algorithms agree } else if
			 * (log.isEnabled(Log.LogLevel.DEBUG) && !polyARun.isNoMatch() &&
			 * polyARun.endPos() < readString.length() -
			 * ADAPTER_MATCH_LENGTH_DEBUG_THRESHOLD) { // Print if at least 4
			 * bases of adapter match System.out.println("\nREADNAME: " +
			 * r.getReadName()); System.out.println("READ:" + readString);
			 * System.out.print("ADAP:"); indent(System.out, withAdapterRun);
			 * System.out.print("^"); indent(System.out, newPolyARun.length -
			 * 1); System.out.println(ADAPTER.getAdapterSequence(r)); int
			 * positionAfterNewAdapter = newPolyARun.endPos() + 1 +
			 * newAdapter.length(); if (positionAfterNewAdapter <
			 * r.getReadLength()) { System.out.print("ANA: ");
			 * indent(System.out, positionAfterNewAdapter);
			 * System.out.println(readString.substring(positionAfterNewAdapter))
			 * ; }
			 */
		}

	}

	private void indent(final PrintStream out, final int amount) {
		for (int i = 0; i < amount; ++i)
			out.print(" ");
	}

	/**
	 * Hard clip out reads that have polyA runs. Finds the longest sequence of
	 * polyAs, and clips all bases that occur in or after that.
	 *
	 * @param r
	 *            The read to trim
	 * @param polyAStart
	 *            Where to clip, or -1 if no clipping. and with at most some
	 *            number of errors.
	 */
	void hardClipPolyAFromRecord(final SAMRecord r, final int polyAStart) {

		int readLength = r.getReadLength();

		// short circuit if the template wasn't found.
		if (polyAStart == SimplePolyAFinder.NO_MATCH)
			return;
		this.readsTrimmed++;

		this.numBasesTrimmed.increment(polyAStart);
		// terrible luck. your read is just a wad of A's.
		if (polyAStart == 0) {
			// attempt a work around for reads that would be 0 length after
			// trimming.
			// instead of trimming the barcode to a 0 length read, set the base
			// qualities to be low.
			byte[] value = new byte[readLength];
			Arrays.fill(value, (byte) 3);
			r.setBaseQualities(value);
			++this.readsCompletelyTrimmed;
			return;
		}

		byte[] read = r.getReadBases();
		read = Arrays.copyOfRange(read, 0, polyAStart);
		r.setReadBases(read);
		// String after = r.getReadString();

		byte[] quality = r.getBaseQualities();
		quality = Arrays.copyOfRange(quality, 0, polyAStart);
		r.setBaseQualities(quality);
		r.setAttribute(TRIM_TAG, polyAStart + 1);
	}

	private void writeSummary(final Histogram<Integer> h) {

		MetricsFile<TrimMetric, Integer> mf = new MetricsFile<>();
		mf.addHistogram(h);
		TrimMetric tm = new TrimMetric(h);
		mf.addMetric(tm);
		mf.write(this.OUTPUT_SUMMARY);
	}

	public class TrimMetric extends MetricBase {
		public Double mean;
		public Double stdev;

		public TrimMetric(final Histogram<Integer> h) {
			mean = h.getMean();
			stdev = h.getStandardDeviation();
		}

		public Double getMean() {
			return mean;
		}

		public Double getStdev() {
			return stdev;
		}
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new PolyATrimmer().instanceMain(args));
	}

}
