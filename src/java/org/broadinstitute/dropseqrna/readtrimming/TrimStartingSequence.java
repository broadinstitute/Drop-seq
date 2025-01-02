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
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

@CommandLineProgramProperties(summary = "Trim the given sequence from the beginning of reads",
        oneLineSummary = "Trim the given sequence from the beginning of reads",
        programGroup = DropSeq.class)
public class TrimStartingSequence extends AbstractTrimmerClp {
	public static final String DEFAULT_TRIM_TAG = "ZS";

	private final Log log = Log.getInstance(TrimStartingSequence.class);

	@Argument(doc="The sequence to look for at the start of reads.")
	public String SEQUENCE;

	@Argument(doc="How many mismatches are acceptable in the sequence.  " +
			"If neither MISMATCHES nor MISMATCH_RATE is specified, default behavior is MISMATCHES=0",
			optional = true, mutex = {"MISMATCH_RATE"})
	public Integer MISMATCHES;

	@Argument(doc="What fraction of bases the matched sequence can mismatch.  Must be >=0 and <1.  " +
			"In contrast to MISMATCHES, this matcher will match the full sequence even if it is preceded by something " +
			"else in the read.",
			optional = true, mutex = {"MISMATCHES", "LEGACY"})
	public Double MISMATCH_RATE;

	@Argument(doc="Enable the old trim algorithm (release <= 2.4.0), which did not match if bases precede the sequence, " +
			"and had bugs if sequence was longer than read length.  ", mutex = {"MISMATCH_RATE", "LENGTH_TAG"})
	public boolean LEGACY = false;

	@Argument(doc="How many bases at the beginning of the sequence must match before trimming occurs.")
	public Integer NUM_BASES=0;

	@Argument (doc="The tag to set for trimmed reads.  This tags the first base to keep in the read.  6 would mean to trim the first 5 bases.")
	public String TRIM_TAG= DEFAULT_TRIM_TAG;

	@Argument(doc="Tag containing the length of sequence matched.  If using MISMATCHES algorithm, this will be the " +
			"same value as stored in TRIM_TAG.  If using MISMATCH_RATE, full-length sequence will match even if " +
			"something precedes it, so this may be different than TRIM_TAG value.  Not stored if not set.",
			optional = true, mutex = {"LEGACY"})
	public String LENGTH_TAG;

	@Override
	protected int doWork() {
		if (MISMATCHES == null && MISMATCH_RATE == null) {
			MISMATCHES = 0;
		}

		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		final ProgressLogger progress = new ProgressLogger(log);

		SamReader bamReader = SamReaderFactory.makeDefault().open(this.INPUT);
		SAMFileHeader header = bamReader.getFileHeader();
		SamHeaderUtil.addPgRecord(header, this);
        SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);

        if (LEGACY) {
			TrimSequenceTemplate t = new TrimSequenceTemplate(this.SEQUENCE);

			for (SAMRecord r: bamReader) {
				SAMRecord rr =  shouldTrim(r)? hardClipBarcodeFromRecordLegacy(r, t, this.NUM_BASES, this.MISMATCHES): r;
				writer.addAlignment(rr);
				progress.record(r);
				this.numReadsTotal++;
			}

		} else {
			final StartingSequenceTrimmer trimmer;
			if (MISMATCHES != null) {
				trimmer = new FixedMismatchStartingSequenceTrimmer(SEQUENCE, NUM_BASES, MISMATCHES);
			} else {
				trimmer = new MismatchRateStartingSequenceTrimmer(SEQUENCE, NUM_BASES, MISMATCH_RATE);
			}

			for (SAMRecord r : bamReader) {
				SAMRecord rr = shouldTrim(r)? hardClipBarcodeFromRecord(r, trimmer): r;
				writer.addAlignment(rr);
				progress.record(r);
				this.numReadsTotal++;
			}
		}

        CloserUtil.close(bamReader);

        writer.close();
        log.info("Number of reads trimmed: " + this.readsTrimmed, " total reads: " + this.numReadsTotal);
        if (this.OUTPUT_SUMMARY!=null) writeSummary(this.numBasesTrimmed);

		return 0;
	}

	@Override
	protected String[] customCommandLineValidation() {
		final ArrayList<String> list = new ArrayList<>(1);
		if (MISMATCH_RATE != null && (MISMATCH_RATE < 0 || MISMATCH_RATE >= 1)) {
			list.add("MISMATCH_RATE must be >= 0 and < 1");
		}
		if (!VALID_WHICH_READ.containsAll(WHICH_READ)) {
			list.add("WHICH_READ must be one of " + VALID_WHICH_READ);
		}
		if (WHICH_READ.isEmpty()) {
			list.add("WHICH_READ must be specified");
		}
		return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
	}

	private void writeSummary (final Histogram<Integer> h) {

		MetricsFile<TrimMetric, Integer> mf = new MetricsFile<>();
		mf.addHistogram(h);
		TrimMetric tm=new TrimMetric(h);
		mf.addMetric(tm);
		mf.write(this.OUTPUT_SUMMARY);
	}

	public static class TrimMetric extends MetricBase {
		public Double mean;
		public Double stdev;

		public TrimMetric (final Histogram<Integer> h) {
			mean=h.getMean();
			stdev=h.getStandardDeviation();
		}


		public Double getMean() {
			return mean;
		}

		public Double getStdev() {
			return stdev;
		}


	}




	/**
	 * Hard clip out a sequence from the start of a record.
	 * Looks for bases at the start of a read that match part of the template, and removes those bases.
	 * At least <minMatch> bases at the start of a read must match the template to be trimmed
	 * Mismatches depends on algorithm selected.
	 * @param r The read to trim
	 * @param trimmer Trimming implementation
	 */
	SAMRecord hardClipBarcodeFromRecord (final SAMRecord r, final StartingSequenceTrimmer trimmer) {
		final StartingSequenceTrimmer.TrimResult trimResult = trimmer.trim(r.getReadString());
		// short circuit if the template wasn't found.
		if (!trimResult.wasTrim())
			return (r);
		this.readsTrimmed++;

		this.numBasesTrimmed.increment(trimResult.endPosition);

		if (trimResult.completelyTrimmed) {
			// attempt a work around for reads that would be 0 length after trimming.
			// instead of trimming the barcode to a 0 length read, set the base qualities to be low.
			byte [] value= new byte [r.getReadLength()];
			Arrays.fill(value, (byte) 3);
			r.setBaseQualities(value);
		} else {

			byte[] read = r.getReadBases();
			read = Arrays.copyOfRange(read, trimResult.endPosition, read.length);
			r.setReadBases(read);

			byte[] quality = r.getBaseQualities();
			quality = Arrays.copyOfRange(quality, trimResult.endPosition, quality.length);
			r.setBaseQualities(quality);
			r.setAttribute(this.TRIM_TAG, trimResult.endPosition);
		}
		// This tag will be set regardless of whether completely or partially trimmed
		if (this.LENGTH_TAG != null) {
			r.setAttribute(this.LENGTH_TAG, trimResult.endPosition - trimResult.startPosition);
		}
		return (r);
	}

	/**
	 * Old implementation for backward compatibility.
	 * Hard clip out a sequence from the start of a record.
	 * Looks for bases at the start of a read that match part of the template, and removes those bases.
	 * At least <minMatch> bases at the start of a read must match the template to be trimmed
	 * At most <mismatchesAllowed> errors can exist between the template and the read for a trim to occur.
	 * @param r The read to trim
	 * @param t The template sequence to search for a trim
	 * @param minMatch
	 * @param mismatchesAllowed
	 * @return
	 */
	SAMRecord hardClipBarcodeFromRecordLegacy (final SAMRecord r, final TrimSequenceTemplate t, final int minMatch, final int mismatchesAllowed) {
		int templateLength= t.getSequence().length();
		int readLength=r.getReadLength();

		String readString = r.getReadString();
		int templatePosition = t.getPositionInTemplate(readString, minMatch, mismatchesAllowed);
		// short circuit if the template wasn't found.
		if (templatePosition==-1)
			return (r);
		this.readsTrimmed++;

		int firstBaseToKeep = (templateLength - templatePosition);
		this.numBasesTrimmed.increment(firstBaseToKeep);

		if (templateLength>=readLength) {
			// attempt a work around for reads that would be 0 length after trimming.
			// instead of trimming the barcode to a 0 length read, set the base qualities to be low.
			byte [] value= new byte [readLength];
			Arrays.fill(value, (byte) 3);
			r.setBaseQualities(value);
			return (r);
		}

		byte [] read = r.getReadBases();
		read=Arrays.copyOfRange(read, firstBaseToKeep, read.length);
		r.setReadBases(read);
		readString = r.getReadString();

		byte [] quality = r.getBaseQualities();
		quality=Arrays.copyOfRange(quality, firstBaseToKeep, quality.length);
		r.setBaseQualities(quality);
		r.setAttribute(this.TRIM_TAG, firstBaseToKeep);
		return (r);
	}


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new TrimStartingSequence().instanceMain(args));
	}
}
