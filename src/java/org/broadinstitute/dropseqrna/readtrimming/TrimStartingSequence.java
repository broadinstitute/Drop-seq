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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.util.Arrays;

import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(usage = "Trim the given sequence from the beginning of reads",
        usageShort = "Trim the given sequence from the beginning of reads",
        programGroup = DropSeq.class)
public class TrimStartingSequence extends CommandLineProgram {
	private final Log log = Log.getInstance(TrimStartingSequence.class);

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.")
	public File INPUT;

	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM file")
	public File OUTPUT;

	@Option(doc = "The output summary statistics", optional=true)
	public File OUTPUT_SUMMARY;

	@Option(doc="The sequence to look for at the start of reads.")
	public String SEQUENCE;

	@Option(doc="How many mismatches are acceptable in the sequence.")
	public Integer MISMATCHES=0;

	@Option(doc="How many bases at the begining of the sequence must match before trimming occurs.")
	public Integer NUM_BASES=0;

	@Option (doc="The tag to set for trimmed reads.  This tags the first base to keep in the read.  6 would mean to trim the first 5 bases.")
	public String TRIM_TAG="ZS";

	private Integer readsTrimmed=0;
	private int numReadsTotal=0;
	private Histogram<Integer> numBasesTrimmed= new Histogram<Integer>();

	@Override
	protected int doWork() {

		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		final ProgressLogger progress = new ProgressLogger(log);

		SamReader bamReader = SamReaderFactory.makeDefault().open(this.INPUT);
		SAMFileHeader header = bamReader.getFileHeader();
		SamHeaderUtil.addPgRecord(header, this);
        SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);

        TrimSequenceTemplate t = new TrimSequenceTemplate(this.SEQUENCE);

        for (SAMRecord r: bamReader) {
        	SAMRecord rr =  hardClipBarcodeFromRecord(r, t, this.NUM_BASES, this.MISMATCHES);
        	writer.addAlignment(rr);
        	progress.record(r);
        	this.numReadsTotal++;
        }

        CloserUtil.close(bamReader);

        writer.close();
        log.info("Number of reads trimmed: " + this.readsTrimmed, " total reads: " + this.numReadsTotal);
        if (this.OUTPUT_SUMMARY!=null) writeSummary(this.numBasesTrimmed);

		return 0;
	}

	private void writeSummary (final Histogram<Integer> h) {

		MetricsFile<TrimMetric, Integer> mf = new MetricsFile<TrimMetric, Integer>();
		mf.addHistogram(h);
		TrimMetric tm=new TrimMetric(h);
		mf.addMetric(tm);
		mf.write(this.OUTPUT_SUMMARY);
	}

	public class TrimMetric extends MetricBase {
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




	// get the length of the template. (41)
	// get the position in the template. (20) (base 0)
	// bases to trim = 41-20=21  starting base to keep is 21 (base 0)
	/**
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
	SAMRecord hardClipBarcodeFromRecord (final SAMRecord r, final TrimSequenceTemplate t, final int minMatch, final int mismatchesAllowed) {
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
