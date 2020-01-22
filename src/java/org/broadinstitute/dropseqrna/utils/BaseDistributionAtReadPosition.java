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
package org.broadinstitute.dropseqrna.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(summary = "Reads each base and generates a composition per-position matrix",
        oneLineSummary = "Reads each base and generates a composition per-position matrix",
        programGroup = DropSeq.class)
public class BaseDistributionAtReadPosition extends CommandLineProgram {

	private final Log log = Log.getInstance(BaseDistributionAtReadPosition.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.")
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output report")
	public File OUTPUT;

	@Argument(doc="Read to gather statistics on [1/2].  If this is set, the tag is ignored.", optional=true)
	public Integer READ_NUMBER=null;

	@Argument(doc="Tag to gather statistics on.  If this is set, the read number is ignored.", optional=true)
	public String TAG = null;

	//@Argument(doc = "Minimum mapping quality to consider the read")
	// public int MINIMUM_MAPPING_QUALITY = 0;

	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		BaseDistributionMetricCollection result;

		if (this.TAG != null)
			result = gatherBaseQualities(INPUT, this.TAG);
		else
			result = gatherBaseQualities(INPUT, this.READ_NUMBER);

        result.writeOutput(OUTPUT);

		return (0);
	}

	BaseDistributionMetricCollection gatherBaseQualities (final File input, final int readNumber) {
		ProgressLogger p = new ProgressLogger(this.log);

		SamReader inputSam = SamReaderFactory.makeDefault().open(input);
		BaseDistributionMetricCollection c = new BaseDistributionMetricCollection();

		MAIN_LOOP:
		for (final SAMRecord samRecord : inputSam) {

			p.record(samRecord);
			if (samRecord.isSecondaryOrSupplementary()) continue;
			boolean readPaired = samRecord.getReadPairedFlag();

			boolean firstRead=false;
			if (!readPaired & readNumber==2)
				continue;
			else if (!readPaired & readNumber==1)
				firstRead=true;
			else
				firstRead = samRecord.getFirstOfPairFlag();

			// if you're looking for the first read and this isn't, or looking for the 2nd read and this isn't, then go to the next read.
			if ((firstRead && readNumber!=1) || (!firstRead && readNumber==1)) continue MAIN_LOOP;

			byte [] bases = samRecord.getReadBases();

			for (int i=0; i<bases.length; i++) {
				char base = (char) (bases[i]);
				int idx=i+1;
				c.addBase(base, idx);
			}
		}

		CloserUtil.close(inputSam);
		return (c);


	}

	BaseDistributionMetricCollection gatherBaseQualities (final File input, final String tag) {
		ProgressLogger pl = new ProgressLogger(this.log);
		SamReader inputSam = SamReaderFactory.makeDefault().open(input);

		BaseDistributionMetricCollection c = new BaseDistributionMetricCollection();
		// MAIN_LOOP:
		for (final SAMRecord samRecord : inputSam) {
			pl.record(samRecord);
			if (samRecord.isSecondaryOrSupplementary()) continue;
			String b = samRecord.getStringAttribute(tag);
			if (b==null) continue;

			byte [] bases = b.getBytes();
			for (int i=0; i<bases.length; i++) {
				char base = (char) (bases[i]);
				int idx=i+1;
				c.addBase(base, idx);
			}
		}

		CloserUtil.close(inputSam);
		return (c);

	}




	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new BaseDistributionAtReadPosition().instanceMain(args));
	}
}
