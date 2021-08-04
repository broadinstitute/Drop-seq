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
package org.broadinstitute.dropseqrna.sbarro;

import java.io.File;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.Sbarro;
import org.broadinstitute.dropseqrna.utils.BaseRange;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

/**
 * For a BAM file full of rabies virus tags, subset the BAM into reads that either pass or fail the filter set, then write them into the appropriate BAMs.
 * @author nemesh
 *
 */
@CommandLineProgramProperties(
        summary = "Filters rabies virus tags into pass/fail BAMs based on metrics associated with the extracted viral barcodes.",
        oneLineSummary = "Filters rabies virus tags into pass/fail BAMs",
        programGroup = Sbarro.class)
public class FilterValidRabiesBarcodes extends CommandLineProgram {

	private final Log log = Log.getInstance(FilterValidRabiesBarcodes.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "A BAM file with rabies barcode tags")
	public File INPUT;

	@Argument(doc="A BAM with accepted rabies barcode tags.", optional=true)
	public File OUTPUT_ACCEPTED;

	@Argument(doc="A BAM with rejected rabies barcode tags.", optional=true)
	public File OUTPUT_REJECTED;

	@Argument (doc="The BAM tag defining the edit distance to the discovered GFP anchor")
	public String GFP_ANCHOR_ED_TAG="ga";
	@Argument (doc="The BAM tag defining the edit distance to the discovered cassette anchor")
	public String CASSETTE_ANCHOR_ED_TAG="ca";
	@Argument (doc="The rabies barcode stop codon tag")
	public String STOP_CODON_BARCODE_TAG="sb";
	@Argument (doc="The rabies barcode polyA tag")
	public String POLY_A_BARCODE_TAG="pb";

	// the values of the filtering parameters.
	@Argument (doc="Stop barcode length range.  This field is dash seperated from the shortest to the longest acceptable length.  For example, 9-11 would accept all barcodes between 9 and 11 bases in length.",optional=true)
	public String STOP_BARCODE_LENGTH_RANGE;

	@Argument (doc="PolyA barcode length range.  This field is dash seperated from the shortest to the longest acceptable length.  For example, 9-11 would accept all barcodes between 9 and 11 bases in length.", optional=true)
	public String POLYA_BARCODE_LENGTH_RANGE;

	@Argument (doc="The maximum edit distance of the GFP anchor.", optional=true)
	public Integer MAX_GFP_ANCHOR_EDIT_DISTANCE;

	@Argument (doc="The maximum edit distance of the GFP anchor.", optional=true)
	public Integer MAX_CASSETTE_ANCHOR_EDIT_DISTANCE;
	
	@Argument (doc="The maximum number of N's in the sequence of the STOP_CODON_BARCODE_TAG or POLY_A_BARCODE_TAG.  Reads with more N bases in their barcodes than this threshold will be rejected")
	public Integer MAX_N_IN_BARCODE=null;

	@Override
	protected int doWork() {
		// set up input.
		IOUtil.assertFileIsReadable(this.INPUT);
		SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
		SAMFileHeader header = reader.getFileHeader();
		SamHeaderUtil.addPgRecord(header, this);

		SAMFileWriter writerAccepted=null;
		SAMFileWriter writerRejected=null;

		if (OUTPUT_ACCEPTED!=null) {
			IOUtil.assertFileIsWritable(this.OUTPUT_ACCEPTED);
			writerAccepted= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT_ACCEPTED);
		}

		if (OUTPUT_REJECTED!=null) {
			IOUtil.assertFileIsWritable(this.OUTPUT_REJECTED);
			writerRejected= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT_REJECTED);
		}

		BaseRange stopBCLengthRange = BaseRange.parseSingleBaseRange(STOP_BARCODE_LENGTH_RANGE);
		BaseRange polyABCLengthRange = BaseRange.parseSingleBaseRange(POLYA_BARCODE_LENGTH_RANGE);

		ProgressLogger pl = new ProgressLogger(log);

		for (SAMRecord r: reader) {
			pl.record(r);
			boolean acceptRead = acceptRead(r, this.GFP_ANCHOR_ED_TAG, this.CASSETTE_ANCHOR_ED_TAG, this.STOP_CODON_BARCODE_TAG, this.POLY_A_BARCODE_TAG, stopBCLengthRange, polyABCLengthRange, this.MAX_GFP_ANCHOR_EDIT_DISTANCE, this.MAX_CASSETTE_ANCHOR_EDIT_DISTANCE, this.MAX_N_IN_BARCODE);
			if (acceptRead && writerAccepted!=null) writerAccepted.addAlignment(r);
			if (!acceptRead && writerRejected!=null) writerRejected.addAlignment(r);
		}

		CloserUtil.close(reader);
		if (writerAccepted!=null) CloserUtil.close(writerAccepted);
		if (writerRejected!=null) CloserUtil.close(writerRejected);

		return 0;
	}

	boolean acceptRead (final SAMRecord r, final String gfpAnchorEDTag, final String cassetteAnchorEDTag, final String stopCodonBCTag, final String polyABCTag,
			final BaseRange stopBCRange, final BaseRange polyABCRange, final Integer maxGFPED, final Integer maxCassetteED, final Integer maxN) {

		// can test the maxN in sequence ahead of the other tags.
		if (maxN!=null) {
			Integer countStopCodonBC=getNBasesInBarcode(r, stopCodonBCTag);
			Integer countPolyABC = getNBasesInBarcode(r, polyABCTag);
			if (countStopCodonBC!=null && countStopCodonBC>maxN) return false;
			if (countPolyABC!=null && countStopCodonBC>maxN) return false;
		}
		// get all the info off the read.  Fields can be null!
		Integer readGFPAnchorED=r.getIntegerAttribute(gfpAnchorEDTag);
		Integer readCassetteAnchorED = r.getIntegerAttribute(cassetteAnchorEDTag);
		Integer stopCodonBCSize=getBarcodeLength(r, stopCodonBCTag);
		Integer polyABCSize = getBarcodeLength(r, polyABCTag);
		boolean flag = acceptRead (readGFPAnchorED, readCassetteAnchorED, stopCodonBCSize, polyABCSize, stopBCRange, polyABCRange, maxGFPED, maxCassetteED);
		return flag;
	}
	
	private Integer getNBasesInBarcode (final SAMRecord r, final String tag) {
		String seq = r.getStringAttribute(tag);
		if (seq==null) return null;
		return StringUtils.countMatches(seq, "N");		
	}

	private Integer getBarcodeLength (final SAMRecord r, final String tag) {
		String v = r.getStringAttribute(tag);
		if (v==null) return null;
		// figure out how long the string is.
		return new Integer (v.length());
	}

	/**
	 * Should this read be accepted?
	 * For a read to be accepted, the parameters with values (readGFPAnchorED,readCassetteAnchorED,  stopCodonBCSize, polyABCSize) must be non-null.
	 * The parameters can be null, in which case they are ignored.
	 * @param readGFPAnchorED
	 * @param readCassetteAnchorED
	 * @param stopCodonBCSize
	 * @param polyABCSize
	 * @param stopBCRange
	 * @param polyABCRange
	 * @param maxGFPED
	 * @param maxCassetteED
	 * @return
	 */
	public static boolean acceptRead (final Integer readGFPAnchorED, final Integer readCassetteAnchorED, final Integer stopCodonBCSize, final Integer polyABCSize,
			final BaseRange stopBCRange, final BaseRange polyABCRange, final Integer maxGFPED, final Integer maxCassetteED) {

		// validate all parameters are set.
		if (readGFPAnchorED==null || readCassetteAnchorED==null || stopCodonBCSize==null || polyABCSize==null) return false;

		// test GFP anchor edit distance.
		if (maxGFPED!=null && readGFPAnchorED>maxGFPED) return false;
		if (maxCassetteED!=null && readCassetteAnchorED>maxCassetteED) return false;
		// test barcode length
		if (stopBCRange!=null && (stopCodonBCSize<stopBCRange.getStart() || stopCodonBCSize>stopBCRange.getEnd())) return false;
		if (polyABCRange!=null && (polyABCSize<polyABCRange.getStart() || polyABCSize>polyABCRange.getEnd())) return false;

		return true;
	}


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new FilterValidRabiesBarcodes().instanceMain(args));
	}



}
