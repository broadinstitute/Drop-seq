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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;

import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.readpairs.ReadPair;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(usage="Adds a BAM tag to every read of the defined range of bases of the sequence of the 1st or 2nd read.  " +
        "Reads must be paired for this program to run.",
        usageShort = "Moves specified bases of each read pair into a tag",
        programGroup = DropSeq.class)
public class TagBamWithReadSequenceExtended extends CommandLineProgram {


	// ADD ability to change the tag name.
	// ADD ability to clip bases off the front of the read and make them a tag instead of discarding the read

	private final Log log = Log.getInstance(TagBamWithReadSequenceExtended.class);

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.")
	public File INPUT;

	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output bam")
	public File OUTPUT;

	@Option(doc = "Summary of barcode base quality",optional=true)
	public File SUMMARY;

	@Option(doc="Base range to extract, seperated by a dash.  IE: 1-4.  Can extract multiple ranges by seperating them by a colon.  For example 1-4:17-22 extracts the first 4 bases, then the 17-22 bases, and glues the sequence together into a single sequence for a tag.")
	public String BASE_RANGE;

	@Option(doc = "The sequence can be from the first or second read [1/2].  ")
	public Integer BARCODED_READ;

	@Option(doc="Add the tag to the sequence the read came from? If false, the read that does not have the barcode gets the tag.  If true, set the tag on the barcoded read.")
	public Boolean TAG_BARCODED_READ=false;

	@Option(doc = "Discard the read the sequence came from?.  If this is true, then the remaining read is marked as unpaired.  If the read is unpaired, then you can't discard a read.")
	public Boolean DISCARD_READ=false;

	@Option (doc="Should the bases selected for the tag be hard clipped from the read?  BE VERY CAREFUL WITH THIS FEATURE, FOR EXPERTS ONLY.  NOT NEEDED FOR STANDARD DROPSEQ DATA PROCESSING." +
	"Don't use on aligned data, does NOT change cigar strings")
	public Boolean HARD_CLIP_BASES=false;

	@Option (doc="Minimum base quality required for barcode")
	public Integer BASE_QUALITY=10;

	@Option (doc="Number of bases below minimum base quality to fail the barcode.")
	public Integer NUM_BASES_BELOW_QUALITY=1;

	@Option (doc="Barcode tag.  This is typically X plus one more capitalized alpha.  For example, 'XS', which is the default.")
	public String TAG_NAME="XS";

	private String TAG_QUALITY="XQ";

	private FailedBaseMetric metric = null;



	@Override
	protected int doWork() {
		if (this.TAG_BARCODED_READ && this.DISCARD_READ) {
			log.error("If TAG_BARCODED_READ=true and DISCARD_READ=true, you're throwing away the tag with the read. Stopping");
			System.exit(1);
		}
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);

		// get the header.
		SamReader inputSam = SamReaderFactory.makeDefault().open(INPUT);
		SAMFileHeader h= inputSam.getFileHeader();
		PeekableIterator<SAMRecord> iter = new PeekableIterator<SAMRecord>(CustomBAMIterators.getQuerynameSortedRecords(inputSam));

		SamHeaderUtil.addPgRecord(h, this);
		SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(h, true, OUTPUT);

		List<BaseRange> baseRanges = BaseRange.parseBaseRange(this.BASE_RANGE);

		this.metric = new FailedBaseMetric(BaseRange.getTotalRangeSize(this.BASE_RANGE));

		ProgressLogger progress = new ProgressLogger(this.log);

		while (iter.hasNext()) {
			//log.info(count);
			SAMRecord r1 = iter.next();
			SAMRecord r2 = iter.peek();

			// if the 2nd read is null, then you're just about done, so don't test the non-existent read name
			boolean sameName=false;
			if (r2!=null)
				sameName=r1.getReadName().equals(r2.getReadName());

			if (!sameName) {
				processSingleRead(r1, baseRanges, writer, this.HARD_CLIP_BASES);
				continue;
			}

			// check to see if the two reads are properly paired if they have the same name
			ReadPair p = new ReadPair(r1, r2);
			if (p.testProperlyPaired()==false) {
				log.error(("Reads not properly paired! R1: " + r1.getReadName() + " R2: " + r2.getReadName()));
				System.exit(1);

			}
			// since you're in paired end land, make the 2nd read a real read and not a peeked read.
			r2=iter.next();
			r1=p.getRead1();
			r2=p.getRead2();
			if (BARCODED_READ==1)
				processReadPair(r1, r2, baseRanges, writer, this.DISCARD_READ);
			if (BARCODED_READ==2)
				processReadPair(r2, r1, baseRanges, writer, this.DISCARD_READ);
			progress.record(r1);
			progress.record(r2);

		}
		writer.close();
		if (this.SUMMARY!=null) writeOutput (this.metric, this.SUMMARY);
		CloserUtil.close(inputSam);
		CloserUtil.close(iter);
		return (0);
	}

	void processSingleRead(final SAMRecord barcodedRead, final List<BaseRange> baseRanges, final SAMFileWriter writer, final boolean hardClipBases) {
		int numBadBases=scoreBaseQuality(barcodedRead, baseRanges);
		String seq = barcodedRead.getReadString();
		// does this have an off by 1 error?  I think it's 0 based so should be ok.
		//seq=seq.substring(baseRange.get, numBases);
		seq=BaseRange.getSequenceForBaseRange(baseRanges, seq);

		if (numBadBases>=this.NUM_BASES_BELOW_QUALITY) {

			// if there's an old quality setting you need to add to it instead of overwriting it.
			Object o = barcodedRead.getAttribute(this.TAG_QUALITY);
			if (o!=null) {
				int oldNumBadBases=(Integer) o;
				numBadBases+=oldNumBadBases;
			}

			barcodedRead.setAttribute(this.TAG_QUALITY, numBadBases);
		}
		barcodedRead.setAttribute(TAG_NAME, seq);
		SAMRecord result = barcodedRead;
		if (hardClipBases) result = hardClipBasesFromRead(barcodedRead, baseRanges);
		writer.addAlignment(result);
	}

	static SAMRecord hardClipBasesFromRead (final SAMRecord r, final List<BaseRange> baseRanges) {

		int readLength=r.getReadLength();

		List<BaseRange> basesToKeep = BaseRange.invert(baseRanges, readLength);
		byte [] newSequence = BaseRange.getBytesForBaseRange(basesToKeep,  r.getReadBases());
		byte [] newBaseQualities = BaseRange.getBytesForBaseRange(basesToKeep,  r.getBaseQualities());

		r.setReadBases(newSequence);
		r.setBaseQualities(newBaseQualities);
		return r;
	}

	private SAMRecord setTagsOnRead (final SAMRecord r, int numBadBases, final String seq) {
		if (numBadBases>=this.NUM_BASES_BELOW_QUALITY) {
			// if there's an old quality setting you need to add to it instead of overwriting it.
			Object o = r.getAttribute(this.TAG_QUALITY);
			if (o!=null) {
				int oldNumBadBases=(Integer) o;
				numBadBases+=oldNumBadBases;
			}
			r.setAttribute(this.TAG_QUALITY, numBadBases);
		}
		r.setAttribute(TAG_NAME, seq);
		return (r);
	}

	void processReadPair (final SAMRecord barcodedRead, final SAMRecord otherRead, final List<BaseRange> baseRanges, final SAMFileWriter writer, final boolean discardRead) {
		int numBadBases=scoreBaseQuality(barcodedRead, baseRanges);
		String seq = barcodedRead.getReadString();
		// does this have an off by 1 error?  I think it's 0 based so should be ok.

		//seq=seq.substring(0, numBases);
		seq=BaseRange.getSequenceForBaseRange(baseRanges, seq);
		if (this.TAG_BARCODED_READ)
			setTagsOnRead(barcodedRead, numBadBases, seq);
		else
			setTagsOnRead(otherRead, numBadBases, seq);

		if (discardRead) {
			int flag =otherRead.getFlags();
			if (otherRead.getMateUnmappedFlag()) flag-=8;
			if (otherRead.getMateNegativeStrandFlag()) flag-=32;
			if (otherRead.getReadPairedFlag()) flag-=1;
			if (otherRead.getFirstOfPairFlag()) flag-=64;
			if (otherRead.getSecondOfPairFlag()) flag-=128;
			otherRead.setFlags(flag);

		} else
			writer.addAlignment(barcodedRead);
		writer.addAlignment(otherRead);
	}

	private int scoreBaseQuality(final SAMRecord barcodedRead, final List<BaseRange> baseRanges) {
		int numBasesBelowQuality=0;
		byte [] qual= barcodedRead.getBaseQualities();
		char [] seq = barcodedRead.getReadString().toUpperCase().toCharArray();
		for (BaseRange b: baseRanges)
			for (int i=b.getStart()-1; i<b.getEnd(); i++) {
				byte q = qual[i];
				char s = seq[i];

				if (q < this.BASE_QUALITY || s=='N')
					numBasesBelowQuality++;
			}

		this.metric.addFailedBase(numBasesBelowQuality);
		return (numBasesBelowQuality);
	}

	private class FailedBaseMetric {
		List<Integer> data = null;

		public FailedBaseMetric (final Integer length){
			data=new ArrayList<Integer>(length+1);
			for (int i=0; i<=length; i++)
				data.add(new Integer(0));
		}

		public void addFailedBase(final int numBasesFailed) {
			Integer i = data.get(numBasesFailed);
			i++;
			data.set(numBasesFailed, i);
		}

		public int getNumFailedBases(final int position) {
			return (data.get(position));
		}

		public int getLength() {
			return (data.size());
		}

	}

	private void writeOutput (final FailedBaseMetric result, final File output) {
		BufferedWriter writer = OutputWriterUtil.getWriter(output);
		String [] header = {"num_failed_bases", "num_barcodes"};
		String h = StringUtils.join(header, "\t");
		OutputWriterUtil.writeResult(h, writer);
		for (int i=0; i<result.getLength(); i++) {
			int count=result.getNumFailedBases(i);
			String [] l={i+"", count+""};
			String line = StringUtils.join(l, "\t");
			OutputWriterUtil.writeResult(line, writer);
		}
		OutputWriterUtil.closeWriter(writer);
	}


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new TagBamWithReadSequenceExtended().instanceMain(args));
	}
}

