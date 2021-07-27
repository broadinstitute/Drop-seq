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

import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.util.*;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropNet;
import org.broadinstitute.dropseqrna.sbarro.utils.ConsensusSequence;
import org.broadinstitute.dropseqrna.sbarro.utils.ConsensusSequenceFactory;
import org.broadinstitute.dropseqrna.sbarro.utils.ExtractBarcodeSequences;
import org.broadinstitute.dropseqrna.sbarro.utils.ExtractedSequenceGroup;
import org.broadinstitute.dropseqrna.utils.BaseQualityFilter;
import org.broadinstitute.dropseqrna.utils.BaseRange;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistance;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readpairs.ReadPair;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Tags an unaligned BAM with rabies virus sequences",
        oneLineSummary = "Tags an unaligned BAM with rabies virus sequences and associated metrics.",
        programGroup = DropNet.class)

public class TagReadWithRabiesBarcodes extends CommandLineProgram {
	private final Log log = Log.getInstance(TagReadWithRabiesBarcodes.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "A BAM file.  If GENERATE_CONSENSUS is enabled, the BAM must be queryname sorted.")
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="A BAM with rabies barcode tags.")
	public File OUTPUT;

	@Argument(doc="A rather verbose report containing the barcode extraction information for each read.", optional=true)
	public File REPORT;

	@Argument(doc="When there are 2 reads that are paired, generate the consensus of the two reads.")
	public boolean GENERATE_CONSENSUS=false;

	@Argument (doc="When GENERATE_CONSENSUS is true, should we assume the second read should be reverse complimented?  This is the normal expectation," +
			"and will save a siginficant amount of computation time by avoiding a second alignment.")
	public boolean ASSUME_REVERSE_COMPLIMENT=true;

	@Argument(doc="Include these BAM tags as columns in the report.")
	public List<String> BAM_TAG_LIST=null;

	@Argument(doc="Reads (1 or both depending on if you're generating consensus) must be this long to process")
	public Integer MIN_READ_LENGTH=20;

	@Argument (doc="The GFP anchor sequence to search for")
	public String GFP_ANCHOR_SEQUENCE="CGGCATGGACGAGCTGTACAAGTAAGCTA";

	@Argument (doc="The cassette anchor sequence to search for")
	public String CASSETTE_ANCHOR_SEQUENCE="CCGGTGGCGCCACTGC";

	@Argument (doc="The rabies barcode stop codon tag")
	public String STOP_CODON_BARCODE_TAG="sb";
	@Argument (doc="The rabies barcode polyA tag")
	public String POLY_A_BARCODE_TAG="pb";
	@Argument (doc="The rabies full barcode tag.  This is the concatonation of the stop codon barcode and polyA barcode.")
	public String RABIES_BARCODE="rb";
	@Argument (doc="The edit distance to the discovered GFP anchor")
	public String GFP_ANCHOR_ED_TAG="ga";
	@Argument (doc="The edit distance to the discovered cassette anchor")
	public String CASSETTE_AHCHOR_ED_TAG="ca";
	@Argument (doc="The edit distance of the read consensus, if running the BAM contains paried reads and GENERATE_CONSENSUS=true")
	public String READ_CONSENSUS_ED_TAG="cd";

	@Argument(doc="If you've adapter marked your reads, this is the tag the reads will use to indicate the first trimmed position.  XT is the standard Picard tag.")
	public String ADAPTER_TAG=ReservedTagConstants.XT;

	@Argument (doc="Minimum base quality required for barcode")
	public Integer BASE_QUALITY=10;

	@Argument (doc="Number of bases below minimum base quality to fail the barcode.")
	public Integer NUM_BASES_BELOW_QUALITY=1;

	@Argument (doc="A tag indicating if the fails base quality metrics.")
	public String TAG_QUALITY="XQ";

	@Argument (doc="Write a histogram of the number of bases in each rabies barcode below the BASE_QUALITY threshold to this file.", optional=true)
	public File BASE_QUALITY_REPORT=null;

	@Argument(doc="Number of threads to use.  Defaults to 1.")
	public int NUM_THREADS=1;

	private ObjectCounter<Integer> failedBaseHistogram;
	// how many reads are processed to report the result.
	private int BATCH_REPORT_SIZE=100000;

	private ForkJoinPool forkJoinPool=null;

	@Override
	protected int doWork() {
		if (this.NUM_THREADS>1) this.forkJoinPool = new ForkJoinPool(this.NUM_THREADS);
		if (this.NUM_THREADS>1 & this.GENERATE_CONSENSUS) log.warn("Consensus generating process not yet multi-threaded.");

		IOUtil.assertFileIsReadable(this.INPUT);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		if (BASE_QUALITY_REPORT!=null) IOUtil.assertFileIsWritable(this.BASE_QUALITY_REPORT);

		failedBaseHistogram = new ObjectCounter<>();

		PrintStream reportOut = null;
		if (this.REPORT!=null) {
			IOUtil.assertFileIsWritable(this.REPORT);
			reportOut =  new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.REPORT));
			writeHeader(reportOut, this.BAM_TAG_LIST);
		}

		SamReader reader = SamReaderFactory.makeDefault().open(INPUT);

		if (GENERATE_CONSENSUS && !reader.getFileHeader().getSortOrder().equals(SortOrder.queryname)) {
			log.error("This program only accepts queryname sorted BAMs if GENERATE_CONSENSUS is true");
			return 1;
		}

		// set up the BAM writer.
		SAMFileHeader header = reader.getFileHeader();
		SamHeaderUtil.addPgRecord(header, this);

		ExtractBarcodeSequences extractor = new ExtractBarcodeSequences(this.GFP_ANCHOR_SEQUENCE, this.CASSETTE_ANCHOR_SEQUENCE);

		if (GENERATE_CONSENSUS) {
			SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);
			processSingleThreaded (reader, writer, extractor, reportOut, this.BASE_QUALITY, this.TAG_QUALITY, this.NUM_BASES_BELOW_QUALITY);
		}
		else {
			SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);
			processMultiThreaded (reader, writer, extractor, reportOut, this.BASE_QUALITY, this.TAG_QUALITY, this.NUM_BASES_BELOW_QUALITY, this.NUM_THREADS);
			//
		}

		if (this.BASE_QUALITY_REPORT!=null) writeBaseQualityReport(this.failedBaseHistogram, this.BASE_QUALITY_REPORT);
		return 0;

	}

	private void writeBaseQualityReport(final ObjectCounter<Integer> failedBaseHistogram, final File outFile) {
		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outFile));
		String [] header = {"NUM_ERRORS", "NUM_READS"};
		out.println(StringUtil.join("\t", header));
		int max = Collections.max(failedBaseHistogram.getKeys());
		for (int i=0; i<=max; i++) {
			String [] body = {Integer.toString(i),Integer.toString(failedBaseHistogram.getCountForKey(i))};
			out.println(StringUtil.join("\t", body));
		}
		CloserUtil.close(out);
	}

	// EXPERIMENTAL NON_WORKING STUFF.
	/**
	 * Multi-threaded version of single read (non-consensus) processing.
	 * @param reader
	 * @param writer
	 * @param extractor
	 * @param reportOut
	 * @param baseQualityThreshold
	 * @param baseQualityTagName
	 * @param numBasesBelowQuality
	 */
	private void processMultiThreaded (final SamReader reader, final SAMFileWriter writer, final ExtractBarcodeSequences extractor,
			final PrintStream reportOut, final int baseQualityThreshold, final String baseQualityTagName, final int numBasesBelowQuality, final int numThreads) {

		// if there's only one thread, fall back to the single threaded approach.
		if (numThreads==1) {
			processSingleThreaded (reader, writer, extractor, reportOut, this.BASE_QUALITY, this.TAG_QUALITY, this.NUM_BASES_BELOW_QUALITY);
			return;
		}
		Iterator<SAMRecord> iter = reader.iterator();
		int batchSize = 1000;
		ProgressLogger pl = new ProgressLogger(log, this.BATCH_REPORT_SIZE);

		// grab a batch of reads, process, then write the results.
		// since we can't do anything about the exceptions, we just catch and print the stack.  Hopefully we never see them...
		while (true) {
			List<SAMRecord> batch = getBatch(iter, batchSize, pl);
			if (batch.isEmpty()) break;
			List<SAMRecord> result=Collections.EMPTY_LIST;
			try {
				result = forkJoinPool.submit(() -> batch.parallelStream().map(x-> processRead (x, extractor, this.BAM_TAG_LIST, reportOut, baseQualityThreshold, baseQualityTagName, numBasesBelowQuality)).collect(Collectors.toList())).get();
			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
			result.stream().forEach(writer::addAlignment);
		}
		CloserUtil.close(iter);
		CloserUtil.close(writer);
		CloserUtil.close(reportOut);

	}

	private List<SAMRecord> getBatch (final Iterator<SAMRecord> iter, final int batchSize, final ProgressLogger pl) {
		List<SAMRecord> recs = new ArrayList<>();
		int count=0;
		while (iter.hasNext() & count<batchSize) {
			count++;
			SAMRecord r1 = iter.next();
			recs.add(r1);
			pl.record(r1);
		}
		return recs;

	}

	private void processSingleThreaded (final SamReader reader, final SAMFileWriter writer, final ExtractBarcodeSequences extractor, final PrintStream reportOut, final int baseQualityThreshold, final String baseQualityTagName, final int numBasesBelowQuality) {
		PeekableIterator<SAMRecord> iter = new PeekableIterator<>(reader.iterator());

		ProgressLogger pl = new ProgressLogger(log, this.BATCH_REPORT_SIZE);

		// walk through the reads and peek one ahead to find paired reads.
		while (iter.hasNext()) {
			SAMRecord r1 = iter.next();
			pl.record(r1);
			if (iter.hasNext()) {
				SAMRecord r2 = iter.peek();
				// process paired read if reads are paired and you want to generate a consensus sequence.
				if (r1.getReadName().equals(r2.getReadName()) && this.GENERATE_CONSENSUS) {
					r2=iter.next();
					pl.record(r2);
					r1=processReadPair(r1,r2, extractor, this.BAM_TAG_LIST, reportOut);
					writer.addAlignment(r1);
					writer.addAlignment(r2);
				}
				else { // only write out the single read if not paired.
					r1= processRead (r1, extractor, this.BAM_TAG_LIST, reportOut, baseQualityThreshold, baseQualityTagName, numBasesBelowQuality);
					writer.addAlignment(r1);
				}
			}
			else { // either you're on the last read, or you're not generating consensus.
				r1= processRead (r1, extractor, this.BAM_TAG_LIST, reportOut, baseQualityThreshold, baseQualityTagName, numBasesBelowQuality);
				writer.addAlignment(r1);
			}
		}

		// close all the inputs/outputs.
		CloserUtil.close(iter);
		CloserUtil.close(writer);
		CloserUtil.close(reportOut);
	}

	private void writeHeader (final PrintStream out, final List<String> bamTags) {
		if (out==null) return;
		String [] line = {"Read_name", "Paired_read_ed","Stopcodon_bc", "PolyA_BC", "R1_Cassette", "R2_Cassette", "ED_R1_Casette_R2_Cassette",
				"ED_GFP_Anchor", "ED_Cassette_Anchor", "Stop_codon_read12_ed", "PolyA_read12_ed"};
		List<String> l = new ArrayList<>(Arrays.asList(line));
		if (bamTags!=null)
			for (String tag: bamTags)
				l.add(tag);
		String h = StringUtils.join(l, "\t");
		out.println(h);
	}

	private void writeBody (final PrintStream out, final SAMRecord r1, final ExtractedSequenceGroup esg, final List<String> bamTags) {
		if (out==null) return;
		List<String> result = new ArrayList<>();
		result.add(r1.getReadName());
		result.add("NA");

		result.add(esg.getStopCodonBarcode().getSequence());
		result.add(esg.getPolyABarcode().getSequence());
		// the read 1 cassette sequence

		String r1Cassette = esg.getStopCodonBarcode().getSequence()+esg.getCassetteAnchor().getSubSequence()+esg.getPolyABarcode().getSequence();
		result.add(r1Cassette);

		result.add("NA");
		result.add("NA");

		int gfpED = esg.getGfpAnchor().getEditDistance().getEditDistance();
		result.add(Integer.toString(gfpED));
		int cassetteED = esg.getCassetteAnchor().getEditDistance().getEditDistance();
		result.add(Integer.toString(cassetteED));
		result.add("NA");
		result.add("NA");

		for (String bamTag: bamTags) {
			Object o = r1.getAttribute(bamTag);
			if (o==null)
				result.add("NA");
			else
				result.add(o.toString());

		}

		String h = StringUtils.join(result, "\t");
		out.println(h);
	}

	private void writeBody (final PrintStream out, final SAMRecord r1, final SAMRecord r2, final ConsensusSequence consensus, final ExtractedSequenceGroup esg, final List<String> bamTags) {
		if (out==null) return;
		List<String> result = new ArrayList<>();
		result.add(r1.getReadName());

		/*
		log.info(r1.getReadName());
		if (r1.getReadName().equals("M00282:312:000000000-AUK2C:1:1101:10001:14455"))
			log.info("STOP");
		*/

		int readED = LevenshteinDistance.getDistance(consensus.getOriginalReadOne(), consensus.getOriginalReadTwo());
		result.add(Integer.toString(readED));
		result.add(esg.getStopCodonBarcode().getSequence());
		result.add(esg.getPolyABarcode().getSequence());
		// get the cassette location (36 bases is standard) from the start of the stop codon BC to the end of the polyA BC.
		// this can be empty if the coordinates are -1 for the start or end.
		int cassetteStart = esg.getStopCodonBarcode().getStart();
		int cassetteEnd = esg.getPolyABarcode().getEnd();
		if (cassetteStart>-1 && cassetteEnd>-1) {
			String r1Cas = consensus.getOriginalSequenceAtConsensusLocation(1, cassetteStart, cassetteEnd);
			String r2Cas = consensus.getOriginalSequenceAtConsensusLocation(2, cassetteStart, cassetteEnd);
			int casED = LevenshteinDistance.getDistance(r1Cas, r2Cas);
			result.add(r1Cas);
			result.add(r2Cas);
			result.add(Integer.toString(casED));
		} else {
			result.add("NA");
			result.add("NA");
			result.add("NA");
		}


		int gfpED = esg.getGfpAnchor().getEditDistance().getEditDistance();
		result.add(Integer.toString(gfpED));
		int cassetteED = esg.getCassetteAnchor().getEditDistance().getEditDistance();
		result.add(Integer.toString(cassetteED));

		int stopCodonStart = esg.getStopCodonBarcode().getStart();
		int stopCodonEnd = esg.getStopCodonBarcode().getEnd();
		if (stopCodonStart>-1 && stopCodonEnd > -1) {
			String stopCodonR1 = consensus.getOriginalSequenceAtConsensusLocation(1, stopCodonStart, stopCodonEnd);
			String stopCodonR2 = consensus.getOriginalSequenceAtConsensusLocation(2, stopCodonStart, stopCodonEnd);
			int stopCodonR1R2ED = LevenshteinDistance.getDistance(stopCodonR1, stopCodonR2);
			result.add(Integer.toString(stopCodonR1R2ED));
		} else
			result.add("NA");

		int polyAStart = esg.getPolyABarcode().getStart();
		int polyAEnd = esg.getPolyABarcode().getEnd();
		if (polyAStart>-1 && polyAEnd > -1) {
			String polyAR1 = consensus.getOriginalSequenceAtConsensusLocation(1, polyAStart, polyAEnd);
			String polyAR2 = consensus.getOriginalSequenceAtConsensusLocation(2, polyAStart, polyAEnd);
			int polyAR1R2ED = LevenshteinDistance.getDistance(polyAR1, polyAR2);
			result.add(Integer.toString(polyAR1R2ED));
		} else
			result.add("NA");

		for (String bamTag: bamTags) {
			Object o = r1.getAttribute(bamTag);
			if (o==null)
				result.add("NA");
			else
				result.add(o.toString());

		}

		String h = StringUtils.join(result, "\t");
		out.println(h);
	}



	/**
	 * Generate consensus sequence on paired records, find barcodes, tag first read of pair and return.
	 * @param r1
	 * @param r2
	 * @param extractor
	 * @return
	 */
	private SAMRecord processReadPair (final SAMRecord r1, final SAMRecord r2, final ExtractBarcodeSequences extractor, final List<String> bamTags, final PrintStream reportOut) {

		ReadPair rp = new ReadPair(r1, r2);
		if (!rp.testProperlyPaired()) {
			log.error(("Reads not properly paired! R1: " + r1.getReadName() + " R2: " + r2.getReadName()));
			System.exit(1);
		}

		ConsensusSequence consensus =  generateAdapterTrimmedConsensus(rp.getFirstRead(), rp.getSecondRead());
		// short circuit if the read(s) are too short.
		if (consensus==null) return r1; // do I need to do some reporting for a too-short read?

		String consensusSequence = consensus.getConsensusSequence();
		r1.setAttribute(this.READ_CONSENSUS_ED_TAG, consensus.getLocalAlignmentEditDistance().getEditDistance());
		// log.info(r1.getReadName());
		ExtractedSequenceGroup result = extractor.findRabiesBarcode(consensusSequence);
		SAMRecord rec = setTagsOnBAM(r1, result);
		if (reportOut !=null) writeBody(reportOut, r1, r2, consensus, result, bamTags);
		return rec;
	}

	private ConsensusSequence generateAdapterTrimmedConsensus (final SAMRecord r1, final SAMRecord r2) {
		String trimmedReadOne = getTrimmedReadSequence(r1);
		String trimmedReadTwo = getTrimmedReadSequence(r2);
		if (trimmedReadOne.length()<this.MIN_READ_LENGTH || trimmedReadTwo.length() < this.MIN_READ_LENGTH) return null;
		ConsensusSequence consensus =  ConsensusSequenceFactory.getInstance().getConsensusSequence(trimmedReadOne, trimmedReadTwo, ASSUME_REVERSE_COMPLIMENT);
		consensus.addReadBaseQualities(getTrimmedReadQuality(r1), getTrimmedReadQuality(r2));
		return consensus;
	}


	private String getTrimmedReadSequence (final SAMRecord r) {
		String seq = r.getReadString();
		Integer clipFrom = r.getIntegerAttribute(this.ADAPTER_TAG);
		if (clipFrom==null) return seq;
		seq=seq.substring(0, clipFrom-1);
		return seq;
	}

	private String getTrimmedReadQuality (final SAMRecord r) {
		String qual = r.getBaseQualityString();
		Integer clipFrom = r.getIntegerAttribute(this.ADAPTER_TAG);
		if (clipFrom==null) return qual;
		qual=qual.substring(0, clipFrom-1);
		return qual;
	}



	/**
	 * Find barcodes on read, tag read, return.
	 * @param r1
	 * @param extractor
	 * @return
	 */
	private SAMRecord processRead (final SAMRecord r1, final ExtractBarcodeSequences extractor, final List<String> bamTags, final PrintStream reportOut, final int baseQualityThrehsold, final String tagName, final int numBasesBelowQuality) {
		String sequence = getTrimmedReadSequence(r1);
		// skip tagging because read is too short!
		if (sequence.length()<this.MIN_READ_LENGTH) return r1;
		ExtractedSequenceGroup result = extractor.findRabiesBarcode(sequence);
		SAMRecord rec =setTagsOnBAM(r1, result);
		int score = getBaseQualityExtractedSequences(r1, result, baseQualityThrehsold);
		failedBaseHistogram.increment(score);
		if (score>=numBasesBelowQuality) r1.setAttribute(tagName, score);
		if (reportOut !=null) writeBody(reportOut, r1, result, bamTags);
		return rec;
	}

	private int getBaseQualityExtractedSequences (final SAMRecord r1, final ExtractedSequenceGroup result, final int baseQualityThrehsold) {
		List<BaseRange> baseRanges = new ArrayList<>();
		// baseRanges.add(new BaseRange(result.getCassetteAnchor().getStart(), result.getCassetteAnchor().getEnd()));
		// baseRanges.add(new BaseRange(result.getGfpAnchor().getStart(), result.getGfpAnchor().getEnd()));
		if (result.getPolyABarcode().isValid())
			baseRanges.add(new BaseRange (result.getPolyABarcode().getStart(), result.getPolyABarcode().getEnd()));
		if (result.getStopCodonBarcode().isValid())
			baseRanges.add(new BaseRange (result.getStopCodonBarcode().getStart(), result.getStopCodonBarcode().getEnd()));
		BaseQualityFilter f = new BaseQualityFilter(baseRanges, baseQualityThrehsold);
		int score = f.scoreBaseQuality(r1);
		return score;

	}

	/**
	 * Given an extracted rabies barcode group and a read, transfer the information from the group to the read.
	 * @param r1
	 * @param esg
	 * @return
	 */
	private SAMRecord setTagsOnBAM (final SAMRecord r1, final ExtractedSequenceGroup esg) {
		r1.setAttribute(this.CASSETTE_AHCHOR_ED_TAG, esg.getCassetteAnchor().getEditDistance().getEditDistance());
		r1.setAttribute(this.GFP_ANCHOR_ED_TAG, esg.getGfpAnchor().getEditDistance().getEditDistance());
		r1.setAttribute(this.STOP_CODON_BARCODE_TAG, esg.getStopCodonBarcode().getSequence());
		r1.setAttribute(this.POLY_A_BARCODE_TAG, esg.getPolyABarcode().getSequence());
		r1.setAttribute(this.RABIES_BARCODE, esg.getStopCodonBarcode().getSequence()+esg.getPolyABarcode().getSequence());
		return r1;
	}

	/**
	 * Remove all the paired flag info from a read.
	 * @param read
	 * @return
	 */
	/*
	private SAMRecord makeReadUnpaired (final SAMRecord read) {
		int flag =read.getFlags();
		if (read.getMateUnmappedFlag()) flag-=8;
		if (read.getMateNegativeStrandFlag()) flag-=32;
		if (read.getReadPairedFlag()) flag-=1;
		if (read.getFirstOfPairFlag()) flag-=64;
		if (read.getSecondOfPairFlag()) flag-=128;
		read.setFlags(flag);
		return read;
	}


	private class TagResultPaired {
		private final SAMRecord r1;
		private final SAMRecord r2;
		private final ConsensusSequence consensus;
		private final ExtractedSequenceGroup g;

		public TagResultPaired (final SAMRecord r1, final SAMRecord r2, final ConsensusSequence consensus, final ExtractedSequenceGroup g) {
			this.r1=r1;
			this.r2=r2;
			this.consensus=consensus;
			this.g= g;
		}

		public SAMRecord getRead1() {
			return r1;
		}

		public SAMRecord getRead2() {
			return r2;
		}

		public ConsensusSequence getConsensus() {
			return consensus;
		}

		public ExtractedSequenceGroup getExtractedSequenceGroup() {
			return g;
		}
	}
	/**
	 * Open up the BAM briefly, test if the first two reads have the same name.
	 * @param input
	 * @return
	 */
	/*
	private boolean testReadsPaired (final File input) {
		SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
		Iterator<SAMRecord> iter = reader.iterator();
		String rn1="";
		String rn2="";
		if (iter.hasNext())
			rn1= iter.next().getReadName();
		if (iter.hasNext())
			rn2= iter.next().getReadName();
		CloserUtil.close(rn2);
		return (rn1.equals(rn2));
	}
	*/
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new TagReadWithRabiesBarcodes().instanceMain(args));
	}
}
