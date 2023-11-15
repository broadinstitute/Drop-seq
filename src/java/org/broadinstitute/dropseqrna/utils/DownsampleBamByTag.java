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

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.metrics.BamTagHistogram;
import org.broadinstitute.dropseqrna.utils.readiterators.MissingTagFilteringIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;

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
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.BasicInputParser;

@CommandLineProgramProperties(summary = "",
        oneLineSummary = "",
        programGroup = DropSeq.class)
public class DownsampleBamByTag extends CommandLineProgram {

	private final Log log = Log.getInstance(DownsampleBamByTag.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<File> INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BAM.")
	public File OUTPUT;

	@Argument(doc="The tag to organize reads around for downsampling")
	public String TAG="ZC";

	@Argument(doc="The map quality of the read to be included.")
	public int READ_MQ=10;

	@Argument(doc="The number of reads to include per TAG.", mutex={"TAG_FILE"})
	public Integer NUM_READS;

	@Argument(doc="Filter PCR Duplicates.  Defaults to false")
	public boolean FILTER_PCR_DUPLICATES=false;

	@Argument(doc="Instead of specifying a set number of reads to downsample all tags by, instead input a file that has values for the tag, with numbers of reads to downsample by.  If this is used then the NUM_READS is ignored." +
			"Has 2 columns.  The value of the tag, and the number of reads.  Tab seperated.", mutex={"NUM_READS"})
	public File TAG_FILE;

	@Argument (doc="If true, use a probabilistic strategy that will not generate the exact number of reads asked for (though very close), but will run much faster.")
	public Boolean USE_PROBABILISTIC_STRATEGY=true;

	@Argument (doc="If true, assume reads are paired, and maintain reads as pairs.  If unpaired reads are detected, they will be discarded. If paired reads that do"
			+ "not have matching mate information are detected, they are excluded from downsampling."
			+ "This option is more expensive to run, as the data needs to be sorted first.  Only available in probabilistic mode")
	public Boolean PAIRED_READS=false;

    @Argument(doc = "Random seed to use if deterministic behavior is desired.  " +
            "Setting to null will cause multiple invocations to produce different results.")
    public Integer RANDOM_SEED = 1;

	private int numUnpairedReads=0;
	private int invalidPairs=0;
	@Override
	protected int doWork() {
		if (PAIRED_READS && !USE_PROBABILISTIC_STRATEGY)
			log.error("PAIRED READS mode only available if USE_PROBABILISTIC_STRATEGY=true");

		ObjectCounter<String> tc = null;

		if (TAG_FILE!=null) {
			tc=parseTagFile(TAG_FILE);
		}

		if (USE_PROBABILISTIC_STRATEGY)
			downsampleBAMByTagProbabilistic (this.INPUT, this.TAG, this.READ_MQ, this.OUTPUT, tc, this.PAIRED_READS, this.NUM_READS, this.RANDOM_SEED);
		else
			downsampleBAMByTag(this.INPUT, this.TAG, this.READ_MQ, this.OUTPUT, tc, this.NUM_READS, new Random (this.RANDOM_SEED));
		return (0);
	}
	
	@Override
	protected String[] customCommandLineValidation() {

		final ArrayList<String> list = new ArrayList<>(1);
		
		this.INPUT = FileListParsingUtils.expandFileList(INPUT);
				
		IOUtil.assertFileIsWritable(OUTPUT);
		
		if (TAG_FILE!=null) 
			IOUtil.assertFileIsReadable(TAG_FILE);				
		
		return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
	}

	/**
	 * Run downsampling, but do it in a probabilistic way instead of exact.
	 * This requires fewer passes/sorting of the data, so while it's slightly less perfect, it's a lot faster.
	 * @param input
	 * @param tag
	 * @param readQuality
	 * @param output
	 * @param desiredReads
	 * @param numReads
	 */
	public void downsampleBAMByTagProbabilistic (final List<File> input, final String tag, final int readQuality, final File output, ObjectCounter<String> desiredReads, final boolean pairedReadMode, final Integer numReads, final int randomSeed) {
		Random random = new Random (randomSeed);

		// get the count of the number of reads for each tag.
		BamTagHistogram bth = new BamTagHistogram();
		ObjectCounter<String> observedReads = bth.getBamTagCounts (input, tag, readQuality, this.FILTER_PCR_DUPLICATES);
		// if there's no input desired reads per tag, generate a desired number of reads per tag set to size numReads. 
		if (desiredReads==null) desiredReads=getDesiredReads(observedReads, numReads);
		
		PeekableIterator<SAMRecord> iter=null;

		SamHeaderAndIterator headerAndIter = SamFileMergeUtil.mergeInputs(input, false, SamReaderFactory.makeDefault());
		
		if (pairedReadMode)
			iter=new PeekableIterator<>(CustomBAMIterators.getQuerynameSortedRecords(headerAndIter));
		else
			iter=new PeekableIterator<>(headerAndIter.iterator);

		SAMFileHeader h= headerAndIter.header;
		// when constructing the writer, the reads will be out of order in paired read mode since we're query name ordering the reads and tag sorting them.
		SAMFileWriter writer=null;
		if (pairedReadMode) writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(h, false, output);
		else writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(h, true, output);

		ProgressLogger pl = new ProgressLogger(this.log);

		Map<String, Double> tagProbabilityMap = getTagProbabilities(desiredReads, observedReads);

		// do all the filtering of reads and downsampling on this pass.
		while (iter.hasNext()) {
			if (!pairedReadMode) processSingleProbabilistic (iter, tagProbabilityMap, tag, writer, pl, random);
			if (pairedReadMode)  processDoubleProbabilistic (iter, tagProbabilityMap, tag, writer, pl, random);
		}
		CloserUtil.close(iter);
		
		if (pairedReadMode) log.info("Done adding downsampled reads, putting data back into original sort order.");
		writer.close();
		if (pairedReadMode) {
			log.info("Found [" + this.numUnpairedReads + "] unpaired reads that were rejected.");
			log.info("Found [" + this.invalidPairs + "] improper read pairs that were rejected.");
		}
	}

	private ObjectCounter<String> getDesiredReads (ObjectCounter<String> observedReads, Integer numReads ) {
		ObjectCounter<String> result = new ObjectCounter<String>();
		for (String k: observedReads.getKeys()) {
			result.incrementByCount(k, numReads);
		}
		return result;
	}
	
	private void processSingleProbabilistic (final PeekableIterator<SAMRecord> iter, final Map<String, Double> tagProbabilityMap, final String tag, final SAMFileWriter writer, final ProgressLogger pl, final Random random) {
		SAMRecord r = iter.next();
		pl.record(r);
		if (!filterRecord(r)) {
			String recordTag = r.getStringAttribute(tag);
			// ignore the read if you don't have a tag with a value, or if the tag isn't set on the read.
			if (recordTag==null) return;
			if (!tagProbabilityMap.containsKey(recordTag)) return;
			double expectedRate = tagProbabilityMap.get(recordTag);
			double randomNumber = random.nextDouble();
			// if below the probability, add the read.
			if (randomNumber<=expectedRate)
				writer.addAlignment(r);
		}
	}

	/**
	 * Like processSingleProbabilistic, but looks at reads in pairs, and the iterator is queryname sorted.
	 * @param iter
	 * @param tagProbabilityMap
	 * @param pl
	 */
	private void processDoubleProbabilistic (final PeekableIterator<SAMRecord> iter, final Map<String, Double> tagProbabilityMap, final String tag, final SAMFileWriter writer, final ProgressLogger pl, final Random random) {
		// get the first read, and gather all the reads that belong with this one.
		SAMRecord r = iter.next();
		pl.record(r);
		List<SAMRecord> reads = new ArrayList<>();
		reads.add(r);

		while (true) {
			SAMRecord r2 = iter.peek();
			if (r2==null)
				break;
			if (r.getReadName().equals(r2.getReadName())) {
				r2=iter.next();
				reads.add(r2);
			} else
				break;
		}
		// read is unpaired.
		if (reads.size()==1) {
			numUnpairedReads++; // an unpaired read.
			return;
		}

		List<SAMRecord> filteredReads =reads.stream().filter(x -> !filterRecord(x)).collect(Collectors.toList());

		// if there aren't at least 2 reads that pass filters, then reject both.
		if (filteredReads.size()!=2)
			return;

		// you made it this far, there are only 2 reads.
		r=filteredReads.get(0);
		SAMRecord r2= filteredReads.get(1);

		// validate that the reads are paired by comparing all the mate info.
		boolean validPair = validateMatePairingInfo(r, r2);
		if (!validPair) {
			invalidPairs++;
			if (invalidPairs<100)
				log.info("Invalid read pair: " + r.toString() + r2.toString());
			return;
		}

		String recordTag = r.getStringAttribute(tag);
		String recordTag2 = r2.getStringAttribute(tag);

		// ignore the read if you don't have a tag with a value, or if the tag isn't set on the read.
		if (recordTag==null || recordTag2==null) return;

		// the tags don't match on the read pair, this should never happen...
		if (!recordTag.equals(recordTag2))
			throw new IllegalStateException("Tags don't match on paired read!");

		if (!tagProbabilityMap.containsKey(recordTag)) return;
		double expectedRate = tagProbabilityMap.get(recordTag);
		double randomNumber = random.nextDouble();
		// if below the probability, add the read.
		if (randomNumber<=expectedRate) {
			writer.addAlignment(r);
			writer.addAlignment(r2);
		}
	}

	private boolean validateMatePairingInfo (final SAMRecord r1, final SAMRecord r2) {
		if (r1.getAlignmentStart()!=r2.getMateAlignmentStart()) return false;
		if (r1.getReadNegativeStrandFlag() !=r2.getMateNegativeStrandFlag()) return false;
		if (r1.getReferenceIndex()!=r2.getMateReferenceIndex()) return false;
		if (r1.getReadUnmappedFlag()!= r2.getMateUnmappedFlag()) return false;
		return true;
	}

	/**
	 *
	 * @param desiredReads The counts of the number of reads desired for each tag
	 * @param observedReads The number of reads on each tag that exist in the BAM
	 * @return The probability that we need to downsample each tag at.
	 */
	private Map<String, Double> getTagProbabilities (final ObjectCounter<String> desiredReads, final ObjectCounter<String> observedReads) {
		Map<String, Double> result = new HashMap<>(desiredReads.getKeys().size());

		for (String key: desiredReads.getKeys()) {
			int obs = observedReads.getCountForKey(key);
			int d = desiredReads.getCountForKey(key);
			double prob = (double) d / (double) obs;
			if (prob>1) prob=1d;
			result.put(key, prob);
		}
		return (result);
	}

	public void downsampleBAMByTag (final List<File> input, final String tag, final int readQuality, final File output, final ObjectCounter<String> desiredReads, final ObjectCounter<String> observedReads, final Integer numReads, final Random random) {
		SamHeaderAndIterator headerAndIter = SamFileMergeUtil.mergeInputs(input, false, SamReaderFactory.makeDefault());

		SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(headerAndIter.header, false, output);

		ProgressLogger pl = new ProgressLogger(this.log);

        final Iterator<SAMRecord> filteringIterator = new MissingTagFilteringIterator(headerAndIter.iterator, tag);
        final PeekableIterator<SAMRecord> toi =
                new PeekableIterator<>(SamRecordSortingIteratorFactory.create(headerAndIter.header, filteringIterator, new StringTagComparator(tag), pl));

		SAMRecord r = toi.peek();

		//set up logging variables
		int tagNumber=0;
		int totalTags=observedReads.getSize();

		// set up tracking variables
		String currentTag = "";

		Integer numExpected = 0;
		int recsSeen=0;
		Set<Integer> index = null;
		pl = new ProgressLogger(this.log);

		log.info("Writing downsampled BAM");
		while (toi.hasNext()) {
			r = toi.peek();
			// if the record is junk, move on.
			if (filterRecord(r)) {
				r=toi.next();
				pl.record(r);
				continue;
			}
			// get the TAG for this record.
			String recordTag = r.getStringAttribute(tag);

			// if the tag doesn't match the current tag, swap onto the new tag, start the loop over.
			if (!recordTag.equals(currentTag)) {
				tagNumber++;
				currentTag = r.getStringAttribute(tag);
				numExpected = observedReads.getCountForKey(currentTag);
				recsSeen=0;
				int numReadsToSample=getNumRecordsToSample(currentTag, observedReads, desiredReads, numReads);
				index = getIndex (numExpected, numReadsToSample, random);
				if (tagNumber%10000==0) log.info("Processing tag " + currentTag + " [" + tagNumber + "] of [" + totalTags +"]");
				continue;
			}

			// you're a real record now.
			r=toi.next();
			pl.record(r);
			// if the read is in the index selected to use, write out the alignment.
			if (index.contains(Integer.valueOf(recsSeen)))
				writer.addAlignment(r);

			recsSeen++;

		}
		CloserUtil.close(toi);
		writer.close();

	}


	public void downsampleBAMByTag (final List<File> input, final String tag, final int readQuality, final File output, final ObjectCounter<String> tc, final Integer numReads, final Random random) {
		BamTagHistogram bth = new BamTagHistogram();
		ObjectCounter<String> histogram = bth.getBamTagCounts (input, tag, readQuality, this.FILTER_PCR_DUPLICATES);
		downsampleBAMByTag(input, tag, readQuality, output, tc, histogram, numReads, random);
	}

	/**
	 *
	 * For a barcode, figure out how many reads exist, and how many reads should be selected.
	 * If the desiredCount is populated, this CAN supply the number of reads if the barcode exists in the object.
	 * Otherwise, use the default count (which can be null), or the expected count.
	 * @param barcode The barcode to select reads from.
	 * @param expected The number of reads for this tag in the BAM.
	 * @param desiredCount The number of reads to downsample to if the TAG_FILE option is invoked.  Not every barocde will be included here, so results from this object for a barcode can be null.
	 * @param defaultCount
	 * @return
	 */
	private int getNumRecordsToSample (final String barcode, final ObjectCounter<String> expected, final ObjectCounter<String> desiredCount, final Integer defaultCount) {
		int result = 0;
		if (desiredCount!=null && desiredCount.hasKey(barcode))
			result = desiredCount.getCountForKey(barcode);
		else if (defaultCount!=null)
			result=defaultCount;
		return (result);
	}

	boolean filterRecord (final SAMRecord r) {
		if (this.FILTER_PCR_DUPLICATES && r.getDuplicateReadFlag()) return true;
		if (r.getMappingQuality()<this.READ_MQ) return true;
		if (r.isSecondaryOrSupplementary()) return true;
		return false;
	}

	/**
	 * My best guess to downsampling without replacement is figure out how many reads there are in a tag,
	 * then make an index that's a sample of the full set.
	 *
	 * @param totalNumReads
	 * @return
	 */
	// TODO: change this to an object wrapping a bitset of length total num reads.  Bitset would have randomly turned on positions if inclusion rate was below 50%, or randomly turned off if above 50%, to minimize the number of collisions.
	// Should be far more memory efficient than holding a list of Integers of <totalNumReads> in memory.
	private Set<Integer> getIndex (final int totalNumReads, final int numToSelect, final Random random) {
		Set<Integer> result = new HashSet<>();
		if (numToSelect==0) return result;

		List<Integer> temp = new ArrayList<>();
		for (int i=0; i<totalNumReads; i++)
			temp.add(i);
		// if you can't downsample, return the entire set.
		if (totalNumReads<=numToSelect)
			return new HashSet<>(temp);

		// otherwise, select out a slice.
		Collections.shuffle(temp, random);

		for (int i=0; i<numToSelect; i++)
			result.add(temp.get(i));
		return (result);
	}

	public static ObjectCounter<String> parseTagFile(final File f) {
		ObjectCounter<String> result = new ObjectCounter<>();
		BasicInputParser parser = new BasicInputParser(false, 2, f);
		while(parser.hasNext()) {
			String [] line =parser.next();
			String tagValue = line[0];
			int count = Integer.parseInt(line[1]);
			result.setCount(tagValue, count);

		}
		parser.close();
		return (result);
	}


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new DownsampleBamByTag().instanceMain(args));
	}


}
