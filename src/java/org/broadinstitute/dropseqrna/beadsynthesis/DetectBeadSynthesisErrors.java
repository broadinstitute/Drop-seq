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
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.barnyard.BarcodeListRetrieval;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.BaseDistributionMetric;
import org.broadinstitute.dropseqrna.utils.BaseDistributionMetricCollection;
import org.broadinstitute.dropseqrna.utils.Bases;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.utils.editdistance.MapBarcodesByEditDistance;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

/**
 *
 * @author nemesh
 *
 */

@CommandLineProgramProperties(
        summary = "For each cell, gather up all the UMIs.  An error in synthesis will result in the last base of the synthesis being fixed in >90% of the UMIs for that cell, across all genes." +
			"This fixed base is T.  For cell barcodes where this occurs, output the cell barcode in a file, as well as (optionally) pad the cell barcodes with N for the error bases.  In cases where we don't know how to fix"
			+ "the error - when there are too many missing bases in the synthesis, or the synthesis error isn't one of the repairable types, we remove the read from the output BAM.",
        oneLineSummary = "Detect barcode synthesis errors where the final base of a UMI is fixed across all UMIs of a cell.",
        programGroup = DropSeq.class
)

public class DetectBeadSynthesisErrors extends CommandLineProgram  {

	private static final Log log = Log.getInstance(DetectBeadSynthesisErrors.class);

	@Argument(doc="Output of detailed information on each cell barcode analyzed.  Each row is a single cell barcode.  "
			+ "The data has multiple columns: the cell barcode, the number of UMIs, then one column per UMI base position containing the count of the reads, with a | "
			+ "delimiter between bases.  Bases are ordered A,C,G,T for these columns.  An example output with a single base UMI would be:"
			+ "AAAAAA	20		5|4|6|5.")
	public File OUTPUT_STATS;

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM files to analyze.  They must all have the same sort order")
	public List<File> INPUT;

	@Argument(doc="Output a summary of the error types and frequencies detected")
	public File SUMMARY;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output BAM, with the synthesis error barcodes removed", optional=true)
	public File OUTPUT;

	@Argument(doc="A report of which barcodes where collapsed, what the intended sequences were, and what base positions were changed.  There are 2 modes of errors reported. "
			+ "The first is an insertion/deletion changes at bases 1-11 in the intended sequence, where a base in the intended sequence is partially incorperated into the related sequence "
			+ "resulting in a deletion in the intended sequence generating a new sequence, and an 'insertion' at the last base of the sequence, which is the first base of the UMI. "
			+ "The second is a substitution event at base 12, where the last base is partially incorperated.  Both changes are the same mechanism, but the signature of the change [indel vs subsitution] differs."
			+ "Not all repaired barcodes will have an intended sequence - if the base is incorperated at a very low rate compared to the number of UMIs in the repaired cell barcode library size, "
			+ "the intended sequence may be too small to find.", optional=true)
	public File REPORT=null;

	@Argument(doc="The sequence of the primer.", optional=true)
	public String PRIMER_SEQUENCE=null;

	@Argument(doc="When looking at fixed UMIs, see if the edit distance from the UMI to the primer is within this threshold.  0 indicates a perfect match between the primer and the UMI.")
	public Integer EDIT_DISTANCE=0;

	@Argument(doc="The cell barcode tag.")
	public String CELL_BARCODE_TAG="XC";

	@Argument(doc="The molecular barcode tag.")
	public String MOLECULAR_BARCODE_TAG="XM";

	@Argument(doc="The Gene/Exon tag")
	public String GENE_EXON_TAG="GE";

	@Argument(doc="The strand of the gene(s) the read overlaps.  When there are multiple genes, they will be comma-separated.")
	public String STRAND_TAG="GS";

	@Argument(doc="The map quality of the read to be included when calculating the barcodes in <NUM_BARCODES>")
	public Integer READ_MQ=10;

	// HAS TO BE FILLED IN.
	@Argument (doc="The minimum number of UMIs required to report a cell barcode.  If this parameter is specified and [NUM_BARCODES or CELL_BC_FILE] are not specified, then all cell barcodes with at least this many UMIs will be analyzed., optional=false")
	public Integer MIN_UMIS_PER_CELL=20;

	// OPTIONAL
	@Argument (doc="Find the top set of <NUM_BARCODES> most common barcodes by HQ reads and only use this set for analysis.", optional=true)
	public Integer NUM_BARCODES;

	// OPTIONAL
	@Argument(doc="Override NUM_BARCODES, and process reads that have the cell barcodes in this file instead.  The file has 1 column with no header.", optional=true)
	public File CELL_BC_FILE;

	// @Argument(doc="Repair Synthesis errors with at most this many missing bases detected.", optional=true)
	private Integer MAX_NUM_ERRORS=1;

	@Argument(doc="Number of threads to use for edit distance collapse.  Defaults to 1.")
	public int NUM_THREADS=1;

	Double EXTREME_BASE_RATIO=0.8;
	DetectPrimerInUMI detectPrimerTool=null;

	SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE);

	private Character PAD_CHARACTER='N';
	private static DecimalFormat df2 = new DecimalFormat("#.##");

	@Override
	protected int doWork() {
		// primer detection if requested.
		if (this.PRIMER_SEQUENCE!=null)
			this.detectPrimerTool = new DetectPrimerInUMI(this.PRIMER_SEQUENCE);

		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT_STATS));

		UMIIterator iterator = prepareUMIIterator();
		BiasedBarcodeCollection biasedBarcodeCollection = findBiasedBarcodes(iterator, out, this.SUMMARY);
		Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions = biasedBarcodeCollection.getBiasedBarcodes();

		ObjectCounter<String> umisPerCell = biasedBarcodeCollection.getUMICounts();
		Map<String, Double> umiBias = biasedBarcodeCollection.getUMIBias();

		// Group T-biased barcodes into groups of neighbors.  The key is the padded sequence the neighbors share
		Map<String, BarcodeNeighborGroup> barcodeNeighborGroups = buildBarcodeNeighborGroups(errorBarcodesWithPositions.values(), this.EXTREME_BASE_RATIO);

		// find intended sequences for these neighbors
		Map<String, String> intendedSequenceMap = findIntendedBarcodes(barcodeNeighborGroups, umisPerCell.getKeys(), umiBias);

		// write report about which sequences were collapsed and why
		writeReport(umisPerCell, barcodeNeighborGroups.values(), intendedSequenceMap, umiBias, this.REPORT);


		// clean up the BAM if desired.
		if (this.OUTPUT!=null)
			cleanBAM(errorBarcodesWithPositions, intendedSequenceMap);
		return 0;
	}

	void writeReport (final ObjectCounter<String> umisPerCell, final Collection<BarcodeNeighborGroup> neighborGroups, final Map<String, String> intendedSequenceMap, final Map<String, Double> umiBias, final File outFile) {
		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outFile));

		IntendedSequenceBuilder b = new IntendedSequenceBuilder(umisPerCell, umiBias);

		// write header.
		String [] header = {"intended_sequence", "related_sequences", "deleted_base", "deleted_base_pos", "non_incorperation_rate", "intended_UMIs", "related_median_UMIs", "intended_TBias", "related_median_TBias"};
		out.println(StringUtils.join(header, "\t"));

		// for each neighbor group, write out a report line.
		for (BarcodeNeighborGroup ng: neighborGroups) {
			// this is slightly janky, but the intended sequence map holds raw neighbor barcodes to the intended sequence, so you only need the first barcode in the list.
			String intendedSequence = intendedSequenceMap.get(ng.getNeighborCellBarcodes().iterator().next());
			IntendedSequence is = b.build(intendedSequence, ng);
			writeLine(is, out);
		}
		CloserUtil.close(out);
	}

	void writeLine (final IntendedSequence is, final PrintStream out) {
		List<String> line = new ArrayList<>();
		line.add(getNullableString(is.getIntendedSequence()));
		line.add(StringUtils.join(is.getRelatedSequences(), ":"));
		line.add(getNullableString(is.getDeletedBase()));
		line.add(getNullableString(is.getDeletedBasePos()));

		if (is.getDeletionRate()==null)
			line.add("NA");
		else
			line.add(df2.format(is.getDeletionRate()));

		if (is.getIntendedSequenceUMIs()==null)
			line.add("NA");
		else
			line.add(is.getIntendedSequenceUMIs().toString());

		line.add(df2.format(is.getMedianRelatedSequenceUMIs()));

		if (is.getIntendedSequenceUMIBias()==null)
			line.add("NA");
		else
			line.add(df2.format(is.getIntendedSequenceUMIBias()));
		line.add(df2.format(is.getRelatedMedianUMIBias()));

		String l = StringUtils.join(line, "\t");
		out.println(l);

	}

	private String getNullableString (final Object o) {
		if (o==null) return "NA";
		return o.toString();
	}


	/**
	 * For each BarcodeNeighborGroup, try to find a sequence that gave rise to these sequences.  This will not always be successful!
	 * This provides a map from the biased barcodes to the intended sequences where it can be found, but does not include all biased barcodes.
	 * @return a map from the biased barcode to the intended sequence.  There will be multiple biased barcodes linked to a single intended sequence.
	 */
	private Map<String, String> findIntendedBarcodes (final Map<String, BarcodeNeighborGroup> biasedGroups, final Collection<String> collection, final Map<String, Double> umiBias) {

		Set<String> potentialIntendedBarcodes = new HashSet<>(collection);
		// gather up biased barcocdes
		Set<String> toRemove = new HashSet<>();
		for (BarcodeNeighborGroup bng: biasedGroups.values())
			toRemove.addAll(bng.getNeighborCellBarcodes());
		// remove biased barcodes from potential intended barcodes.
		potentialIntendedBarcodes.removeAll(toRemove);

		List<String> potentialIntendedBarcodesList = new ArrayList<>(potentialIntendedBarcodes);

		List<String> repairedCellBarcodes = new ArrayList<>(biasedGroups.keySet());

		MapBarcodesByEditDistance mbed = new MapBarcodesByEditDistance(true, this.NUM_THREADS);
		Map<String, String> intendedSequenceMap = mbed.findIntendedIndelSequences (repairedCellBarcodes, potentialIntendedBarcodesList, 1);

		Map<String,String> result = new HashMap<>();
		// for each intended sequence, add each of the neighbors to the map.
		for (String paddedSequence: intendedSequenceMap.keySet()) {
			BarcodeNeighborGroup bng = biasedGroups.get(paddedSequence);
			String intendedSeq = intendedSequenceMap.get(paddedSequence);
			if (intendedSeq!=null) {
				bng.setIntendedSequence(intendedSeq);
				// map all neighbors to the intended sequence.
				bng.getNeighborCellBarcodes().stream().forEach(x-> result.put(x, intendedSeq));
			}
		}

		//TODO: are we missing BarcodeNeighborGroup?  Probably.  which ones?
		// ones where no there's no intended sequence.
		/*
		Set<String> unrepairedBarcodes = new HashSet<> (biasedGroups.keySet());
		unrepairedBarcodes.removeAll(intendedSequenceMap.keySet());
		for (String paddedSequence: unrepairedBarcodes) {
			log.info("STOP");
			BarcodeNeighborGroup ng = biasedGroups.get(paddedSequence);
			// if there's big bias in all but one of the barcodes of the BarcodeNeighborGroup we could call that the intended sequence.
			for (String n: ng.getNeighborCellBarcodes()) {
				double umiBiasN = umiBias.get(n);
			}
		}
		*/
		return result;
	}




	/**
	 * Find all the cell barcodes that are biased.
	 *
	 * This walks through many (perhaps all?) of the cell barcodes, and for cell barcodes that have a sufficient number of UMIs, looks to see if there's an error
	 * @param iter The BAM iterator to walk through
	 * @param out the verbose output stream.
	 * @param outSummary The summary output stream.
	 *
	 * TODO: A less memory-hog version of this would write out the summary file as it runs.  Could even write this out to a SortingIteratorFactory by implementing a codec...
	 * Only need to hang onto errors that are SYNTH_MISSING_BASE, and leave the rest null (and check for that when running barcode repair.)
	 *
	 * @return A collection of biased cell barcodes.
	 */
	private BiasedBarcodeCollection findBiasedBarcodes (final UMIIterator iter, final PrintStream out, final File outSummary) {
		log.info("Finding Cell Barcodes with UMI errors");
		// Group the stream of UMICollections into groups with the same cell barcode.
        GroupingIterator<UMICollection> groupingIterator = new GroupingIterator<>(iter,
                new Comparator<UMICollection>() {
                    @Override
                    public int compare(final UMICollection o1, final UMICollection o2) {
                        return o1.getCellBarcode().compareTo(o2.getCellBarcode());
                    }
                });


		// for holding barcodes results.  The key is the cell barcode, the value is the first base to pad.
		// Used for cleanup of BAMs.
		Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions = new HashMap<>();

 		// for holding UMI Strings efficiently
		// TODO: evaluate if this is needed this anymore now that the report is written out disk immediately.
		StringInterner  umiStringCache = new StringInterner();

		// a sorting collection so big data can spill to disk before it's sorted and written out as a report.
		SortingCollection<BeadSynthesisErrorData> sortingCollection= SortingCollection.newInstance(BeadSynthesisErrorData.class, new BeadSynthesisErrorDataCodec(), new BeadSynthesisErrorData.SizeComparator(), this.MAX_RECORDS_IN_RAM);

        // gather up summary stats
     	BeadSynthesisErrorsSummaryMetric summary = new BeadSynthesisErrorsSummaryMetric();

     	// track the list of cell barcodes with sufficient UMIs to process.  We'll need them later to find intended sequences.
     	ObjectCounter<String> umisPerCellBarcode = new ObjectCounter<>();

     	// track the UMI Bias at the last base
     	Map<String, Double> umiBias = new HashMap<>();

     	// a log for processing
     	ProgressLogger prog = new ProgressLogger(log, 1000000, "Processed Cell/Gene UMIs");

     	// main data generation loop.
     	// to ease memory usage, after generating the BeadSynthesisErrorData object, use its cell barcode string for registering additional data.
        for (final List<UMICollection> umiCollectionList : groupingIterator) {
            BeadSynthesisErrorData bsed = buildBeadSynthesisErrorData(umiCollectionList, umiStringCache, prog);
            // if the cell has too few UMIs, then go to the next cell and skip all processing.
            if (bsed.getUMICount() < this.MIN_UMIS_PER_CELL)
				// not sure I even want to track this...
            	// summary.LOW_UMI_COUNT++;
            	continue;

            // add the result to the summary
            summary=addDataToSummary(bsed, summary);
            umisPerCellBarcode.incrementByCount(bsed.getCellBarcode(), bsed.getNumTranscripts()); // track the cell barcode if it's sufficiently large to process.

            // explicitly call getting the error type.
            BeadSynthesisErrorType errorType=bsed.getErrorType(this.EXTREME_BASE_RATIO, this.detectPrimerTool, this.EDIT_DISTANCE);
            // gather up the UMI bias at the last base.

            double barcodeUMIBias = bsed.getPolyTFrequencyLastBase();
            umiBias.put(bsed.getCellBarcode(), barcodeUMIBias);

            // finalize object so it uses less memory.
            bsed.finalize();
            // only add to the collection if you have UMIs and a repairable error.
            if (bsed.getUMICount()>=this.MIN_UMIS_PER_CELL && errorType==BeadSynthesisErrorType.SYNTH_MISSING_BASE)
            	errorBarcodesWithPositions.put(bsed.getCellBarcode(), bsed);

            // add to sorting collection if you have enough UMIs.
            sortingCollection.add(bsed);
        }

        PeekableIterator<BeadSynthesisErrorData> bsedIter = new PeekableIterator<>(sortingCollection.iterator());

        log.info("Writing Biased UMI reports");
        // write out the records from the sorting collection.
        writeFile(bsedIter, out);
        // write out the summary
        writeSummary(summary, outSummary);
        CloserUtil.close(bsedIter);

        // the error barcodes we want to fix.
        BiasedBarcodeCollection result = new BiasedBarcodeCollection(errorBarcodesWithPositions, umisPerCellBarcode, umiBias);
        return result;
	}

	/**
	 * For a single cell barcode, gather up all the reads/UMIs to test for UMI errors.
	 * @param umiCollectionList A collection of UMIs for a cell.
	 * @param umiStringCache A cache to save space on UMIs - there are only 4^length possible UMIs, and length is usually 8-9, so UMIs are often repeated.
	 * @param prog A progress logger.
	 * @return A BeadSynthesisErrorData object for a single cell barcode across all UMIs.
	 */
	private BeadSynthesisErrorData buildBeadSynthesisErrorData (final List<UMICollection> umiCollectionList, final StringInterner umiStringCache, final ProgressLogger prog) {
		final String cellBarcode = umiCollectionList.get(0).getCellBarcode();
		BeadSynthesisErrorData bsed = new BeadSynthesisErrorData(cellBarcode);
		for (final UMICollection umis : umiCollectionList) {
			int transcriptCounts = umis.getDigitalExpression(1, 1, false);
			int readCounts = umis.getDigitalExpression(1, 1, true);
			Collection<String> umiCol = umis.getMolecularBarcodes();
			umiCol = getUMIsFromCache(umiCol, umiStringCache);
			bsed.addUMI(umiCol);
			bsed.incrementReads(readCounts);
			bsed.incrementTranscripts(transcriptCounts);
			prog.record(null, 0);
		}
		return bsed;
	}


	/**
	 * Group errors into related neighbor sequences.  Only looks at synthesis errors with 1 base missing.
	 *
	 * @param errorBarcodesWithPositions
	 * @param umiBiasThreshold
	 * @return
	 */
	Map<String, BarcodeNeighborGroup> buildBarcodeNeighborGroups (final Collection<BeadSynthesisErrorData> errorBarcodesWithPositions, final double umiBiasThreshold) {
    	Map<String, BarcodeNeighborGroup> result = new HashMap<>();
    	for (BeadSynthesisErrorData b: errorBarcodesWithPositions) {
    		int polyTErrorPosition = b.getPolyTErrorPosition(umiBiasThreshold);
    		int umiLength = b.getBaseLength();
    		if (polyTErrorPosition==umiLength && b.getErrorType(umiBiasThreshold)==BeadSynthesisErrorType.SYNTH_MISSING_BASE) {
    			String cellBCRoot = padCellBarcode(b.getCellBarcode(), polyTErrorPosition, umiLength);
        		BarcodeNeighborGroup bng = result.get(cellBCRoot);
        		// if the neighbor group is null create and add it.
        		if (bng==null)
    				bng=new BarcodeNeighborGroup(cellBCRoot);
        		bng.addNeighbor(b);
        		result.put(cellBCRoot, bng);
    			}
    	}
    	// List<BarcodeNeighborGroup> neighborGroups= new ArrayList<> (result.values());
    	return result;
    }

	/**
	 * For each problematic cell, replace cell barcodes positions with N.
	 * Take the replaced bases and prepend them to the UMI, and trim the last <X> bases off the end of the UMI.
	 */
	private void cleanBAM (final Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions, final Map<String, String> intendedSequenceMap) {
		log.info("Cleaning BAM");
        final SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(INPUT, true);
		SamHeaderUtil.addPgRecord(headerAndIterator.header, this);

		SAMFileWriter writer= new SAMFileWriterFactory().setCreateIndex(CREATE_INDEX).makeSAMOrBAMWriter(headerAndIterator.header, true, OUTPUT);
		ProgressLogger pl = new ProgressLogger(log);
		for (SAMRecord r: new IterableAdapter<>(headerAndIterator.iterator)) {
			pl.record(r);
			r=padCellBarcodeFix(r, errorBarcodesWithPositions, intendedSequenceMap, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.EXTREME_BASE_RATIO);
			if (r!=null)
				writer.addAlignment(r);
		}
		CloserUtil.close(headerAndIterator.iterator);
		writer.close();
	}


	/**
	 * @return null if the read should not be included in the output BAM.
	 */
	SAMRecord padCellBarcodeFix (final SAMRecord r, final Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions, final Map<String, String> intendedSequenceMap, final String cellBarcodeTag, final String molecularBarcodeTag, final double extremeBaseRatio) {
		String cellBC=r.getStringAttribute(cellBarcodeTag);

		BeadSynthesisErrorData bsed = errorBarcodesWithPositions.get(cellBC);
		if (bsed==null) return (r); // no correction data, no fix.

		// we're only going to fix cells where there's one or more synthesis errors
		BeadSynthesisErrorType bset = bsed.getErrorType(extremeBaseRatio, this.detectPrimerTool, this.EDIT_DISTANCE);
		if (bset==BeadSynthesisErrorType.NO_ERROR) return (r); // no error, return.
		// has an error, not a synthesis error...
		if (bset!=BeadSynthesisErrorType.SYNTH_MISSING_BASE)
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

		// if there's an intended sequence, use that instead of the default padded cell barcode.
		String intendedSeq = intendedSequenceMap.get(cellBC);
		if (intendedSeq!=null)
			cellBCFixed=intendedSeq;

		r.setAttribute(cellBarcodeTag, cellBCFixed);
		r.setAttribute(molecularBarcodeTag, umiFixed);
		return r;
	}


	private void writeSummary(final BeadSynthesisErrorsSummaryMetric summary, final File out) {
		MetricsFile<BeadSynthesisErrorsSummaryMetric, Integer> outFile = new MetricsFile<>();
		outFile.addMetric(summary);
		outFile.addHistogram(summary.getHistogram());
		outFile.write(out);
	}

	private BeadSynthesisErrorsSummaryMetric addDataToSummary (final BeadSynthesisErrorData bsde, final BeadSynthesisErrorsSummaryMetric summary) {
		BeadSynthesisErrorType t = bsde.getErrorType(EXTREME_BASE_RATIO, this.detectPrimerTool, this.EDIT_DISTANCE);
		if (t==null) {
			log.info("STOP");
			t = bsde.getErrorType(EXTREME_BASE_RATIO, this.detectPrimerTool, this.EDIT_DISTANCE);
		}
		summary.NUM_BEADS++;
		switch (t) {
			case SYNTH_MISSING_BASE: summary.SYNTHESIS_MISSING_BASE++; summary.incrementSynthesisMissingBase(bsde.getPolyTErrorPosition(this.EXTREME_BASE_RATIO)); break;
			case PRIMER: summary.PRIMER_MATCH++; break;
			case SINGLE_UMI: summary.SINGLE_UMI_ERROR++; break;
			case FIXED_FIRST_BASE: summary.FIXED_FIRST_BASE++; break;
			case OTHER_ERROR: summary.OTHER_ERROR_COUNT++; break;
			default: summary.NO_ERROR++; break;
		}
		return summary;
	}


	private void writeFile (final PeekableIterator <BeadSynthesisErrorData> iter, final PrintStream out) {
		if (!iter.hasNext()) {
			out.close();
			return;
		}

		BeadSynthesisErrorData first = iter.peek();
		int umiLength = first.getBaseLength();
		writeBadBarcodeStatisticsFileHeader(umiLength, out);
		while (iter.hasNext()) {
			BeadSynthesisErrorData bsed = iter.next();
			writeBadBarcodeStatisticsFileEntry(bsed, out);
		}
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
		header.add(BeadSynthesisErrorType.SYNTH_MISSING_BASE.toString());
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
		line.add(data.getErrorType(EXTREME_BASE_RATIO, this.detectPrimerTool, this.EDIT_DISTANCE).toString());


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


	//TODO: maybe make this part of BaseDistributionMetric
	private String format (final BaseDistributionMetric bdm) {
		List<String> d = new ArrayList<>();

		for (Bases b: Bases.values()) {
			char bb = b.getBase();
			int count = bdm.getCount(bb);
			d.add(Integer.toString(count));
		}
		return StringUtils.join(d, "|");
	}


	/**
	 * Gets a reference to the UMI strings from the cache.  Has the side effect of populating the cache with additional
	 * strings.  This reduces total memory footprint by returning references to repeated strings instead of
	 * holding new objects for the same UMI over and over.
	 * @param umis A list of strings to get references to
	 * @param umiStringCache The cache of strings holding references.
	 */
	private Collection<String> getUMIsFromCache (final Collection<String> umis, final StringInterner umiStringCache) {
		List<String> result = new ArrayList<>(umis.size());
		for (String umi: umis)
			result.add(umiStringCache.intern(umi));
		return (result);
	}

	/**
	 * Take the original cell barcode and UMI, and move bases from the end of the cell barcode to the start of the UMI,
	 * then trim an equal number of bases off the end of the UMI so the length is the same.
	 * Example:
	 * Cell barcode: 		ACGCTCATACAG
	 * UMI: 				TCCTTATT
	 * errorPosition: 		2
	 * New Cell Barcode:	ACGCTCATACNN
	 * New UMI:				AGTCCTTA
	 *
	 * @param cellBarcode The original cell barcode
	 * @param umi The original UMI
	 * @param errorPosition The position in the UMI where the error occurred.
	 */
	public String fixUMI (final String cellBarcode, final String umi, final int errorPosition) {
		// 0 based, from end of cell barcode.
		int badBasesUMI=umi.length()-errorPosition;
		int lastBase = cellBarcode.length();
		int firstBaseToPad = lastBase-badBasesUMI-1;
		String cellBCBases=cellBarcode.substring(firstBaseToPad, cellBarcode.length());

		String umiRemaining=umi.substring(0, errorPosition-1);
		return cellBCBases+umiRemaining;
	}

	/**
	 * Set up the UMI Iterator.
	 * @return
	 */
	public UMIIterator prepareUMIIterator() {
		List<String> barcodes=getCellBarcodes();
		return new UMIIterator(SamFileMergeUtil.mergeInputs(INPUT, false, samReaderFactory),
                this.GENE_EXON_TAG, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.STRAND_TAG, this.READ_MQ,
                false, false, barcodes, true);
	}

	public List<String> getCellBarcodes () {
		if (this.CELL_BC_FILE==null & this.NUM_BARCODES==null) {
			log.info("Selecting cells based on number of transcripts [" + this.MIN_UMIS_PER_CELL +"]");
			return null;
		}

		if (this.CELL_BC_FILE!=null) {
			IOUtil.assertFileIsReadable(this.CELL_BC_FILE);
			List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(this.CELL_BC_FILE);
			log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
			return (cellBarcodes);
		}
		log.info("Gathering barcodes for the top [" + this.NUM_BARCODES +"] cells");
        return new BarcodeListRetrieval().getListCellBarcodesByReadCount(
                SamFileMergeUtil.mergeInputs(INPUT, false, samReaderFactory).iterator,
                this.CELL_BARCODE_TAG, this.READ_MQ, null, this.NUM_BARCODES);
	}


	/**
	 * Picks a number of bases to pad.
	 * If errorPosition =-1, then don't pad any bases.
	 */
	String padCellBarcode (final String cellBarcode, final int errorPosition, final int umiLength) {
		if (errorPosition==-1) return (cellBarcode);

		// 0 based, from end of cell barcode.
		int badBasesUMI=umiLength-errorPosition;
		int lastBase = cellBarcode.length();
		int firstBaseToPad = lastBase-badBasesUMI-1;

		char [] charAr = cellBarcode.toCharArray();
		for (int i=firstBaseToPad; i<lastBase; i++)
			charAr[i]=this.PAD_CHARACTER;
		return new String (charAr);
	}


	 @Override
	    protected String[] customCommandLineValidation() {
	        for (final File input : INPUT)
				IOUtil.assertFileIsReadable(input);
	        IOUtil.assertFileIsWritable(this.OUTPUT_STATS);
	        IOUtil.assertFileIsWritable(this.SUMMARY);
	        if (this.REPORT!=null) IOUtil.assertFileIsWritable(this.REPORT);

	        if (this.CELL_BARCODE_TAG!=null & this.NUM_BARCODES!=null)
	        	throw new TranscriptomeException("Must specify at most one of CELL_BARCODE_TAG or  NUM_BARCODES");

	        if (OUTPUT!=null) IOUtil.assertFileIsWritable(this.OUTPUT);
	        return super.customCommandLineValidation();
	    }


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new DetectBeadSynthesisErrors().instanceMain(args));
	}
}
