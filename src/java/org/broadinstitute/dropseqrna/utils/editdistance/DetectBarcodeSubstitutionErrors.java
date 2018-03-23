/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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


package org.broadinstitute.dropseqrna.utils.editdistance;

import java.io.File;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.beadsynthesis.BeadSynthesisErrorData;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.StringUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(summary = "Collapses umambiguously related small barcodes into larger neighbors.  Unambiguously related barcodes are situations where a smaller barcode"
		+ "only has 1 neighbor within the edit distance threshold, so the barcode can not be collapsed to the wrong neighbor.  These sorts of errors can be due to problems with barcode synthesis."
		+ "Ambiguous barcodes are situations where a smaller neighbor has multiple larger neighbors.  These barcodes can be optionally filtered.)",
oneLineSummary = "Collaps umambiguously related small barcodes into larger neighbors.)",
programGroup = DropSeq.class)

public class DetectBarcodeSubstitutionErrors extends CommandLineProgram{

	private final Log log = Log.getInstance(DetectBarcodeSubstitutionErrors.class);
	private ProgressLogger pl = new ProgressLogger(this.log);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input DropSeq BAM file to analyze", minElements = 1)
	public List<File> INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BAM file with cell barcodes collapsed.")
	public File OUTPUT;

	@Argument(doc="Output report detailing which barcodes were merged, and what the position of the substitution and intended/changed bases were for each pair of barcordes merged.", optional=true)
	public File OUTPUT_REPORT;

	@Argument(doc="The output barcode tag for the newly collapsed barcodes.  Defaults to the CELL_BARCODE_TAG if not set.", optional=true)
	public String OUT_CELL_BARCODE_TAG;

	@Argument(doc="Remove smaller barcodes that map at the edit distance to multiple larger barcodes.")
	public Boolean FILTER_AMBIGUOUS=true;

	@Argument(doc="The cell barcode tag.")
	public String CELL_BARCODE_TAG="XC";

	@Argument(doc="The molecular barcode tag.")
	public String MOLECULAR_BARCODE_TAG="XM";

	@Argument(doc="The Gene/Exon tag")
	public String GENE_EXON_TAG="GE";

	@Argument(doc="The strand of the gene(s) the read overlaps.  When there are multiple genes, they will be comma-separated.")
	public String STRAND_TAG="GS";

	@Argument (doc="The minimum number of UMIs required to consider a cell barcode for collapse.  Setting this number higher speeds up cleanup.  Very small barcodes will not contribute many UMIs, so are not a useful return on investment.  Suggested values range from 20 to 200.")
	public Integer MIN_UMIS_PER_CELL=20;

	@Argument (doc="The amount of bias (all UMIs for a cell have the same base) at which a cell barcode is considered biased?", optional=true)
	public Double UMI_BIAS_THRESHOLD=0.8;

	@Argument (doc="Which base to scan for UMI bias.  This is typically the last base of the UMI.  If set to null, program will use the last base of the UMI.", optional=true)
	public Integer UMI_BIAS_BASE=null;

	@Argument(doc="The edit distance to collapse barcodes")
	public Integer EDIT_DISTANCE=1;

	@Argument(doc="Read quality filter.  Filters all reads lower than this mapping quality.  Defaults to 10.  Set to 0 to not filter reads by map quality.")
	public Integer READ_MQ=10;

	@Argument(doc="Number of threads to use.  Defaults to 1.")
	public int NUM_THREADS=1;

	@Override
	protected int doWork() {
		for (final File input : INPUT)
			IOUtil.assertFileIsReadable(input);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (this.OUTPUT_REPORT!=null) IOUtil.assertFileIsWritable(this.OUTPUT_REPORT);

        // set the cell barcode tag for the output if not set.
        if (this.OUT_CELL_BARCODE_TAG==null) this.OUT_CELL_BARCODE_TAG=this.CELL_BARCODE_TAG;

        // build up the UMI per cell data set.
        Collection <String> cellBarcodes = null;
        UMIIterator umiIterator = new UMIIterator(SamFileMergeUtil.mergeInputs(this.INPUT, false),
                this.GENE_EXON_TAG, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.STRAND_TAG,
                this.READ_MQ, false, true, cellBarcodes, true);

        // get list of barcodes that have enough UMIs, and are not polyT biased.
        ObjectCounter<String> umiCounts=getUMIsPerCell(umiIterator, this.MIN_UMIS_PER_CELL, this.UMI_BIAS_BASE, this.UMI_BIAS_THRESHOLD);

        // how do they collapse bottom up?
        MapBarcodesByEditDistance med = new MapBarcodesByEditDistance(true, this.NUM_THREADS, 10000);
        log.info("Starting Barcode Collapse of [" + umiCounts.getSize()+ "] barcodes");
        BottomUpCollapseResult result= med.bottomUpCollapse(umiCounts, this.EDIT_DISTANCE);
        log.info("Barcode Collapse Complete - ["+ result.getUnambiguousSmallBarcodes().size() + "] barcodes collapsed");
        umiIterator.close();

        // write report on which substitutions were found
        if (this.OUTPUT_REPORT!=null) writeReport(result, umiCounts);

        // perform repair/filtering.
        repairBAM(result);
        log.info("Finished");
		return 0;
	}

	private void writeReport (final BottomUpCollapseResult result, final ObjectCounter<String> umiCounts) {
		PrintStream outReport = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.OUTPUT_REPORT));
		/// write header
		String [] header= {"intended_barcode", "neighbor_barcode", "intended_size", "neighbor_size", "position", "intended_base", "neighbor_base"};
		outReport.println(StringUtil.join("\t", header));


		Iterator<String> smalls = result.getUnambiguousSmallBarcodes().iterator();
		while (smalls.hasNext()) {
			String small=smalls.next();
			String large = result.getLargerRelatedBarcode(small);
			int [] pos = HammingDistance.getHammingDistanceChangePositions(small, large);
			if (pos.length!=1)
				new IllegalArgumentException("Strings don't have edit distance 1!");
			String intendedBase = large.substring(pos[0], pos[0]+1);
			String neighborBase = small.substring(pos[0], pos[0]+1);
			// output position is one based!
			String [] body = {large, small, Integer.toString(umiCounts.getCountForKey(large)), Integer.toString(umiCounts.getCountForKey(small)),
					Integer.toString(pos[0]+1), intendedBase, neighborBase};

			outReport.println(StringUtil.join("\t", body));
		}


		CloserUtil.close(outReport);
	}


	private void repairBAM (final BottomUpCollapseResult result) {

		final SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(this.INPUT, false);
		SamHeaderUtil.addPgRecord(headerAndIterator.header, this);
		headerAndIterator.header.addComment("Bottom-up edit distance collapse tag " + this.CELL_BARCODE_TAG +" with edit distance " + this.EDIT_DISTANCE+ " filtering ambiguous neighbors=" + this.FILTER_AMBIGUOUS);
		SAMFileWriter writer= new SAMFileWriterFactory().setCreateIndex(CREATE_INDEX).makeSAMOrBAMWriter(headerAndIterator.header, true, OUTPUT);

		ProgressLogger pl = new ProgressLogger(log);
		log.info("Repairing BAM");
		for (SAMRecord r: new IterableAdapter<>(headerAndIterator.iterator)) {
			pl.record(r);
			r=repairBarcode(r, result);
			if (r!=null)
				writer.addAlignment(r);
		}
		CloserUtil.close(headerAndIterator.iterator);
		writer.close();
		log.info("Repair Complete");
	}


	/**
	 * Repairs the cell barcode, which can mean one of the following:
	 * 1) If the cell barcode is ambiguous and you're filtering ambiguous, return null
	 * 2) If the cell barcode is the smaller neighbor of the pair, set the tag to be the larger neighbor
	 * 2a) If using a different output cell barcode, tag each read with the new tag and either the unchanged cell barcode or the new neighbor.
	 * @param r The input read
	 * @param result The cell Barcode collapse result so we can replace smaller barcodes with larger ones to merge them.
	 * @return The modified read, or null if the read should be filtered from the data.
	 */
	private SAMRecord repairBarcode (final SAMRecord r, final BottomUpCollapseResult result) {
		String cellBarcode=r.getStringAttribute(this.CELL_BARCODE_TAG);
		// check filter first as it's faster.
		boolean ambiguous=result.isAmbiguousBarcode(cellBarcode);
		if (this.FILTER_AMBIGUOUS && ambiguous) return null;
		String newCellBarcode=result.getLargerRelatedBarcode(cellBarcode);
		// if not null, we have something to repair.
		if (newCellBarcode!=null)
			r.setAttribute(this.OUT_CELL_BARCODE_TAG, newCellBarcode);
		else // add the cell barocode tag with the old value.  If you're writing out to a new tag this is important.
			r.setAttribute(this.OUT_CELL_BARCODE_TAG, cellBarcode);
		return r;
	}


	/**
	 * Get the count of UMIs per cell.  Filter out barcodes that have barcode synthesis errors (more polyT than allowed at a specific position, usually the last base of the UMI.)
	 * @param iter The UMIterator
	 * @param minUMIsPerCell Cell must have at least this many transcripts to be added
	 * @param polyTPosition What position should be checked for a barcode synthesis error.  If set to null, uses the last base of the UMI.
	 * @param polyTThreshold How much bias excludes a cell [0-1].
	 * @return
	 */
	public ObjectCounter<String> getUMIsPerCell (final UMIIterator iter, final int minUMIsPerCell, Integer polyTPosition, final double polyTThreshold) {
		if (polyTThreshold > 1 | polyTThreshold<0)
			throw new IllegalArgumentException("PolyT Threshold must be between 0 and 1.");

		GroupingIterator<UMICollection> groupingIterator = new GroupingIterator<>(iter,
                new Comparator<UMICollection>() {
                    @Override
                    public int compare(final UMICollection o1, final UMICollection o2) {
                        return o1.getCellBarcode().compareTo(o2.getCellBarcode());
                    }
                });

		ObjectCounter<String> result = new ObjectCounter<>();
		int counter=0;
		int polyTBiasedBarcodes=0;
		log.info("Gathering UMI counts per cell and filtering out UMI biased barcodes as appropriate");
        for (final List<UMICollection> umiCollectionList : groupingIterator) {
            final String cellBarcode = umiCollectionList.get(0).getCellBarcode();
            BeadSynthesisErrorData bsed = new BeadSynthesisErrorData(cellBarcode);
            for (final UMICollection umis : umiCollectionList) {
            	int transcriptCounts= umis.getDigitalExpression(1, 1, false);
            	int readCounts = umis.getDigitalExpression(1, 1, true);
                Collection<String> umiCol = umis.getMolecularBarcodes();
            	if (polyTPosition==null && umiCol.size()>0) {
            		polyTPosition=umiCol.iterator().next().length();
            		log.info("Auto-discovered last base of UMI ["+polyTPosition +"]");
            	}

                bsed.addUMI(umiCol);
                bsed.incrementReads(readCounts);
                bsed.incrementTranscripts(transcriptCounts);
                counter++;
                if (counter%1000000==0) log.info("Processed [" + counter + "] Cell/Gene UMIs.");
            }
            // check if the barcode is polyT biased at the last base.
            boolean polyTBiased = bsed.isPolyTBiasedAtPosition(polyTPosition, polyTThreshold);
            if (polyTBiased) polyTBiasedBarcodes++;
            if (bsed.getNumTranscripts() >= minUMIsPerCell && !polyTBiased)
				result.incrementByCount(bsed.getCellBarcode(), bsed.getNumTranscripts());
        }
        log.info("Finished gathering a list of cell barcodes to collapse");
        log.info("Discovered ["+polyTBiasedBarcodes+"] barcodes that were excluded due to incomplete synthesis and had T bias at last base of UMI [" + polyTPosition+"]");
		return result;
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new DetectBarcodeSubstitutionErrors().instanceMain(args));
	}

}
