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
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.beadsynthesis.BeadSynthesisErrorData;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(summary = "Collapses umambiguously related small barcodes into larger neighbors.  Unambiguously related barcodes are situations where a smaller barcode"
		+ "only has 1 neighbor within the edit distance threshold, so the barcode can not be collapsed to the wrong neighbor.  These sorts of errors can be due to problems with barcode synthesis."
		+ "Ambiguous barcodes are situations where a smaller neighbor has multiple larger neighbors.  These barcodes can be optionally filtered.)",
oneLineSummary = "Collaps umambiguously related small barcodes into larger neighbors.)",
programGroup = DropSeq.class)

public class BottonUpCollapse extends CommandLineProgram{

	private final Log log = Log.getInstance(CollapseBarcodesInPlace.class);
	private ProgressLogger pl = new ProgressLogger(this.log);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input DropSeq BAM file to analyze", minElements = 1)
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BAM file with cell barcodes collapsed.")
	public File OUTPUT;

	@Argument(doc="The output barcode tag for the newly collapsed barcodes.  Defaults to the CELL_BARCODE_TAG if not set.", optional=true)
	public String OUT_CELL_BARCODE_TAG;

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
		// TODO Auto-generated method stub
		IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        Collection <String> cellBarcodes = null;

        UMIIterator umiIterator = new UMIIterator(SamFileMergeUtil.mergeInputs(Collections.singletonList(this.INPUT), false),
                this.GENE_EXON_TAG, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.STRAND_TAG,
                this.READ_MQ, false, true, cellBarcodes);


		MapBarcodesByEditDistance med = new MapBarcodesByEditDistance(true, this.NUM_THREADS, 10000);


		return 0;
	}


	/**
	 * Get the count of UMIs per cell.  Filter out barcodes that have barcode synthesis errors (more polyT than allowed at a specific position, usually the last base of the UMI.)
	 * @param iter The UMIterator
	 * @param minUMIsPerCell Cell must have at least this many transcripts to be added
	 * @param polyTPosition What position should be checked for a barcode synthesis error
	 * @param polyTThreshold How much bias excludes a cell [0-1].
	 * @return
	 */
	public ObjectCounter<String> getUMIsPerCell (final UMIIterator iter, final int minUMIsPerCell, final int polyTPosition, final double polyTThreshold) {
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

        for (final List<UMICollection> umiCollectionList : groupingIterator) {
            final String cellBarcode = umiCollectionList.get(0).getCellBarcode();
            BeadSynthesisErrorData bsed = new BeadSynthesisErrorData(cellBarcode);
            for (final UMICollection umis : umiCollectionList) {
            	int transcriptCounts= umis.getDigitalExpression(1, 1, false);
            	int readCounts = umis.getDigitalExpression(1, 1, true);
                Collection<String> umiCol = umis.getMolecularBarcodes();
                bsed.addUMI(umiCol);
                bsed.incrementReads(readCounts);
                bsed.incrementTranscripts(transcriptCounts);
                counter++;
                if (counter%1000000==0) log.info("Processed [" + counter + "] Cell/Gene UMIs.");
            }
            // check if the barcode is polyT biased at the last base.
            boolean polyTBiased = bsed.isPolyTBiasedAtPosition(polyTPosition, polyTThreshold);
            if (bsed.getUMICount() >= minUMIsPerCell & !polyTBiased)
				result.incrementByCount(bsed.getCellBarcode(), bsed.getNumTranscripts());
        }
		return result;
	}

}
