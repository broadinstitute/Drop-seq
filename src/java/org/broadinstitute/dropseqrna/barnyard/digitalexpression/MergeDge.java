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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.AbstractTripletDgeWriterClp;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeader;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderMerger;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools.DGEMatrix;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Merges multiple Digital Gene Expression files into a single file",
        oneLineSummary = "Merge Digital Gene Expression files",
        programGroup = DropSeq.class
)
public class MergeDge
		extends AbstractTripletDgeWriterClp {

	private final Log log = Log.getInstance(MergeDge.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The DGE file to merge.  Can be used an arbirary number of times.")
	public List<File> INPUT;

	@Argument(doc = "The Prefix to add to every cell barcode so that cell barcodes seen "
			+ "in multiple DGEs can be differentiated.  A typical setting would be 'EXP1_', yielding barcodes like EXP1_ACTGACCGTTTG.  Must be"
			+ "invoked as many times as INPUT." )
	public List<String> PREFIX;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output DGE file")
	public File OUTPUT;

	@Argument(doc="The expression value to set for missing values.  When two DGEs are merged, if one contains a gene annotation and the other does not, then all cells"
			+ "missing that gene have this expression value set.")
	public Float MISSING_VALUE=0f;

	@Argument(doc="Should the output DGE be formatted as integers?  Set this to true when merging normal DGE files.  Set to false "
			+ "if you're merging DGE files that you have performed some other transformations that change the expression results from "
			+ "quantized integers to a continuous variable")
	public Boolean INTEGER_FORMAT=true;

	@Argument(shortName = "FORMAT", doc="Write in standard DGE tabular format, or Drop-seq Matrix Market format")
	public DGEMatrix.FileFormat OUTPUT_FORMAT = DGEMatrix.FileFormat.DENSE;

	@Argument(shortName = "H", doc="If true, write a header in the DGE file.  If true, OUTPUT_FORMAT must == DENSE")
	public boolean OUTPUT_HEADER=true;

	@Argument(doc="How strict to be when merging DGE headers")
	public DgeHeaderMerger.Stringency HEADER_STRINGENCY = DgeHeaderMerger.Stringency.STRICT;

	@Override
	protected String[] customCommandLineValidation() {
		final String[] superErrors = super.customCommandLineValidation();
		if (OUTPUT_FORMAT == DGEMatrix.FileFormat.DENSE || !OUTPUT_HEADER) {
			return superErrors;
		} else {
			final ArrayList<String> list = new ArrayList<>(1);
			if (superErrors != null) {
				for (final String msg: superErrors) {
					list.add(msg);
				}
			}
			list.add("OUTPUT_HEADER==true is not supported if OUTPUT_FORMAT!=DENSE");
			return list.toArray(new String[list.size()]);
		}
	}

	@Override
	protected int doWork() {
        if (PREFIX == null)
            PREFIX = new ArrayList<String>();

        INPUT = FileListParsingUtils.expandFileList(INPUT);

		if (!PREFIX.isEmpty() && INPUT.size()!=PREFIX.size()) {
			log.error("Must have the same number of INPUT and PREFIX parameters, or else no PREFIXes can be specified");
			return 1;
		}
		for (final File f: INPUT)
			IOUtil.assertFileIsReadable(f);
		IOUtil.assertFileIsWritable(OUTPUT);

		if (INPUT.size()==0) return 0;

		mergeSparseMatrix();

		return 0;
	}

	private void mergeSparseMatrix() {
		DGEMatrix dge = readMatrix(0);
		for (int i=1; i<INPUT.size(); i++) {
            DGEMatrix dgeOther = readMatrix(i);
			log.info("Merging " + INPUT.get(i).getName());
			dge=dge.merge(dgeOther);
		}

        final DgeHeader mergedHeader;
        if (OUTPUT_HEADER) {
            mergedHeader = DgeHeaderMerger.mergeDgeHeaders(INPUT, PREFIX, HEADER_STRINGENCY);
			mergedHeader.addCommand(getCommandLine());
        } else {
            mergedHeader = null;
        }
		log.info("Writing output");
		dge.writeFile(this.OUTPUT, INTEGER_FORMAT, OUTPUT_FORMAT, mergedHeader);

		maybeWriteNamesFiles(dge);
	}

	private DGEMatrix readMatrix(final int which) {
        final String cellBarcodePrefix;
        if (PREFIX.isEmpty())
			cellBarcodePrefix = "";
		else
			cellBarcodePrefix = PREFIX.get(which);
        return DGEMatrix.parseFile(INPUT.get(which), cellBarcodePrefix);
    }

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new MergeDge().instanceMain(args));
	}

	protected void maybeWriteNamesFiles(DGEMatrix dge) {
		maybeWriteNamesFiles(dge.getGenes(), dge.getCellBarcodes());
	}
}
