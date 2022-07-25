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
package org.broadinstitute.dropseqrna.barnyard;

import java.io.File;
import java.util.List;

import org.broadinstitute.barclay.argparser.Argument;

import picard.cmdline.StandardOptionDefinitions;

public abstract class DGECommandLineBase extends GeneFunctionCommandLineBase {

	@Argument(doc="The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
	public String CELL_BARCODE_TAG="XC";
	
	@Argument(doc="The molecular barcode tag.")
	public String MOLECULAR_BARCODE_TAG="XM";
	
	@Argument(doc="The edit distance that molecular barcodes should be combined at within a gene.")
	public Integer EDIT_DISTANCE=1;
	
	@Argument(doc="The map quality of the read to be included.")
	public Integer READ_MQ=10;
	
	@Argument(doc="The minimum number of reads a molecular barcode should have to be considered.  This is done AFTER edit distance collapse of barcodes.")
	public Integer MIN_BC_READ_THRESHOLD=0;
	
	@Argument(doc="Gather up all cell barcodes that have more than some number of reads.", optional=true)
	public Integer MIN_NUM_READS_PER_CELL=null;
	
	@Argument(doc="The minumum number of genes for a cell barcode to be reported.", optional=true)
	public Integer MIN_NUM_GENES_PER_CELL=null;
	
	@Argument(doc="The minumum number of transcripts for a cell barcode to be reported.", optional=true)
	public Integer MIN_NUM_TRANSCRIPTS_PER_CELL=null;
	
	@Argument(doc="Number of cells that you think are in the library.  The largest <X> barcodes are used.", optional=true)
	public Integer NUM_CORE_BARCODES=null;
	
	@Argument(doc="Override CELL_BARCODE and MIN_NUM_READS_PER_CELL, and process reads that have the cell barcodes in this file instead.  The file has 1 column with no header.", optional=true)
	public File CELL_BC_FILE=null;
	
	@Argument (doc="Drop UMIs within a cell/gene pair that occur less than the average number of reads*<FILTER_FREQ> for all UMIs in the cell/gene pair.  " +
			"For example, if you had on average 1000 reads per UMI and a UMI with 1-10 reads, those small UMIs would be removed when this frequency was set to 0.01." +
			"This is off by default.  A sensible value might be 0.01.")
	public double RARE_UMI_FILTER_THRESHOLD=0;


}
