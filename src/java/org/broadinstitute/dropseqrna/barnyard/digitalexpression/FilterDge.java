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

import java.io.BufferedInputStream;
import java.io.File;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.metrics.MetricsFile;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.DigitalExpression;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.RetainRemoveList;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintWriter;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedInputParser;

@CommandLineProgramProperties(
        summary = "Filter a DGE to remove genes or cells.",
        oneLineSummary = "Filter a DGE to remove genes or cells.",
        programGroup = DropSeq.class
)
public class FilterDge extends CommandLineProgram{

	private final Log log = Log.getInstance(FilterDge.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The DGE file to process.")
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output DGE file")
	public File OUTPUT;

	@Argument(doc="A file with a single column containing the names of genes to remove.", optional=true)
	public File GENES_REMOVE;

	@Argument(doc="A file with a single column containing the names of cell barcodes to remove.", optional=true)
	public File CELLS_REMOVE;

	@Argument(doc="A file with a single column containing the names of genes to retain. Elements in this list trump GENES_REMOVE", optional=true)
	public File GENES_RETAIN;

	@Argument(doc="A file with a single column containing the names of cell barcodes to retain.  Elements in this list trump CELLS_REMOVE.", optional=true)
	public File CELLS_RETAIN;

	@Argument(shortName = "H", doc="If true, write a header in the DGE file")
	public boolean OUTPUT_HEADER=true;

	@Argument(doc="If OUTPUT_SUMMARY is specified, read this file to find values that can't be computed from DGE. " +
			"If not specified, but OUTPUT_SUMMARY is specified, then NUM_GENIC_READS will be 0.", optional = true)
	public File INPUT_SUMMARY;

	@Argument(doc="If specified, write a digital expression summary file corresponding to the DGE.", optional = true)
	public File OUTPUT_SUMMARY;

	private final String GENE_HEADER="GENE";
	private final int PROGRESS_INTERVAL=1000;

	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(this.INPUT);
		IOUtil.assertFileIsWritable(this.OUTPUT);

		if (OUTPUT_SUMMARY != null) {
			IOUtil.assertFileIsWritable(this.OUTPUT_SUMMARY);
			if (INPUT_SUMMARY != null) {
				IOUtil.assertFileIsReadable(INPUT_SUMMARY);
			}
		}

		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));

		final BufferedInputStream in = new BufferedInputStream(IOUtil.openFileForReading(INPUT));

        final DgeHeaderCodec headerCodec = new DgeHeaderCodec();
        final DgeHeader dgeHeader = headerCodec.decode(in, INPUT.getAbsolutePath());
		dgeHeader.addCommand(getCommandLine());

		TabbedInputParser parser = new TabbedInputParser(false, in);
		if (!parser.hasNext()) {
			log.error("No lines in input file" + this.INPUT);
			parser.close();
			return 1;
		}
		// get a map of the header names to their positions.
		String [] header = parser.next();
		// header mapped to positions, before any cells are filtered out.
		Map<String, Integer> headerPositionMap = getSamplePositionMap(header);

		// what barcodes are in the output?  The first entry in the array is the gene name, so remove that.
		String [] cells = Arrays.copyOfRange(header, 1, header.length);
		// this is the list of cells to iterate over. calculate this just once.
		List<String> cellsRetain = getElementsToRetain(Arrays.asList(cells), this.CELLS_REMOVE, this.CELLS_RETAIN);
		log.info("Retaining [", cellsRetain.size(), "] cells");

        if (OUTPUT_HEADER) {
            final PrintWriter writer = new ErrorCheckingPrintWriter(out);
            headerCodec.encode(writer, dgeHeader);
            writer.flush();
        }

        final Map<String, DigitalExpression.DESummary> summaries = maybePopulateSummaries(cellsRetain);

		// write the header of the output file.
		writeLine(cellsRetain, out, true);

		// parse the bulk of the file, write out lines
		Set<String> genesRetain = getRetainRemoveSet(this.GENES_RETAIN);
		Set<String> genesRemove = getRetainRemoveSet(this.GENES_REMOVE);
		int counter=0;
		while (parser.hasNext()) {
			String [] line = parser.next();
			List<String> lineParsed = buildBodyLine(line, genesRetain, genesRemove, cellsRetain, headerPositionMap, summaries);
			if (lineParsed!=null)
				writeLine(lineParsed, out, false);
			if (counter%PROGRESS_INTERVAL==0) log.info("Processed [" + counter +"] lines");
			counter++;
		}
		// cleanup.
		parser.close();
		out.close();

		if (this.OUTPUT_SUMMARY!=null) {
			MetricsFile<DigitalExpression.DESummary, Integer> summaryMetricsFile = getMetricsFile();
			DigitalExpression.writeSummary(summaries.values(), summaryMetricsFile, this.OUTPUT_SUMMARY);
		}


		return 0;
	}

	/**
	 * Builds one line of output - a single gene with the subset of cells that are retained.
	 * Can return null if the gene should not be retained.
	 */
	private List<String> buildBodyLine (final String [] line, final Set<String> genesRetain,
										final Set<String> genesRemove, final List<String> cellsRetain,
										final Map<String, Integer> headerPositionMap,
										final Map<String, DigitalExpression.DESummary> summaries) {

		String gene = line[0];
		boolean retain = new RetainRemoveList<String>().retainElement(gene, genesRetain, genesRemove);
		if (!retain) return null;
		List<String> result = new ArrayList<>();
		List<Double> values = new ArrayList<>();
		// only write the gene 
		if (cellsRetain.size()>0)
			result.add(gene);
		for (String s: cellsRetain) {
			// get the position in the unfiltered line
			int pos = headerPositionMap.get(s);
			// result element
			String e = line[pos];
			double d = Double.parseDouble(line[pos]);
			if (summaries != null && d > 0) {
				DigitalExpression.DESummary summary = summaries.get(s);
				++summary.NUM_GENES;
				summary.NUM_TRANSCRIPTS += (int)d;
			}
			values.add(d);
			result.add(e);
		}
		// check if any values > 0.  If so, return the entire gene.  If not, return null.
		for (Double d: values)
			if (d>0) return result;

		return null;
	}

	private Map<String, Integer> getSamplePositionMap (final String [] header) {
		Map<String, Integer> result = new HashMap<>();
		for (int i=0; i<header.length; i++)
			result.put(header[i], i);
		return (result);
	}

	/**
	 * Selects the set of elements(cells) to retain.
	 * Maintains the order of the DGE elements(cells)
	 * @param elements The original list of cells
	 * @param remove A file containing a list of cells to remove, can be null.
	 * @param retain A file containing a list of cells to retain, can be null.
	 */
	public List<String> getElementsToRetain (final List<String> elements, final File remove, final File retain) {

		// get final set to remove
		Set <String> toRemove = getRetainRemoveSet(remove);
		Set <String> toRetain = getRetainRemoveSet(retain);

		// in the case where the input retain file is not null, but has 0 entries, there is nothing to retain. 
		if (retain!=null & toRetain.size()==0) 
			return Collections.emptyList();		

		RetainRemoveList<String> rrl = new RetainRemoveList<>();
		return (rrl.getElementsToRetain(elements, toRemove, toRetain));
	}


	private Set<String> getRetainRemoveSet (final File f) {
		if (f==null)
			return (new HashSet<>(0));
		Set<String >result = new HashSet<> (ParseBarcodeFile.readCellBarcodeFile(f));
		return (result);
	}

	private void writeLine (final List<String> line, final PrintStream out, final boolean addGeneHeader) {
		List<String> lineOut = new ArrayList<>(line);

		// short-circuit: only write the line if there is at least one donor. 		
		if (line.size()==0) 
			return;
		
		if (addGeneHeader)
			lineOut.add(0, this.GENE_HEADER);
		String b = StringUtils.join(lineOut, "\t");
		out.println(b);
	}

	/**
	 * Optionally initialize map of cell_barcode => DESummary
	 * @param cellsRetain the cells that will be output
	 * @return null if OUTPUT_SUMMARY is null, else a map of DESummaries, maybe with NUM_GENIC_READS initialized
	 */
	private Map<String, DigitalExpression.DESummary> maybePopulateSummaries(final List<String> cellsRetain) {
		if (OUTPUT_SUMMARY == null) {
			return null;
		} else {
			final Map<String, DigitalExpression.DESummary> summaries = new HashMap<>();
			for (final String cell : cellsRetain) {
				summaries.put(cell, new DigitalExpression.DESummary(cell));
			}
			if (INPUT_SUMMARY != null) {
				final List<DigitalExpression.DESummary> oldSummaries = MetricsFile.readBeans(INPUT_SUMMARY);
				for (final DigitalExpression.DESummary oldSummary : oldSummaries) {
					final DigitalExpression.DESummary summary = summaries.get(oldSummary.CELL_BARCODE);
					if (summary != null) {
						summary.NUM_GENIC_READS = oldSummary.NUM_GENIC_READS;
					}
				}
			}
			return summaries;
		}
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new FilterDge().instanceMain(args));
	}
}
