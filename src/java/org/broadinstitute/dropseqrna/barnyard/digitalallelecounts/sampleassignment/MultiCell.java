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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment;

import htsjdk.samtools.util.IOUtil;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import picard.util.BasicInputParser;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class MultiCell {

	private final List<String> cellBarcodes;
	private final String mergedCellName;

	public MultiCell(final List<String> cellBarcodes) {
		this.cellBarcodes=cellBarcodes;
		mergedCellName=StringUtils.join(cellBarcodes, ":");
	}

	public List<String> getCellBarcodes () {
		return this.cellBarcodes;
	}

	public String getMultiCellName() {
		return mergedCellName;
	}

	public static List<MultiCell> parseFile(final File f) {
		List<MultiCell> result = new ArrayList<>();
		BasicInputParser parser = new BasicInputParser(false, 1, f);
		while (parser.hasNext()) {
			String[] line = parser.next();
			String [] cellIDs=line[0].split(":");
			MultiCell d = new MultiCell(Arrays.asList(cellIDs));
			result.add(d);
		}
		parser.close();
		return (result);
	}

	/**
	 * Get a map from each individual cell name to the multi-cell name.
	 * @return A map with the key of a single cell, and a value a list of the multi cell name(s).
	 * Multiple keys will have the same value.
	 */
	public static Map<String, List<String>> getMap (final Collection<MultiCell> multiCells) {
		Map<String, List<String>> result = new HashMap<>();

		for (MultiCell m: multiCells)
			for (String singleCell: m.cellBarcodes) {
				List<String> l= result.get(singleCell);
				if (l==null)
					l= new ArrayList<>();
				l.add(m.getMultiCellName());
				result.put(singleCell, l);
			}
		return result;
	}

	public static List<MultiCell> generateMultiCells (final Collection<String> cellBarcodes, final int numMultiCells, final int multiplicity) {
		if ((numMultiCells*multiplicity)>cellBarcodes.size())
			throw new IllegalArgumentException("Requested ["+ numMultiCells+" multi-cells with multiplicity ["+ multiplicity +"] "
					+ "from ["+ cellBarcodes.size()+"] You need more cell barcodes or fewer multi-cells / lower multiplicity");


		List<MultiCell> result = new ArrayList<>(numMultiCells);

		// make a new list holding the cell barcodes and shuffle it.
		List <String> cb = new ArrayList<>(cellBarcodes);
		Collections.shuffle(cb);
		Iterator<String> iter = cb.iterator();

		for (int i=0; i<numMultiCells; i++) {
			// generate single multiCell
			List<String> currentCB = new ArrayList<>(multiplicity);
			for (int j=0; j<multiplicity; j++)
				currentCB.add(iter.next());
			MultiCell mc = new MultiCell(currentCB);
			result.add(mc);
		}

		return (result);

	}

	public static void writeToFile (final List<MultiCell> multiCells, final File outFile) {
		PrintStream writer=  new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outFile));
		for (MultiCell mc: multiCells)
			writer.println(mc.getMultiCellName());
		writer.close();
	}

	@Override
	public String toString() {
		return "MultiCell [" + this.mergedCellName +"] individual cell list " + this.cellBarcodes.toString();
	}

}
