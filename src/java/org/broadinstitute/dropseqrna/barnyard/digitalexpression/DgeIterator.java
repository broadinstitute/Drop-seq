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

package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import java.io.BufferedInputStream;
import java.io.File;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeIterator.DgeLine;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import picard.util.TabbedInputParser;

/**
 * Given a DGE file, parse the header and return DGE data a line at a time.
 * Also holds convenience methods to loop up cell barcode information (position, list of cell barcodes)
 * @author nemesh
 *
 */
public class DgeIterator implements CloseableIterator <DgeLine>{

	private final BufferedInputStream inputStream;
	private final DgeHeader dgeHeader;
	private TabbedInputParser parser;
	// map of cell barcode to position.
	private final LinkedHashMap<String, Integer> identifierMap;
	private final String geneColumnLabel;

	public DgeIterator (final File input) {
		this(new BufferedInputStream(IOUtil.openFileForReading(input)), input.getAbsolutePath());
	}

	public DgeIterator (final BufferedInputStream inputStream, final String filename) {
		this.inputStream=inputStream;
		final DgeHeaderCodec headerCodec = new DgeHeaderCodec();
        this.dgeHeader = headerCodec.decode(inputStream, filename);
        // set up a tabbed input parser for the line by line data.
        this.parser = new TabbedInputParser(false, inputStream);

        if (!parser.hasNext())  { // empty file.
            parser.close();
            identifierMap=new LinkedHashMap<>();
            geneColumnLabel="";
            return;
        }

        String [] header = parser.next();
        geneColumnLabel = header[0];

        // provide dependable iteration order for keys.
        LinkedHashMap<String, Integer> tempMap=new LinkedHashMap<>();
        // populate cell barcode map.
        for (int i=1; i<header.length; i++)
        	// position map 0 based.
        	tempMap.put(header[i], i-1);
        this.identifierMap= tempMap;

	}

	/**
	 * Sometimes a DGE can be empty, so there is not even a header line.
	 * @return true if there is nothing but comments in the file.
	 */
	public boolean isEmpty() {
		return this.identifierMap.isEmpty();
	}

	public DgeHeader getDgeHeader () {
		return this.dgeHeader;
	}

	@Override
	public boolean hasNext() {
		return parser.hasNext();
	}

	@Override
	public DgeLine next() {
		if (!hasNext()) return null;
		String [] line = this.parser.next();
		return new DgeLine(this.identifierMap, line);
	}

	public void close() {
		CloserUtil.close(this.parser);
		CloserUtil.close(this.inputStream);
	}

	public void subset (final Set<String> identifiers) {

	}

	/**
	 * Identifiers are ordered in the same way they are originally input.
	 * @return
	 */
	public List<String> getIdentifiers () {
		return new ArrayList<>(this.identifierMap.keySet());
	}
	
	public boolean hasIdentifier (String identifier) {
		return this.identifierMap.containsKey(identifier);
	}

	/**
	 * Get the identifier at the top right hand corner of the matrix that labels the gene column.
	 * This is usually "GENE" for DGE files, but may be different for meta-cells or eQTL data.
	 * @return
	 */
	public String getGeneColumnLabel() {
		return geneColumnLabel;
	}

	public class DgeLine {
		private String gene;
		private double [] expression;
		LinkedHashMap<String, Integer> identifierMap;

		DgeLine (final LinkedHashMap<String, Integer> identifierMap, final String [] line) {
			this.identifierMap=identifierMap;
			expression= new double [line.length-1];
			this.gene=line[0];
			for (int i=1; i<line.length; i++)
				expression[i-1]=Double.parseDouble(line[i]);
		}

		DgeLine (final LinkedHashMap<String, Integer> identifierMap, final String gene, final double [] expression) {
			this.identifierMap=identifierMap;
			this.gene=gene;
			this.expression=expression;
		}

		public String getGene () {
			return this.gene;
		}

		/**
		 * Return the expression of this identifier if the identifier is in the data set, else return null.
		 * @param identifier The donor/sample ID.
		 * @return Return the expression of this identifier if the identifier is in the data set, else return null.
		 */
		public Double getExpression (final String identifier) {
			if (!identifierMap.containsKey(identifier)) return (null);
			Integer pos = identifierMap.get(identifier);			
			return expression[pos];
		}

		public void setExpression (final String identifier, final double value) {
			Integer pos = identifierMap.get(identifier);
			if (pos==null)
				throw new IllegalStateException ("Asked for an identifier ["+identifier+"] that doesn't exist.");
			expression[pos]=value;
		}

		public double [] getExpression () {
			return this.expression;
		}

		public Set<String> getIdentifiers () {
			return identifierMap.keySet();
		}

		LinkedHashMap<String, Integer> getIdentifierMap() {
			return this.identifierMap;
		}

		/**
		 * Is there any identifier with expression greater than 0?
		 * @return return true if at least one identifier is non-zero.
		 */
		public boolean isNonZero () {
			double totalCount = Arrays.stream(this.expression).sum();
			if (BigDecimal.ZERO.equals(new BigDecimal(totalCount))) return false;
			return true;
		}

		public DgeLine subset(final Set<String> identifiers) {
			LinkedHashMap<String, Integer> newMap = new LinkedHashMap<>();

			int index=0;
			double [] exp = new double [identifiers.size()];

			for (String id: this.identifierMap.keySet())
				if (identifiers.contains(id)) {
					newMap.put(id, index);
					exp[index]=this.getExpression(id);
					index++;
				}

			DgeLine l = new DgeLine(newMap, this.getGene(), exp);
			return l;
		}

	}




}
