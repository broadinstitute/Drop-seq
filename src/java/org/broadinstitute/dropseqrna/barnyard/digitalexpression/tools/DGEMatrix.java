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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeader;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderCodec;
import org.broadinstitute.dropseqrna.matrixmarket.MatrixMarketConstants;
import org.broadinstitute.dropseqrna.matrixmarket.MatrixMarketReader;
import org.broadinstitute.dropseqrna.matrixmarket.MatrixMarketWriter;
import org.la4j.Matrix;
import org.la4j.Vectors;
import org.la4j.iterator.MatrixIterator;
import org.la4j.matrix.SparseMatrix;
import org.la4j.matrix.sparse.CRSMatrix;
import org.la4j.vector.functor.VectorAccumulator;
import picard.util.BasicInputParser;
import picard.util.TabbedInputParser;

import java.io.*;
import java.util.*;

/**
 * This class is designed to hold a matrix of gene expression, addressable by gene name and cell barcode
 * This class supports reading and writing in the standard format for DGE data, merging of DGEMatrixes in place,
 * removing cells/genes and associated expression data, and applying arbitrary in-place transforms of the expression data matrix
 * using the MatrixTransformI interface.
 * @author nemesh
 *
 */
public class DGEMatrix {
	public enum FileFormat {
		DENSE, MM_SPARSE, MM_SPARSE_10X
	}

	private static final Log log = Log.getInstance(DGEMatrix.class);

    // holds the row names of the matrix.
	private final Map<String, Integer> geneMap;

	// holds a map of the cell barcode to the position in the expression array
	// this is consistent across all cells.
	private final Map<String, Integer> cellBarcodeMap;
	// Computed as necessary, and invalidate when appropriate
	private List<String> cellBarcodeCache;

	// the expression data in a matrix
	private Matrix expressionMatrix;
	private final double SPARSE_VALUE=0;

	public DGEMatrix(final List<String> cellBarcodes, final List<String> geneNames, final double [] [] expressionMatrix) {
		this(cellBarcodes, geneNames, from2DCtorHelper(cellBarcodes, geneNames, expressionMatrix));
	}

	private static Matrix from2DCtorHelper(final List<String> cellBarcodes,
										   final List<String> geneNames,
										   final double [] [] expressionMatrix) {
		if (expressionMatrix.length!=geneNames.size()) {
			throw new IllegalArgumentException("Rows of expression matrix [" + expressionMatrix.length +
					"] not equal to number of genes [" + geneNames.size() + "]");
		}
		if (expressionMatrix[0].length!=cellBarcodes.size()) {
			throw new IllegalArgumentException("Columns of expression matrix [" + expressionMatrix[0].length +
					"] not equal to number of cells [" + cellBarcodes.size() + "]");
		}
		return CRSMatrix.from2DArray(expressionMatrix);
	}

	/**
	 * Sets up the matrix, but with all expression set to 0.
	 * @param geneNames A list of gene names that are the rows of the matrix
	 * @param cellBarcodes A list of cell barcodes that are columns of the matrix.
	 */
	public DGEMatrix(final List<String> cellBarcodes, final List<String> geneNames) {
		this(cellBarcodes, geneNames, CRSMatrix.zero(geneNames.size(), cellBarcodes.size()));
	}

	private DGEMatrix(final List<String> cellBarcodes, final List<String> geneNames, final Matrix m) {
		this.expressionMatrix = m;
		this.geneMap = listToMap(geneNames);
		this.cellBarcodeMap=listToMap(cellBarcodes);
		this.cellBarcodeCache = new ArrayList<>(cellBarcodes);
	}


	/**
	 * Convert a list to a map from the entry in the list to the index of that position
	 * This avoids indexOf calls which are O(N), and yields O(1).
	 * @param list a list of objects
	 * @return a map with each object, and its position in the list.
	 */
	private <T> Map<T, Integer> listToMap (final List<T> list) {
		Map<T, Integer> result = new HashMap<>();
		for (int i=0; i<list.size(); i++)
			result.put(list.get(i), i);
		return (result);
	}

	Matrix getMatrix () {
		return this.expressionMatrix;
	}

	/**
	 * Apply a change in place function to the data in this matrix.
	 */
	public void applyTransform (final MatrixTransformI function) {
		function.apply(this.expressionMatrix);
	}

	/**
	 * Convert the storage format to a dense matrix.
	 * Useful if you started sparse, but now have (or will soon have) non-zero values in at least 1/2 of the matrix elements
	 */
	public void toDenseMatrix () {
		this.expressionMatrix=this.expressionMatrix.toDenseMatrix();
	}

	/**
	 * Convert the storage format to a dense matrix.
	 * Useful if you started dense, but now have (or will soon have) zero values in at least 1/2 of the matrix elements
	 */
	public void toSparseMatrix() {
		this.expressionMatrix.toSparseMatrix();
	}

	/**
	 * Reads from the map and constructs the ordered list of Strings based on the index positions saved in the map.
	 */
	private List<String> mapToList (final Map<String, Integer> map) {
		int size=map.size();
		String [] result = new String [size];
		// get the keys out by order.
		for (String s: map.keySet()) {
			int pos=map.get(s);
			result[pos]=s;
		}
		return (Arrays.asList(result));
	}
	/**
	 *
	 * @return a list of cell barcodes in this data set.
	 */
	public List<String> getCellBarcodes () {
		if (cellBarcodeCache == null) {
			cellBarcodeCache = mapToList(this.cellBarcodeMap);
		}
		return Collections.unmodifiableList(cellBarcodeCache);
	}

	/**
	 * Removes cell barcodes (columns) from the data.
	 * This removes both the cell barcodes as well as the expression data for that cell.
	 * All remaining cells stay in their original order, but the data matrix shrinks around the removed data.
	 * @param cellBarcodes A list of cell barcodes to remove from all genes.
	 */
	public void removeCellBarcodes (final Collection <String> cellBarcodes) {
		// get the indexes where data will be removed.
		List<Integer> toRemove = new ArrayList<>();
		for (String cell: cellBarcodes) {
			Integer idx = this.cellBarcodeMap.get(cell);
			if (idx!=null) 
				toRemove.add(idx);							
		}
		// remove the indexes from the map
		for (String cell: cellBarcodes) {
			removeFromMap(cell, this.cellBarcodeMap);
		}
		cellBarcodeCache = null;
		
		Collections.sort(toRemove);
		Collections.reverse(toRemove);
		// if I remove elements from the end and work forward, I will not mess up if columns are filling in the gaps behind me.
		for (Integer idx: toRemove)
			this.expressionMatrix=this.expressionMatrix.removeColumn(idx);
	}

	/**
	 * Only retain cells that have greater than <numGenes> genes with positive expression
	 * @param numGenes The Cells must more than this number of genes with  to be retained.
	 */
	public void removeCellsWithLowExpression(final int numGenes) {
		if (numGenes==0) return;  // short circuit.
		VectorAccumulator nonZeroAccumulator = nonZeroAccumulator();
		double [] nonZeroExpression = this.expressionMatrix.foldColumns(nonZeroAccumulator);
		List<String> cells = this.getCellBarcodes();
		List<String> cellsToRemove = new ArrayList<>();
		for (int i=0; i<nonZeroExpression.length; i++)
			if (nonZeroExpression[i]<=numGenes)
				cellsToRemove.add(cells.get(i));
		this.removeCellBarcodes(cellsToRemove);
	}

	/**
	 *
	 * @param numCells Only retain genes that have greater than this many cells with positive expression
	 */
	public void removeGenesWithLowExpression(final int numCells) {
		if (numCells==0) return; // short circuit.
		VectorAccumulator nonZeroAccumulator = nonZeroAccumulator();
		double [] nonZeroExpression = this.expressionMatrix.foldRows(nonZeroAccumulator);
		List<String> genes = this.getGenes();
		List<String> genesToRemove = new ArrayList<>();
		for (int i=0; i<nonZeroExpression.length; i++)
			if (nonZeroExpression[i]<=numCells)
				genesToRemove.add(genes.get(i));
		this.removeGenes(genesToRemove);
	}


	/**
	 * @return a list of genes in this data set.
	 */
	// TODO: should these results be cached?
	public List<String> getGenes () {
		return (mapToList(this.geneMap));
	}

	/**
	 * Remove genes from this data set.  Removes the rows of expression that have these genes.
	 * TODO: For LARGE data sets, this is horribly slow....
	 * @param geneNames A list of genes.
	 */
	public void removeGenes(final List<String> geneNames) {
		// get the indexes where data will be removed.
		List<Integer> toRemove = new ArrayList<>();
		for (String gene : geneNames) {
			Integer idx = this.geneMap.get(gene);
			if (idx != null)
				toRemove.add(idx);
		}

		// remove from map.
		for (String gene : geneNames)
			removeFromMap(gene, this.geneMap);

		Collections.sort(toRemove);
		Collections.reverse(toRemove);

		/*
		 *this is ridiculously, horribly slow.
		if (this.expressionMatrix instanceof CRSMatrix) {
			log.info("Converting Matrix");
			this.expressionMatrix=this.expressionMatrix.toColumnMajorSparseMatrix();
			log.info("Matrix Converted");
		}
		*/
		// if I remove elements from the end and work forward, I will not mess
		// up if columns are filling in the gaps behind me.
		int counter=0;
		for (Integer idx : toRemove) {
			counter++;
			this.expressionMatrix = this.expressionMatrix.removeRow(idx);
			// if (counter%100==0) log.info("Removed " + counter + "genes of [" + toRemove.size() +"]");
		}
	}

	/**
	 * After items are removed from the map, need to shift the indexes to point at the right indexes.
	 * For example, if you have 5 entries from 0-4 and entry 2 is removed, you'd need to shift all indexes
	 * past the removal to the left by 1.
	 */
	private void removeFromMap (final String key, final Map<String, Integer> map) {
		Integer idx = map.get(key);
		if (idx==null) return;
		map.remove(key);
		for (String k: map.keySet()) {
			Integer value = map.get(k);
			if (value>idx) {
				value--;
				map.put(k, value);
			}
		}

	}
	/**
	 * Returns expression for a gene, or null if the gene does not exist.
	 * The cells represented are in the same order as getCellBarcodes() returns.
	 * @param gene the gene name to get expression for
	 * @return An array of expression data.  If there is no gene with this name, return null.
	 */
	// TODO: evaluate if getting a vector is faster than pulling single values out.
	public double [] getExpression (final String gene) {
		if (this.expressionMatrix instanceof CRSMatrix)
			return (getExpressionSparse(gene));
		else
			return (getExpressionDense(gene));
	}

	private double [] getExpressionSparse (final String gene) {
		CRSMatrix thisM = (CRSMatrix) this.expressionMatrix;
		Integer rowIdx = this.geneMap.get(gene);
		if (rowIdx==null) return null;
		double [] result = new double [this.cellBarcodeMap.size()];
		for (int i=0; i<result.length; i++)
			result[i] = thisM.getOrElse(rowIdx, i, this.SPARSE_VALUE);
		return result;
	}

	private double [] getExpressionDense (final String gene) {
		Integer rowIdx = this.geneMap.get(gene);
		if (rowIdx==null) return null;
		double [] result = new double [this.cellBarcodeMap.size()];
		for (int i=0; i<result.length; i++)
			result[i] = this.expressionMatrix.get(rowIdx, i);
		return result;
	}

	/**
	 * Convenience method
	 * If the input DGE row is not null, return it's expression array
	 * If it is null, return a new array with numElements length filled with the missing value.
	 */
	private double [] getExpressionWithNull (final double [] r, final float missingValue, final int numElements) {
		if (r!=null) return r;
		double [] result = new double [numElements];
		Arrays.fill(result, missingValue);
		return (result);
	}

	/**
	 * Get a 2d matrix of expression for this experiment.
	 * Rows are genes (use getGenes to get row names)
	 * Columns are cells (use getCellBarcodes to get cell barcode ordering)
	 * @return a 2d float matrix of expression data.  This is a copy of the original data.
	 */
	public double [] [] getExpressionMatrix () {
		return this.expressionMatrix.toDenseMatrix().toArray();
	}

	public DGEMatrix merge (final DGEMatrix other) {
		if (this.getCellBarcodes().equals(other.getCellBarcodes())) {
			return mergeExpressionForCells(other);
		} else {
			return mergeDisjointCells(other);
		}
	}

	private DGEMatrix mergeExpressionForCells(DGEMatrix other) {
		List<String> thisGenes = this.getGenes();
		List<String> otherGenes = other.getGenes();
		Collection<String> geneOverlap = CollectionUtils.intersection(thisGenes, otherGenes);
		if (!geneOverlap.isEmpty()) {
			throw new IllegalArgumentException("The two matrixes have overlapping genes, this should not happen:" +
					geneOverlap.toString());
		}
		List<String> allGenes = new ArrayList<> (thisGenes);
		CollectionUtils.addAll(allGenes, otherGenes);
		Collections.sort(allGenes);
		CRSMatrix m = new CRSMatrix(allGenes.size(), expressionMatrix.columns());
		final ProgressLogger progressLogger = new ProgressLogger(log, 1000, "processed",
				String.format("genes of %d", allGenes.size()));
		final Set<String> thisGenesSet = new HashSet<>(thisGenes);
		for (int i = 0; i < allGenes.size(); ++i) {
			final String gene = allGenes.get(i);
			final DGEMatrix sourceMatrix = thisGenesSet.contains(gene)? this: other;
			m.setRow(i, sourceMatrix.expressionMatrix.getRow(sourceMatrix.geneMap.get(gene)));
			progressLogger.record(gene, i);
		}
		return (new DGEMatrix(this.getCellBarcodes(), allGenes, m));
	}

	/**
	 * Merges two DGE matrixes together.  This retains all columns in both experiments, and adds rows that are missing in either data set.
	 * The missing values are set to <missingValue>.
	 * If the two experiments have overlapping cell barcodes, throw an IllegalArgumentException.
	 * @param other The other data set to merge with this one.
	 * @return a new DGEMatrix containing the merge of the this with other.
	 */
	public DGEMatrix mergeDisjointCells (final DGEMatrix other) {
		Collection<String>overlap = CollectionUtils.intersection(this.getCellBarcodes(), other.getCellBarcodes());
		if (!overlap.isEmpty()) {
			throw new IllegalArgumentException("The two matrixes have overlapping cell barcodes, this should not happen:" + overlap.toString());
		}

		// merge the cell barcodes.
		List <String> cellBarcodes = new ArrayList<>(this.getCellBarcodes());
		cellBarcodes.addAll(other.getCellBarcodes());
		List<String> allGenes = new ArrayList<> (CollectionUtils.union(this.getGenes(), other.getGenes()));
		Collections.sort(allGenes);

		int numCellsThis=this.getCellBarcodes().size();
		int numCellsOther=other.getCellBarcodes().size();
		CRSMatrix m = new CRSMatrix(allGenes.size(), cellBarcodes.size());

		for (int rowIdx=0; rowIdx<allGenes.size(); rowIdx++) {
			String g = allGenes.get(rowIdx);
			double [] thisExp = getExpressionWithNull(this.getExpression(g), 0, numCellsThis);
			double [] otherExp = getExpressionWithNull(other.getExpression(g), 0, numCellsOther);
			double [] expression  = ArrayUtils.addAll(thisExp, otherExp);
			for (int columnIdx=0; columnIdx<expression.length; columnIdx++)
				if (expression[columnIdx]!=0)
					m.set(rowIdx, columnIdx, expression[columnIdx]);
		}
		return (new DGEMatrix(cellBarcodes, allGenes, m));
	}

	/**
	 * Merges two DGE matrixes together.  This retains all cell barcodes present in both experiments, and adds rows that are missing in either data set.
	 * If cell barcodes are present in both experiments, their UMI counts per gene [rows of the matrix] are summed.
	 * The missing values are set to <missingValue>.
	 * @param other The other data set to merge with this one.
	 * @return a new DGEMatrix containing the merge of the this with other.
	 */
	public DGEMatrix mergeWithCollapse (final DGEMatrix other) {
		// start with this data set's cell barcodes, find barcodes that are new to the 2nd matrix, and add those.
		List<String> cellBarcodes = new ArrayList<>(this.getCellBarcodes());
		// get the other list of cell barcodes, remove cells from the first list.
		List<String> cellBarcodesOther = new ArrayList<>(other.getCellBarcodes());
		cellBarcodesOther.removeAll(cellBarcodes);
		//merge results together
		cellBarcodes.addAll(cellBarcodesOther);

		// merge the genes.
		List<String> allGenes = new ArrayList<> (CollectionUtils.union(this.getGenes(), other.getGenes()));
		Collections.sort(allGenes);

		int numCellsThis=this.getCellBarcodes().size();
		int numCellsOther=other.getCellBarcodes().size();

		// for a donor, get their position in each array [can be null in one array but not both].
		// return the expression value of the donor.
		// need a map from the donor name to the position in the array.

		CRSMatrix m = new CRSMatrix(allGenes.size(), cellBarcodes.size());

		for (int rowIdx=0; rowIdx<allGenes.size(); rowIdx++) {
			String g = allGenes.get(rowIdx);
			double [] thisExp = getExpressionWithNull(this.getExpression(g), 0, numCellsThis);
			double [] otherExp = getExpressionWithNull(other.getExpression(g), 0, numCellsOther);
			double [] expression  = getExpression(cellBarcodes, thisExp, otherExp, this.cellBarcodeMap, other.cellBarcodeMap);
			for (int columnIdx=0; columnIdx<expression.length; columnIdx++)
				if (expression[columnIdx]!=0)
					m.set(rowIdx, columnIdx, expression[columnIdx]);
		}
		return (new DGEMatrix(cellBarcodes, allGenes, m));
	}

	private double [] getExpression (final List<String> cellBarcodes, final double [] expressionOne, final double [] expressionTwo, final Map<String,Integer> mapOne, final Map<String,Integer> mapTwo) {
		double [] expression = new double [cellBarcodes.size()];
		for (int i=0; i<cellBarcodes.size(); i++)
			expression[i]=getExpression(cellBarcodes.get(i), expressionOne, expressionTwo, mapOne, mapTwo);
		return expression;
	}
	private double getExpression (final String cellBarcode, final double [] expressionOne, final double [] expressionTwo, final Map<String,Integer> mapOne, final Map<String,Integer> mapTwo) {
		Integer posOne = mapOne.get(cellBarcode);
		Integer posTwo = mapTwo.get(cellBarcode);
		if (posOne!=null && posTwo==null) return expressionOne[posOne];
		if (posOne==null && posTwo!=null) return expressionTwo[posTwo];
		if (posOne==null && posTwo==null)
			// this should never happen.
			throw new IllegalArgumentException("Requested a cell barcode [" + cellBarcode +"] not present in either data set.");
		return expressionOne[posOne]+expressionTwo[posTwo];
	}

	/**
	 * Sum the columns of the matrix and return one entry per sample.
	 * @return The total expression of each sample.
	 */
	public double [] getTotalExpressionPerSample () {
		final VectorAccumulator sumAccum = Vectors.asSumAccumulator(0);
		final double [] sums = expressionMatrix.foldColumns(sumAccum);
		return sums;
	}

	/**
	 * Parses a DigitalGeneExpression file.  If there are no lines in the file, return an empty DGEMatrix with no rows and no cells.
	 * The format of the DGE file is tab delimited and the rows/columns are as follows:
	 * 1st line header, first column fixed GENE header followed by cell barcodes.
	 * all other lines: first column has the gene name, other columns have the expression for that gene.
	 * @param input The input file to parse
	 * @param cellBarcodePrefix a prefix to add to the cell barcodes.  Useful to distinguish these barcodes from other barcodes in a different DGEExperiment.  If set to null
	 * no prefix is added.
	 */
	public static DGEMatrix parseFile (final File input, final String cellBarcodePrefix) {
        final FileFormat format = detectFileFormat(input);

        if (format == FileFormat.DENSE)
			return parseDenseFile(input, cellBarcodePrefix);
		else
			return parseDropSeqMatrixMarket(input, cellBarcodePrefix);
    }

	public static DgeHeader parseDgeHeader(final File input) {
        BufferedReader reader = null;
        try {
            reader = IOUtil.openFileForBufferedReading(input);
            return new DgeHeaderCodec().decode(reader, input.getAbsolutePath());
        } finally {
            CloserUtil.close(reader);
        }
    }

	public static DGEMatrix parseFile (final File input, final File genes, final File cells, final String cellBarcodePrefix) {
		final FileFormat format = detectFileFormat(input);
		if (format == FileFormat.DENSE && cells==null && genes==null)
			return parseDenseFile(input, cellBarcodePrefix);
		if (format == FileFormat.MM_SPARSE && cells==null && genes==null)
			return parseDropSeqMatrixMarket(input, cellBarcodePrefix);
		if (format == FileFormat.MM_SPARSE_10X && cells!=null && genes!=null)
			return parse10XGenomicsMatrixMarket(input, genes, cells, cellBarcodePrefix);
		else
			throw new IllegalArgumentException("Cell and Gene files passed in, but doesn't look like 10x Genomics data.  Not sure what to do!");
	}

	public static DGEMatrix parseDenseFile(final File input, final String cellBarcodePrefix) {
		return parseDenseFile(input, cellBarcodePrefix, true);
	}

    public static DGEMatrix parseDenseFile(final File input, final String cellBarcodePrefix, final boolean verbose) {
        // get lines in the file to determine the row dimensions.
        if (verbose) log.info("Getting lines in input file");
        int lines = DGEMatrix.countLines(input);

        // The reading of the header here is only to remove the header from the input stream.
        final BufferedInputStream inputStream = new BufferedInputStream(IOUtil.openFileForReading(input));
        new DgeHeaderCodec().decode(inputStream, input.getAbsolutePath());


        TabbedInputParser parser = new TabbedInputParser(false, inputStream);
        if (!parser.hasNext())  {
            parser.close();
            return new DGEMatrix(new ArrayList<String>(), new ArrayList<String>(), CRSMatrix.zero(0, 0));
        }

        String [] header = parser.next();
        int numCells= header.length-1;

        ArrayList<String> cellBarcodes = new ArrayList<> (numCells);
        for (int i=1; i<header.length; i++)
            if (cellBarcodePrefix!=null)
                cellBarcodes.add(cellBarcodePrefix + header[i]);
            else
                cellBarcodes.add(header[i]);
        if (verbose) log.info("Found [" + (lines-1) + "] genes and [" + cellBarcodes.size() +"] cells");

        // initialize the sparse matrix
        CRSMatrix m = CRSMatrix.zero(lines-1, cellBarcodes.size());

        List<String> geneNames = new ArrayList<>();

        int rowIdx=0;
        while(parser.hasNext()) {
            String [] line =parser.next();
            String gene = line[0];
            geneNames.add(gene);
            for (int columnIdx=1; columnIdx<line.length; columnIdx++) {
                double expression=Double.parseDouble(line[columnIdx]);
                if (expression!=0)
                    m.set(rowIdx, columnIdx-1, expression);
            }
            rowIdx++;
            if (rowIdx%1000==0)
            	if (verbose) log.info("Parsed [" + rowIdx +"] lines of DGE File [" + input.getName() +"]");
        }

        parser.close();
		return (new DGEMatrix(cellBarcodes, geneNames, m));
    }

    /**
     * Parse a 10x genomics matrix market format DGE file.
     * @param input The matrix of expression values, which are integers
     * @param genesFile A file containing a list of gene names in the 2nd column corresponding to rows of the matrix market data
     * @param barcodesFile A file containing a list of cell barcodes in the 1st column corresponding to columns of the matrix market data
     * @param cellBarcodePrefix A prefix to add at the front of all cell barcodes.
     */
    public static DGEMatrix parse10XGenomicsMatrixMarket(final File input, final File genesFile, final File barcodesFile, final String cellBarcodePrefix) {
		final MatrixMarketReader matrixReader = new MatrixMarketReader(input);
		if (matrixReader.getElementType() != MatrixMarketConstants.ElementType.integer)
			throw new RuntimeException(input.getAbsolutePath() + " does not start with correct 10x Genomics Matrix Market header");
		//TODO: this should be the 2nd column [1] to get the gene symbol.
		List<String> geneNames=parse10xGenomicsAttributeFile(genesFile, 0, 2);
		List<String> cellBarcodes=parse10xGenomicsAttributeFile(barcodesFile, 0, 1);
		return (parseGenericMatrixMarket(matrixReader, cellBarcodePrefix, geneNames, cellBarcodes));

	}

    /**
     *
     * @param input The file containing either the genes or the cell barcodes.  Neither has headers.
     * @param columnIdx the 0 based column index to extract from the file.
     * @return list of values
     */
    private static List<String> parse10xGenomicsAttributeFile (final File input, final int columnIdx, final int numColumnsInFile) {
    	IOUtil.assertFileIsReadable(input);
		List<String> result = new ArrayList<>();
		BasicInputParser parser = new BasicInputParser(false, numColumnsInFile, input);
		while(parser.hasNext()) {
			String [] line =parser.next();
			result.add(line[columnIdx]);
		}
		parser.close();
		return (result);
    }



    public static DGEMatrix parseDropSeqMatrixMarket(final File input, final String cellBarcodePrefix) {
		final MatrixMarketReader mmReader = new MatrixMarketReader(input, MatrixMarketConstants.GENES, MatrixMarketConstants.CELL_BARCODES);
		return (parseGenericMatrixMarket(
				mmReader, cellBarcodePrefix,
				mmReader.getRowNames(),
				new ArrayList<>(mmReader.getColNames())));
	}

    private static DGEMatrix parseGenericMatrixMarket (final MatrixMarketReader matrixReader, final String cellBarcodePrefix,
                                                       final List<String> geneNames, final List<String> cellBarcodes) {
        if (cellBarcodePrefix != null && !cellBarcodePrefix.isEmpty())
            for (int i = 0; i < cellBarcodes.size(); ++i)
                cellBarcodes.set(i, cellBarcodePrefix + cellBarcodes.get(i));
        int rows = matrixReader.getNumRows();
        int cols = matrixReader.getNumCols();
        if (rows != geneNames.size())
            throw new RuntimeException("Number of rows in matrix does not agree with length of gene list in " + matrixReader.getFilename());
        if (cols != cellBarcodes.size())
            throw new RuntimeException("Number of columns in matrix does not agree with length of cell barcode list in " + matrixReader.getFilename());
        log.info("Found [" + rows + "] genes and [" + cols +"] cells");

        // initialize the sparse matrix
        CRSMatrix m = CRSMatrix.zero(rows, cols);
        for (final MatrixMarketReader.Element element: matrixReader)
			m.set(element.row, element.col, element.realValue());

        CloserUtil.close(matrixReader);
		return (new DGEMatrix(cellBarcodes, geneNames, m));

    }

    /**
	 * Parses a DGE file.  If there are no lines in the file, return an empty {@link DGEMatrix} with no rows and no cells.
	 * The format of the DGE file is tab delimited and the rows/columns are as follows:
	 * 1st line header, first column fixed GENE header followed by cell barcodes.
	 * all other lines: first column has the gene name, other columns have the expression for that gene.
	 * Assumes no added prefix name for cell IDs, and the sparse value is 0.
	 * @param input The input file to parse
	 */
	public static DGEMatrix parseFile (final File input) {
		return parseFile(input, null);
	}

	private static int countLines(final File input) {
        BasicInputParser bip = null;
        try {
            int cnt = 0;
            final BufferedInputStream inputStream = new BufferedInputStream(IOUtil.openFileForReading(input));
            // Skip header if present
            new DgeHeaderCodec().decode(inputStream, input.getAbsolutePath());
            bip = new BasicInputParser(true, inputStream);
            while (bip.hasNext())
				bip.next();
            cnt= bip.getCurrentLineNumber();
            return cnt;
        } finally {
            CloserUtil.close(bip);
        }
    }

	private static FileFormat detectFileFormat(final File input) {
        BufferedReader reader = null;
        try {
            reader = IOUtil.openFileForBufferedReading(input);
            if (MatrixMarketReader.isMatrixMarketReal(reader))
				return FileFormat.MM_SPARSE;
			else if (MatrixMarketReader.isMatrixMarketInteger(reader))
				return FileFormat.MM_SPARSE_10X;
			else
				return FileFormat.DENSE;
        } finally {
            CloserUtil.close(reader);
        }
    }

	/**
	 * Write this object out in DigitalGeneExpression format.
	 * @param output The output file to write to
	 * @param formatAsInteger Should the DGE expression values be truncated to integer values?
	 */
    public void writeFile (final File output, final boolean formatAsInteger, final FileFormat format) {
        writeFile(output, formatAsInteger, format, null);
    }
    public void writeFile (final File output, final boolean formatAsInteger, final FileFormat format, final DgeHeader header) {
		if (format == FileFormat.DENSE)
			writeDenseDgeFile(output, formatAsInteger, header);
		else if (format == FileFormat.MM_SPARSE) {
            if (header != null)
				throw new IllegalArgumentException("DGE header not support for sparse matrix formats");
            writeDropSeqMatrixMarket(output, formatAsInteger, false);
        } else
			throw new IllegalArgumentException(format.name() + " Not yet supported for output.");
	}

	public void writeDenseDgeFile(final File output, final boolean formatAsInteger) {
        writeDenseDgeFile(output, formatAsInteger, null);
    }

    public void writeDenseDgeFile(final File output, final boolean formatAsInteger, final DgeHeader dgeHeader) {
        try {
			IOUtil.assertFileIsWritable(output);
			List<String> cellBarcodes = this.getCellBarcodes();
			List<String> genes = this.getGenes();
			BufferedWriter out = new BufferedWriter(new OutputStreamWriter(IOUtil.openFileForWriting(output)));

            if (dgeHeader != null)
				new DgeHeaderCodec().encode(out, dgeHeader);

			// write header
			List<String> header = new ArrayList<>(cellBarcodes.size()+1);
			header.add("GENE");
			for (String c: cellBarcodes)
                header.add(c);
			String h = StringUtils.join(header, "\t");
			out.write(h);
			out.newLine();

			// write body
			for (String gene: genes) {
                List<String> line = new ArrayList<>(cellBarcodes.size()+1);
                line.add(gene);
                double [] expressionByGene = this.getExpression(gene);
                for (double exp: expressionByGene)
					line.add(formatExpressionValue(exp, formatAsInteger));
                String b = StringUtils.join(line, "\t");
                out.write(b); out.newLine();
            }
			out.close();
		} catch (IOException e) {
			throw new RuntimeException("Trouble writing " + output.getAbsolutePath(), e);
		}
	}

	private String formatExpressionValue(final double exp, final boolean formatAsInteger) {
		if (formatAsInteger)
			return Integer.toString((int)exp);
		else
			return Double.toString(exp);
	}

	public void writeDropSeqMatrixMarket(final File output, final boolean formatAsInteger, final boolean transpose) {
		try {
			IOUtil.assertFileIsWritable(output);
			final int cardinality;
			if (expressionMatrix instanceof SparseMatrix)
				cardinality = ((SparseMatrix) expressionMatrix).cardinality();
			else
				cardinality = expressionMatrix.columns() * expressionMatrix.columns();
			final MatrixMarketWriter writer = new MatrixMarketWriter(output, MatrixMarketConstants.ElementType.real,
					transpose? expressionMatrix.columns(): expressionMatrix.rows(),
					transpose? expressionMatrix.rows(): expressionMatrix.columns(),
					cardinality, this.getGenes(),this.getCellBarcodes(),
					MatrixMarketConstants.GENES, MatrixMarketConstants.CELL_BARCODES);
			if (expressionMatrix instanceof SparseMatrix) {
                final SparseMatrix mat = (SparseMatrix)expressionMatrix;
                final MatrixIterator it = mat.nonZeroIterator();
                while (it.hasNext()) {
					final double val = it.next();
                    final int row = it.rowIndex();
                    final int col = it.columnIndex();
					writeMatrixMarketTriplet(writer, row, col, val, formatAsInteger, transpose);
                }
            } else
				for (int i = 0; i < expressionMatrix.rows(); ++i)
					for (int j = 0; j < expressionMatrix.columns(); ++j) {
						double exp = expressionMatrix.get(i, j);
						if (exp != SPARSE_VALUE)
							writeMatrixMarketTriplet(writer, i, j, exp, formatAsInteger, transpose);
					}
			writer.close();
		} catch (IOException e) {
			throw new RuntimeException("Trouble writing " + output.getAbsolutePath(), e);
		}
	}

	private void writeMatrixMarketTriplet(final MatrixMarketWriter writer, int row, int col, final double val,
										  final boolean formatAsInteger, final boolean transpose) {
    	if (transpose) {
    		int swap = row;
    		row = col;
    		col = swap;
		}
		if (formatAsInteger)
			writer.writeTriplet(row, col, (int)Math.round(val));
		else
			writer.writeTriplet(row, col, val);
	}

	private String formatMatrixMarketLine(int row, int col, final double exp, final boolean formatAsInteger, final boolean transpose) {
		if (transpose) {
			final int swap = row;
			row = col;
			col = swap;
		}
		return String.format("%d %d %s", row+1, col+1, formatExpressionValue(exp, formatAsInteger));
	}

	private String formatFirstMatrixMarketLine(final long cardinality, final boolean transpose) {
		final int nrows;
		final int ncols;
		if (transpose) {
			nrows = expressionMatrix.columns();
			ncols = expressionMatrix.rows();
		} else {
			nrows = expressionMatrix.rows();
			ncols = expressionMatrix.columns();
		}
		return String.format("%d %d %d", nrows, ncols, cardinality);
	}

	private static VectorAccumulator nonZeroAccumulator() {
        return new VectorAccumulator() {

            private int count=0;
            @Override
            public void update(final int i, final double value) {
            	if (value!=0)
					count++;
            }

            @Override
            public double accumulate() {
                double result = count;
                // if you don't reset the variables, further calls to this accumulator will gain the old results
                count=0;
                return result;
            }
        };
    }


}
