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

import java.math.BigDecimal;
import java.math.RoundingMode;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools.MatrixTransformI;
import org.la4j.Matrix;
import org.la4j.Vectors;
import org.la4j.vector.functor.VectorAccumulator;
import org.la4j.vector.functor.VectorFunction;

/**
 * Produces functions that transform a matrix.
 * @author nemesh
 *
 */
public class MatrixTransformFactory {

	/**
	 * For a given matrix, calculate the sum of each column, then divide each
	 * column by it's sum, so each column adds to 1.
	 * In R this would be (assuming a is your matrix of DGE data):
	 * b=sweep(a,2,colSums(a),'/')
	 * This function keeps sparse matrixes sparse.
	 * @param m
	 */
	public static MatrixTransformI normalizeColumns () {
		return new MatrixTransformI() {

			@Override
			public void apply(final Matrix m) {
				final VectorAccumulator sumAccum = Vectors.asSumAccumulator(0);
				final double [] sums = m.foldColumns(sumAccum);
				// divide by the sums.  This is like an apply to each column.
				for (int i=0; i<sums.length; i++) {
					// make a function with the set value to divide by.
					VectorFunction f = Vectors.asDivFunction(sums[i]);
					m.updateColumn(i, f);
				}
			}
		};

	}

	/**
	 * For a given matrix, apply a transform to every cell of the matrix
	 * For example, take the log of (10000*N)+1 where N is the value of the cell in the matrix.
	 * In R this would be equivalent to log(10000*sweep(a,2,colSums(a),'/')+1)
	 * This function keeps sparse matrixes sparse.
	 * @param multiply Multiply by this value
	 * @param add add this value
	 */
	public static MatrixTransformI logOfDGE (final int multiply, final int add) {
		 return new MatrixTransformI() {
			 private final int mult = multiply;
			 private final int a = add;

			 @Override
			 public void apply(final Matrix m) {
				 // Note, I tried this as a MatrixFunction and called m.update(function) on it, but it gave me the wrong results.
				 // Lots of elements were 0?
				 // apply by row/column seems to work ok though.
				 VectorFunction f = logTransform(multiply, add);
				 int numRows = m.rows();
				 for (int i=0; i<numRows; i++)
					m.updateRow(i, f);
			 }
		 };
	}


	/**
	 * For a given matrix, scale the data by rows.
	 * For each row, calculate the mean.  Subtract this from the row.
	 * For each row, calculate the SD.  Divide the row by the SD.
	 * In R this would be equivalent to t(scale (t(x))) where x is your data matrix.
	 * This function does NOT keeps sparse matrixes sparse.  Data size explosion here.
	 * Use scaleDataNonZero to avoid this issue.
	 */
	public static MatrixTransformI scaleData () {
		return scaleData(true);
	}

	/**
	 * For a given matrix, scale the data by rows.
	 * For each row, calculate the mean.  Subtract this from the row.
	 * For each row, calculate the SD.  Divide the row by the SD.
	 * In R this would be equivalent to t(scale (t(x))) where x is your data matrix.
	 * This function keeps sparse matrixes sparse, by only processing non-zero elements of the matrix.
	 * This is experimental - we're not sure if this normalization will be equivalent for downstream analysis.
	 *
	 */
	// TODO: this still under construction.
	// TODO: I think this actually may not make logical sense.
	// If you scale low expression to be lower than 0 for genes that have expression values >0, then
	// you're making the missing data have a higher expression than some of the genes.

	private static MatrixTransformI scaleDataNonZero () {
		return scaleData(false);
	}

	// does the work of the two public methods.
	// retainZeros=false needs more testing.
	private static MatrixTransformI scaleData (final boolean retainZeros) {
		return new MatrixTransformI() {

			@Override
			public void apply(final Matrix m) {
				// this forces the data to be dense instead of sparse.
				// this doesn't work somehow, this m is being changed but the changes aren't occurring to the m that was passed in.
				// side effect not working?
				// m=m.toDenseMatrix();
				// m.toDenseMatrix(); // this returns a dense matrix, but we don't capture that.
				// calculate means
				VectorAccumulator meanAccumulator = meanAccumulator(retainZeros);
				double [] means = m.foldRows(meanAccumulator);
				// subtract means from rows.
				for (int i=0; i<means.length; i++) {
					// make a function with the set value to divide by.
					VectorFunction f = Vectors.asMinusFunction(means[i]);
					m.updateRow(i, f);
				}
				// calculate SD
				VectorAccumulator sdAccumulator = sdAccumulator(retainZeros);
				double [] sd = m.foldRows(sdAccumulator);
				for (int i=0; i<sd.length; i++) {
					// make a function with the set value to divide by.
					VectorFunction f = Vectors.asDivFunction(sd[i]);
					m.updateRow(i, f);
				}
			}
		};
	}

	/**
	 * A transform that multiplies a number by a value, adds some other number, then takes the log
	 *
	 * @param multiply Multiply by this value
	 * @param add add this value
	 * @return
	 */
	private static VectorFunction logTransform(final int multiply, final int add) {
        return new VectorFunction() {
            private final int m = multiply;
            private final int a = add;

            @Override
            public double evaluate(final int i, final double value) {
                double vv = (value*m)+a;
                double vvv=Math.log(vv);
                return vvv;
            }
        };
    }

	private static VectorAccumulator meanAccumulator(final boolean retainZeros) {
        return new VectorAccumulator() {
            private BigDecimal result = new BigDecimal(0);
            private int count=0;
            @Override
            public void update(final int i, final double value) {
            	if (retainZeros || value!=0) {
            		result = result.add(new BigDecimal(value));
            		count++;
            	}
            }

            @Override
            public double accumulate() {
                double value = result.setScale(Vectors.ROUND_FACTOR, RoundingMode.CEILING).doubleValue();
                double mean = value / count;
                // if you don't reset the variables, further calls to this accumulator will gain the old results
                count=0;
                result = new BigDecimal(0);
                return mean;
            }
        };
    }

	private static VectorAccumulator sdAccumulator(final boolean retainZeros) {
        return new VectorAccumulator() {
            StandardDeviation d = new StandardDeviation();

            @Override
            public void update(final int i, final double value) {
            	if (retainZeros || value!=0)
					d.increment(value);

            }

            @Override
            public double accumulate() {
                double sd = d.getResult();
                d=new StandardDeviation();
                return sd;
            }
        };
    }



}
