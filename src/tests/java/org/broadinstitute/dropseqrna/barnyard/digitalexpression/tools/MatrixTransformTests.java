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

import org.la4j.Matrix;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class MatrixTransformTests {

	private final File exampleOne = new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalexpression/dge_example1.txt.gz");
	private final File exampleTwo = new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalexpression/dge_example2.txt.gz");
	private final File preMergedExample = new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalexpression/dge_example_merged.txt.gz");
	// a=read.table("/Users/nemesh/dropseqrna/transcriptome_java/testdata/org/broadinstitute/dropseq/private/barnyard/digitalexpression/dge_example1.txt.gz", header=T, stringsAsFactors=F, row.names=1)
	@Test
	public void testNormalizeRows() {
		DGEMatrix result= DGEMatrix.parseFile(exampleOne);				
		result.applyTransform(MatrixTransformFactory.normalizeColumns());
		Matrix m = result.getMatrix();
		Assert.assertNotNull(m);
		double [] [] actualValues = m.toDenseMatrix().toArray();
		// sweep(a,2,colSums(a),'/')
		double [] [] expectedVals = {{0.609756, 0.731707, 0.765957, 0.666667, 0.828571},
									   {0.024390, 0.000000, 0.000000, 0.000000, 0.000000},
									   {0.195122, 0.146341, 0.085106, 0.190476, 0.114286},
									   {0.000000, 0.000000, 0.000000, 0.000000, 0.028571},
									   {0.121951, 0.000000, 0.000000, 0.000000, 0.000000},
									   {0.024390, 0.073170, 0.063829, 0.095238, 0.028571},
									   {0.024390, 0.048780, 0.021276, 0.047619, 0.000000},
									   {0.000000, 0.000000, 0.063829, 0.000000, 0.000000}};
		
		for (int i=0; i<expectedVals.length; i++) {
			for (int j=0; j<expectedVals[0].length; j++) {
				Assert.assertEquals(expectedVals[i][j], actualValues[i][j], 0.0001);
			}
		}
						
		/**
		 *  		CAATCCGACAAC CACTAAAGCCAG TCCCTTCAAGTA ATGGTCTCAAAC CCTTCCATGCGA
			Arl6ip1   0.60975610   0.73170732   0.76595745   0.66666667   0.82857143
			Prkg1     0.02439024   0.00000000   0.00000000   0.00000000   0.00000000
			Rasa1     0.19512195   0.14634146   0.08510638   0.19047619   0.11428571
			Exosc2    0.00000000   0.00000000   0.00000000   0.00000000   0.02857143
			Gimap5    0.12195122   0.00000000   0.00000000   0.00000000   0.00000000
			Ranbp6    0.02439024   0.07317073   0.06382979   0.09523810   0.02857143
			Tet3      0.02439024   0.04878049   0.02127660   0.04761905   0.00000000
			Nmnat3    0.00000000   0.00000000   0.06382979   0.00000000   0.00000000
		 */
	}
	
	/**
	 * The equivalent of sweep(a,2,colSums(a),'/') followed by log ((10000*(N))+1)
	 */
	@Test(enabled=true)
	public void testNormalizeLogRows() {
		DGEMatrix result= DGEMatrix.parseFile(exampleOne);
		
		result.applyTransform(MatrixTransformFactory.normalizeColumns());
		result.applyTransform(MatrixTransformFactory.logOfDGE(10000, 1));
		Matrix m = result.getMatrix();
		Assert.assertNotNull(m);
		double [] [] actualValues = m.toDenseMatrix().toArray();
		// b=log(10000*sweep(a,2,colSums(a),'/')+1)
		double [] [] expectedVals = {{8.715808, 8.898102, 8.943842, 8.805025, 9.022409},
									 {5.50086, 0, 0, 0, 0},
									   {7.576722, 7.289211, 6.747661, 7.552637, 7.042161},
									   {0, 0, 0, 0, 5.658486},
									   {7.107026, 0, 0, 0, 0},
									   {5.50086, 6.596746, 6.46037, 6.860015, 5.658486},
									   {5.50086, 6.191963, 5.364882, 6.167916, 0},
									   {0, 0, 6.46037, 0, 0}};
		
		for (int i=0; i<expectedVals.length; i++) {
			for (int j=0; j<expectedVals[0].length; j++) {
				Assert.assertEquals(expectedVals[i][j], actualValues[i][j], 0.0001);
			}
		}
						
		/**
		 *  		CAATCCGACAAC CACTAAAGCCAG TCCCTTCAAGTA ATGGTCTCAAAC CCTTCCATGCGA
		Arl6ip1     8.715808     8.898102     8.943842     8.805025     9.022409
		Prkg1       5.500860     0.000000     0.000000     0.000000     0.000000
		Rasa1       7.576722     7.289211     6.747661     7.552637     7.042161
		Exosc2      0.000000     0.000000     0.000000     0.000000     5.658486
		Gimap5      7.107026     0.000000     0.000000     0.000000     0.000000
		Ranbp6      5.500860     6.596746     6.460370     6.860015     5.658486
		Tet3        5.500860     6.191963     5.364882     6.167916     0.000000
		Nmnat3      0.000000     0.000000     6.460370     0.000000     0.000000
		 */
	}
	
	/**
	 * The equivalent of sweep(a,2,colSums(a),'/') followed by log ((10000*(N))+1), followed by t(scale (t(b)))
	 */
	@Test(enabled=true)
	public void test() {
		DGEMatrix result= DGEMatrix.parseFile(exampleOne);
		Matrix m = result.getMatrix();
		result.applyTransform(MatrixTransformFactory.normalizeColumns());
		result.applyTransform(MatrixTransformFactory.logOfDGE(10000, 1));
		result.toDenseMatrix();
		result.applyTransform(MatrixTransformFactory.scaleData());
		double [] [] actualValues  = result.getExpressionMatrix();
		Assert.assertNotNull(m);
		double [] [] expectedVals = {{-1.34803, 0.176123, 0.558552, -0.60209, 1.215444},
				 					 {1.788854, -0.447214, -0.447214, -0.447214, -0.447214},
				 					 {0.952354, 0.135109, -1.404233, 0.883893, -0.567124},
				 					 {-0.447214, -0.447214, -0.447214, -0.447214, 1.788854},
				 					 {1.788854, -0.447214, -0.447214, -0.447214, -0.447214},
				 					 {-1.19003, 0.63538, 0.40822, 1.073904, -0.927473},
				 					 {0.326133, 0.589522, 0.27431, 0.580357, -1.770322},
				 					 {-0.447214, -0.447214, 1.788854, -0.447214, -0.447214}};
		
		for (int i=0; i<expectedVals.length; i++) {
			for (int j=0; j<expectedVals[0].length; j++) {
				Assert.assertEquals(expectedVals[i][j], actualValues[i][j], 0.0001);
			}
		}
		
		/**
		 *         CAATCCGACAAC CACTAAAGCCAG TCCCTTCAAGTA ATGGTCTCAAAC CCTTCCATGCGA
			Arl6ip1   -1.3480298    0.1761233    0.5585525   -0.6020897    1.2154437
			Prkg1      1.7888544   -0.4472136   -0.4472136   -0.4472136   -0.4472136
			Rasa1      0.9523543    0.1351095   -1.4042331    0.8838930   -0.5671236
			Exosc2    -0.4472136   -0.4472136   -0.4472136   -0.4472136    1.7888544
			Gimap5     1.7888544   -0.4472136   -0.4472136   -0.4472136   -0.4472136
			Ranbp6    -1.1900304    0.6353800    0.4082197    1.0739040   -0.9274734
			Tet3       0.3261329    0.5895222    0.2743097    0.5803573   -1.7703221
			Nmnat3    -0.4472136   -0.4472136    1.7888544   -0.4472136   -0.4472136
		 */
	}
	// b=log(10000*sweep(a,2,colSums(a),'/')+1)
	// t(scale (t(b)))
	
}
