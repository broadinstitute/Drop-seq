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

import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.digitalexpression.BarcodeSimulator;
import org.testng.Assert;
import org.testng.annotations.Test;

public class MultiCellTest {

	
	private final static String rootDir="testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts";
	
	private final File doubletFile = new File(rootDir + "/multicell_2.txt");

	private final File tripletFile = new File(rootDir + "/multicell_3.txt");

	@Test
	public void testReadFileDoublet () {
		List<MultiCell> result =  MultiCell.parseFile(this.doubletFile);
		for (int i=0; i<result.size(); i++) {
			MultiCell c = result.get(i);
			if (i==0) {
				Assert.assertEquals(c.getMultiCellName(), "AAAA:BBBB");
				List<String> cells = c.getCellBarcodes();
				for (int j=0; j<cells.size(); j++) {
					if (j==0)
						Assert.assertEquals(cells.get(j), "AAAA");
					if (j==1)
						Assert.assertEquals(cells.get(j), "BBBB");

				}
			}
			if (i==1) {
				Assert.assertEquals(c.getMultiCellName(), "CCCC:DDDD");
				List<String> cells = c.getCellBarcodes();
				for (int j=0; j<cells.size(); j++) {
					if (j==0)
						Assert.assertEquals(cells.get(j), "CCCC");
					if (j==1)
						Assert.assertEquals(cells.get(j), "DDDD");

				}
			}
		}
	}

	@Test
	public void testReadFileTriplet () {
		List<MultiCell> result =  MultiCell.parseFile(this.tripletFile);
		for (int i=0; i<result.size(); i++) {
			MultiCell c = result.get(i);
			if (i==0) {
				Assert.assertEquals(c.getMultiCellName(), "AAAA:BBBB:CCCC");
				List<String> cells = c.getCellBarcodes();
				for (int j=0; j<cells.size(); j++) {
					if (j==0)
						Assert.assertEquals(cells.get(j), "AAAA");
					if (j==1)
						Assert.assertEquals(cells.get(j), "BBBB");
					if (j==2)
						Assert.assertEquals(cells.get(j), "CCCC");
				}
			}
			if (i==1) {
				Assert.assertEquals(c.getMultiCellName(), "DDDD:EEEE:FFFF");
				List<String> cells = c.getCellBarcodes();
				for (int j=0; j<cells.size(); j++) {
					if (j==0)
						Assert.assertEquals(cells.get(j), "DDDD");
					if (j==1)
						Assert.assertEquals(cells.get(j), "EEEE");
					if (j==2)
						Assert.assertEquals(cells.get(j), "FFFF");

				}
			}
			if (i==2) {
				Assert.assertEquals(c.getMultiCellName(), "GGGG:HHHH:IIII");
				List<String> cells = c.getCellBarcodes();
				for (int j=0; j<cells.size(); j++) {
					if (j==0)
						Assert.assertEquals(cells.get(j), "GGGG");
					if (j==1)
						Assert.assertEquals(cells.get(j), "HHHH");
					if (j==2)
						Assert.assertEquals(cells.get(j), "IIII");

				}
			}

		}
	}

	@Test
	public void testMakeMultiCells () {
		BarcodeSimulator s = new BarcodeSimulator(8);
		List<String> barcodes = s.getRandomBarcodes(50);

		// test asking for too many merged cells, which throws an error
		boolean foundError=false;
		try {
			MultiCell.generateMultiCells(barcodes, 30, 2);
		} catch (IllegalArgumentException i) {
			Assert.assertNotNull(i);
			foundError=true;
		}
		Assert.assertTrue(foundError);

		// try to merge cells
		int numMultiCells=10;
		int multiplicity=2;
		// assert barcodes are unique!
		Set<String> uniqueBC=new HashSet<String>(barcodes);
		// int numUniqueBCs=uniqueBC.size();

		List<MultiCell> result = MultiCell.generateMultiCells(barcodes, numMultiCells, multiplicity);
		Assert.assertSame(result.size(), numMultiCells);

		// validate that each MultiCell contains 2 cells, and none of them contain overlapping barcodes.
		for (MultiCell mc: result)
			Assert.assertSame(mc.getCellBarcodes().size(), multiplicity);

		// try to merge cells
		numMultiCells=5;
		multiplicity=3;

		result = MultiCell.generateMultiCells(barcodes, numMultiCells, multiplicity);
		Assert.assertSame(result.size(), numMultiCells);

		// validate that each MultiCell contains 3 cells, and none of them contain overlapping barcodes.
		for (MultiCell mc: result)
			Assert.assertSame(mc.getCellBarcodes().size(), multiplicity);


	}


}
