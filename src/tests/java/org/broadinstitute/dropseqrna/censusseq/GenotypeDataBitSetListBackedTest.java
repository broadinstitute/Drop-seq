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
package org.broadinstitute.dropseqrna.censusseq;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.GenotypeType;
import org.testng.Assert;

import org.testng.annotations.Test;

import java.util.Random;

public class GenotypeDataBitSetListBackedTest {

	private static final Log log = Log.getInstance(GenotypeDataBitSetListBackedTest.class);

	@Test
	public void testAddRetrieve () {
		GenotypeDataBitSetListBacked d = new GenotypeDataBitSetListBacked(3);
		int [] values = {0,1,2,2,2,1,0,1,0};
		int counter=0;
		// insert data.
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++) {
				d.add(i, j, values[counter]);
				counter++;
			}
		// retrieve data.
		counter=0;
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++) {
				int expected = values[counter];
				int actual = d.get(i, j);
				Assert.assertEquals(expected, actual);
				counter++;
			}
		// the data ordered by variant:
		// 0,2,0 -> (0,0,0,1,0,0) -> {3}
		// 1,2,1 -> (1,0,0,1,1,0) -> {0,3,4}
		// 2,1,0 -> (0,1,1,0,0,0) -> {1,2}
		// indeed, this is true if you get a string representation of the list of bitsets:
		// [{3}, {0, 3, 4}, {1, 2}]

	}
	
	@Test
	public void testAddRetrieveGenotypeType () {
		GenotypeDataBitSetListBacked d = new GenotypeDataBitSetListBacked(3);
		int [] values = {0,1,2,2,2,1,0,1,0,-1,-1,-1};
		GenotypeType [] valuesGT = {GenotypeType.HOM_REF, GenotypeType.HET, GenotypeType.HOM_VAR, GenotypeType.HOM_VAR, GenotypeType.HOM_VAR, GenotypeType.HET, 
				GenotypeType.HOM_REF, GenotypeType.HET, GenotypeType.HOM_REF, GenotypeType.NO_CALL, GenotypeType.NO_CALL,GenotypeType.NO_CALL};
		int counter=0;
		// insert data.
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++) {
				d.add(i, j, valuesGT[counter]);
				counter++;
			}
		// retrieve data.
		counter=0;
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++) {
				int expected = values[counter];
				int actual = d.get(i, j);
				Assert.assertEquals(expected, actual);
				GenotypeType gtExpected = valuesGT[counter];
				GenotypeType actualGT = d.getAsGenotypeType(i, j);
				Assert.assertEquals(gtExpected, actualGT);
				counter++;
			}
		// the data ordered by variant:
		// 0,2,0 -> (0,0,0,1,0,0) -> {3}
		// 1,2,1 -> (1,0,0,1,1,0) -> {0,3,4}
		// 2,1,0 -> (0,1,1,0,0,0) -> {1,2}
		// indeed, this is true if you get a string representation of the list of bitsets:
		// [{3}, {0, 3, 4}, {1, 2}]

	}
	

	@Test
	public void testAddRetrieveMissingData () {
		GenotypeDataBitSetListBacked d = new GenotypeDataBitSetListBacked(3);
		int [] values = {-1,1,2,-1,2,1,-1,1,0};
		int counter=0;
		// insert data.
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++) {
				d.add(i, j, values[counter]);
				counter++;
			}
		// retrieve data.
		counter=0;
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++) {
				int expected = values[counter];
				int actual = d.get(i, j);
				Assert.assertEquals(expected, actual);
				counter++;
			}
		// the data ordered by variant:
		// 0,2,0 -> (0,0,0,1,0,0) -> {3}
		// 1,2,1 -> (1,0,0,1,1,0) -> {0,3,4}
		// 2,1,0 -> (0,1,1,0,0,0) -> {1,2}
		// indeed, this is true if you get a string representation of the list of bitsets:
		// [{0, 1, 2, 3, 4, 5}, {0, 3, 4}, {1, 2}]

	}


	private GenotypeDataBitSetListBacked insertData (final int [] values, final int numSamples) {
		GenotypeDataBitSetListBacked d = new GenotypeDataBitSetListBacked(numSamples);
		int counter=0;
		int numVariants=(int) Math.ceil((double) (values.length)/(double) numSamples);
		// insert data.
		for (int i=0; i<numVariants; i++)
			for (int j=0; j<numSamples; j++) {
				d.add(i, j, values[counter]);
				counter++;
			}
		return d;
	}

	private GenotypeDataBitSetListBacked generateFromRandomData (final int numSamples, final int numVariants) {
		GenotypeDataBitSetListBacked d = new GenotypeDataBitSetListBacked(numSamples);
		Random random = new Random();
		// insert data.
		for (int i=0; i<numVariants; i++) {
			if (i%1000000==0) log.info("Processed [" + i + "] variants of [" + numVariants +"]");
			for (int j=0; j<numSamples; j++)
				d.add(i, j, random.nextInt(3));
		}
		return d;
	}

	// Big memory test.
	@Test(enabled=false)
	public void testStorageSize () {
		// generate random ints.
		// test if this fits in 2G of heap memory with -Xmx2g.
		int numSamples=200;
		int numVariants = 30000000; //20M
		long dataSize=(long) numSamples* (long) numVariants*2L;
		double dataSizeGB=dataSize/8e9;
		log.info("Expected memory usage [" + dataSizeGB +"] GB");
		@SuppressWarnings("unused")
		GenotypeDataBitSetListBacked d = generateFromRandomData(numSamples, numVariants);
		log.info("Finished");
	}
}
