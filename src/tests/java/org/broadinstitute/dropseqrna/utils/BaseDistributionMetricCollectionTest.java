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
package org.broadinstitute.dropseqrna.utils;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.ValidationStringency;
import org.testng.Assert;
import org.testng.annotations.Test;


public class BaseDistributionMetricCollectionTest {

  @Test
  public void readBaseDistributionTest() throws IOException {
	  
	  // round trip BaseDistributionMetricCollection
	  File temp =  File.createTempFile("BaseDistributionMetricCollectionTest", ".txt");
	  temp.deleteOnExit();
	  
	  BaseDistributionMetricCollection c = new BaseDistributionMetricCollection();
	  c.addBase ('A', 1, 1);
	  c.addBase ('C', 2, 4);
	  c.addBase ('G', 3, 9);
	  c.addBase ('T', 4, 16);	  	  
	  c.addBases(new byte [] {'A', 'C'});
	  c.addBases(new char [] {'G', 'T'});
	  c.addBases("ACGT");
	  
	  
	  c.writeOutput(temp);
	  
	  BaseDistributionMetricCollection c2 = BaseDistributionMetricCollection.readBaseDistribution(temp);
	  Assert.assertTrue(c.distributionsEqual(c2));	  	  	  
  }

  @Test(expectedExceptions = IllegalArgumentException.class)
  public void readBaseDistributionTestStrict() throws IOException {
	  File temp =  File.createTempFile("BaseDistributionMetricCollectionTestStrict", ".txt");
	  temp.deleteOnExit();

	  BaseDistributionMetricCollection c = new BaseDistributionMetricCollection(ValidationStringency.STRICT);
	  c.addBase ('A', 1, 1);
	  c.addBase ('C', 2, 4);
	  c.addBase ('G', 3, 9);
	  c.addBase ('T', 4, 16);
	  c.addBase ('-', 5, 25);
  }

  @Test
  public void readBaseDistributionTestLenient() throws IOException {
	  // round trip BaseDistributionMetricCollection
	  File temp =  File.createTempFile("BaseDistributionMetricCollectionTestLenient", ".txt");
	  temp.deleteOnExit();

	  BaseDistributionMetricCollection c = new BaseDistributionMetricCollection(ValidationStringency.LENIENT);
	  c.addBase ('A', 1, 1);
	  c.addBase ('C', 2, 4);
	  c.addBase ('G', 3, 9);
	  c.addBase ('T', 4, 16);
	  c.addBase ('-', 5, 25);
	  c.addBases(new byte [] {'A', 'C'});
	  c.addBases(new char [] {'G', 'T', '-'});
	  c.addBases("ACGT-");

	  c.writeOutput(temp);

	  BaseDistributionMetricCollection c2 = BaseDistributionMetricCollection.readBaseDistribution(temp);
	  Assert.assertTrue(c.distributionsEqual(c2));
  }

  @Test
  public void readBaseDistributionTestSilent() throws IOException {
	  // round trip BaseDistributionMetricCollection
	  File temp =  File.createTempFile("BaseDistributionMetricCollectionTestSilent", ".txt");
	  temp.deleteOnExit();

	  BaseDistributionMetricCollection c = new BaseDistributionMetricCollection(ValidationStringency.SILENT);
	  c.addBase ('A', 1, 1);
	  c.addBase ('C', 2, 4);
	  c.addBase ('G', 3, 9);
	  c.addBase ('T', 4, 16);
	  c.addBase ('-', 5, 25);
	  c.addBases(new byte [] {'A', 'C'});
	  c.addBases(new char [] {'G', 'T', '-'});
	  c.addBases("ACGT-");

	  c.writeOutput(temp);

	  BaseDistributionMetricCollection c2 = BaseDistributionMetricCollection.readBaseDistribution(temp);
	  Assert.assertTrue(c.distributionsEqual(c2));
  }

  @Test
  public void testMergeMetricCollections () {
	  BaseDistributionMetricCollection c = new BaseDistributionMetricCollection();
	  c.addBase ('A', 1, 1);
	  c.addBase('T', 1, 1);
	  
	  BaseDistributionMetricCollection c2 = new BaseDistributionMetricCollection();
	  c2.addBase ('C', 2, 4);
	  c2.addBase('T', 1, 1);
	  
	  BaseDistributionMetricCollection all = new BaseDistributionMetricCollection();
	  all.addBase ('A', 1, 1);
	  all.addBase ('C', 2, 4);
	  all.addBase('T', 1, 2);
	  
	  c.mergeMetricCollections(c2);
	  Assert.assertTrue(c.distributionsEqual(all));
	  
	  
	  
  }
  
}
