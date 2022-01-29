package org.broadinstitute.dropseqrna.utils;

import java.io.File;
import java.io.IOException;

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
