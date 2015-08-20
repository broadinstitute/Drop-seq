package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.junit.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * NOTE: TagValueProcessor no longer exists, but this tests the functionality that replaces it
 */
public class TagValueProcessorTest {
	
	File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/testTagSorting.bam");
	
  @Test(enabled=true)
  public void testFilterKeepValues() {

      final Set<String> values = new HashSet<String>();
      values.add("TGGCGAAGAGAT");
      final Iterator<SAMRecord> it = TagOrderIteratorTest.getTagOrderIterator(IN_FILE, "ZC");
      final Iterator<SAMRecord> toi = new FilteredIterator<SAMRecord>(it) {
          @Override
          protected boolean filterOut(SAMRecord rec) {
              return !values.contains(rec.getAttribute("ZC"));
          }
      };
	  int counter=0;
	  
	  String [] cellOrder={"TGGCGAAGAGAT", "TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT"};
	  
	  
	  while (toi.hasNext()) {
		  SAMRecord r = toi.next();
		  r.getReadName();
		  String cellName = r.getStringAttribute("ZC");
		  
		  String expectedName = cellOrder[counter];
		  Assert.assertEquals(cellName, expectedName);
		  counter++;
	  }
	  
  }
  
  @Test(enabled=true)
  public void testFilterRejectValues() {
	 
	  final Set<String> values = new HashSet<String>();
	  values.add("ATCAGGGACAGA");
      final Iterator<SAMRecord> it = TagOrderIteratorTest.getTagOrderIterator(IN_FILE, "ZC");
      final Iterator<SAMRecord> toi = new FilteredIterator<SAMRecord>(it) {
          @Override
          protected boolean filterOut(SAMRecord rec) {
              return values.contains(rec.getAttribute("ZC"));
          }
      };
	  int counter=0;
	  
	  String [] cellOrder={"TGGCGAAGAGAT", "TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT"};
	  
	  
	  while (toi.hasNext()) {
		  SAMRecord r = toi.next();
		  r.getReadName();
		  String cellName = r.getStringAttribute("ZC");
		  
		  String expectedName = cellOrder[counter];
		  Assert.assertEquals(cellName, expectedName);
		  counter++;
	  }
	  
  }
  
}
