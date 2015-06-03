package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Assert;
import org.testng.annotations.Test;

public class TagValueProcessorTest {
	
	File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/testTagSorting.bam");
	
  @Test(enabled=true)
  public void testFilterKeepValues() {
	 
	  Set<String> values = new HashSet<String>();
	  values.add("TGGCGAAGAGAT");
	  TagValueProcessor p = new TagValueProcessor("ZC", values, true);
	  
	  List<String> sortingTags = new ArrayList<String>();
	  sortingTags.add("ZC");
	  
	  TagOrderIterator toi = new TagOrderIterator(IN_FILE, sortingTags, sortingTags, p, true);
	  int counter=0;
	  
	  String [] cellOrder={"TGGCGAAGAGAT", "TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT"};
	  
	  
	  while (toi.hasNext()) {
		  SAMRecord r = toi.next();
		  String readName = r.getReadName();
		  String cellName = r.getStringAttribute("ZC");
		  
		  String expectedName = cellOrder[counter];
		  Assert.assertEquals(cellName, expectedName);
		  counter++;
	  }
	  
  }
  
  @Test(enabled=true)
  public void testFilterRejectValues() {
	 
	  Set<String> values = new HashSet<String>();
	  values.add("ATCAGGGACAGA");
	  TagValueProcessor p = new TagValueProcessor("ZC", values, false);
	  
	  List<String> sortingTags = new ArrayList<String>();
	  sortingTags.add("ZC");
	  
	  TagOrderIterator toi = new TagOrderIterator(IN_FILE, sortingTags, sortingTags, p, true);
	  int counter=0;
	  
	  String [] cellOrder={"TGGCGAAGAGAT", "TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT"};
	  
	  
	  while (toi.hasNext()) {
		  SAMRecord r = toi.next();
		  String readName = r.getReadName();
		  String cellName = r.getStringAttribute("ZC");
		  
		  String expectedName = cellOrder[counter];
		  Assert.assertEquals(cellName, expectedName);
		  counter++;
	  }
	  
  }
  
}
