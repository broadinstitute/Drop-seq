/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils.readiterators;

import java.io.File;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.junit.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;

/**
 * NOTE: TagValueProcessor no longer exists, but this tests the functionality that replaces it
 */
public class TagValueProcessorTest {

	File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/testTagSorting.bam");

  @Test(enabled=true)
  public void testFilterKeepValues() {

      final Set<String> values = new HashSet<>();
      values.add("TGGCGAAGAGAT");
      final Iterator<SAMRecord> it = TagOrderIteratorTest.getTagOrderIterator(IN_FILE, "ZC");
      final Iterator<SAMRecord> toi = new FilteredIterator<SAMRecord>(it) {
          @Override
          public boolean filterOut(final SAMRecord rec) {
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

	  final Set<String> values = new HashSet<>();
	  values.add("ATCAGGGACAGA");
      final Iterator<SAMRecord> it = TagOrderIteratorTest.getTagOrderIterator(IN_FILE, "ZC");
      final Iterator<SAMRecord> toi = new FilteredIterator<SAMRecord>(it) {
          @Override
          public boolean filterOut(final SAMRecord rec) {
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
