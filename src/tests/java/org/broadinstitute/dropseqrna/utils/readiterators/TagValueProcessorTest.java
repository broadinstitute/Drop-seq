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
