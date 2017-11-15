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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.dropseqrna.utils.MultiComparator;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
import org.junit.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

/**
 * NOTE: TagOrderIterator no longer exists, but this tests the functionality of the classes that replace it.
 */
public class TagOrderIteratorTest {
	
	//File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.bam");
	File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/testTagSorting.bam");
	
	/**
	 * Make a BAM file that has a set of reads representing multiple cells/genes/molecular barcodes to test different iteration patterns.
	 * Extract the following reads from 5cell3gene.bam (below is the data in cell/gene/molBC order.)

	 * NS500217:67:H14GMBGXX:3:22612:4003:5921   ZC:Z:ATCAGGGACAGA GE:Z:CHUK    XM:Z:CGAATTTT
	 * NS500217:67:H14GMBGXX:2:22312:14606:19968 ZC:Z:ATCAGGGACAGA GE:Z:CHUK    XM:Z:CGATTTTT
	 * NS500217:67:H14GMBGXX:1:12306:21948:16777 ZC:Z:ATCAGGGACAGA GE:Z:NKTR	XM:Z:GACGAGGG
	 * NS500217:67:H14GMBGXX:1:22312:20727:11418 ZC:Z:ATCAGGGACAGA GE:Z:NKTR	XM:Z:GCCGGAGC
	 * NS500217:67:H14GMBGXX:3:13604:8267:15833  ZC:Z:ATCAGGGACAGA GE:Z:SNRPA1  XM:Z:CGGCCTGG
	 * NS500217:67:H14GMBGXX:3:11507:4353:14582  ZC:Z:ATCAGGGACAGA GE:Z:SNRPA1  XM:Z:TTGGGGCA
	 * 
	 * NS500217:67:H14GMBGXX:4:21505:26726:10699 ZC:Z:TGGCGAAGAGAT GE:Z:CHUK    XM:Z:GGGTGAGA
	 * NS500217:67:H14GMBGXX:1:13105:19089:17127 ZC:Z:TGGCGAAGAGAT GE:Z:CHUK    XM:Z:TGGAGTTT
	 * NS500217:67:H14GMBGXX:1:22207:3769:12483  ZC:Z:TGGCGAAGAGAT GE:Z:NKTR    XM:Z:GCAGGGCG
	 * NS500217:67:H14GMBGXX:3:21507:14155:11331 ZC:Z:TGGCGAAGAGAT GE:Z:NKTR    XM:Z:GGCGGGTT
	 * NS500217:67:H14GMBGXX:1:22306:14005:17334 ZC:Z:TGGCGAAGAGAT GE:Z:SNRPA1  XM:Z:AAGGTTGG
	 * NS500217:67:H14GMBGXX:2:23305:7629:7836   ZC:Z:TGGCGAAGAGAT GE:Z:SNRPA1  XM:Z:TTCTTCTA
	 * 
	 * 
	 */
	
	
	
  @Test(enabled=true)
  public void testGeneSorting() {
      final Iterator<SAMRecord> toi = getTagOrderIterator(IN_FILE, "GE");
	  int counter=0;
	  
	  String [] geneOrder={"CHUK", "CHUK", "CHUK", "CHUK", "NKTR", "NKTR", "NKTR", "NKTR", "SNRPA1", "SNRPA1", "SNRPA1", "SNRPA1"};
	  
	  while (toi.hasNext()) {
		  SAMRecord r = toi.next();
		  r.getReadName();
		  String geneName = r.getStringAttribute("GE");
		  String expectedName = geneOrder[counter];
		  Assert.assertEquals(geneName, expectedName);
		  counter++;
	  }
  }

    public static Iterator<SAMRecord> getTagOrderIterator(final File inFile, final String...tags) {
        final SamReader reader = SamReaderFactory.makeDefault().open(inFile);
        final MissingTagFilteringIterator filteringIterator = new MissingTagFilteringIterator(reader.iterator(), tags);
        final List<Comparator<SAMRecord>> comparators = new ArrayList<>(tags.length);
        for (final String tag : tags) {
            comparators.add(new StringTagComparator(tag));
        }
        final Comparator<SAMRecord> comparator = new MultiComparator<>(comparators);
        return SamRecordSortingIteratorFactory.create(reader.getFileHeader(), filteringIterator, comparator, null);
    }

    @Test(enabled=true)
  public void testCellSorting() {
        final Iterator<SAMRecord> toi = getTagOrderIterator(IN_FILE, "ZC");
	  int counter=0;
	  
	  String [] cellOrder={"ATCAGGGACAGA","ATCAGGGACAGA","ATCAGGGACAGA","ATCAGGGACAGA","ATCAGGGACAGA","ATCAGGGACAGA","TGGCGAAGAGAT", "TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT"};
	  
	  
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
  public void testGeneCellSorting() {
      final Iterator<SAMRecord> toi = getTagOrderIterator(IN_FILE, "GE", "ZC");
	  int counter=0;
	  
	  String [] cellOrder={"ATCAGGGACAGA", "ATCAGGGACAGA", "TGGCGAAGAGAT", "TGGCGAAGAGAT","ATCAGGGACAGA", "ATCAGGGACAGA", "TGGCGAAGAGAT", "TGGCGAAGAGAT","ATCAGGGACAGA", "ATCAGGGACAGA", "TGGCGAAGAGAT", "TGGCGAAGAGAT"};
	  String [] geneOrder={"CHUK", "CHUK", "CHUK", "CHUK", "NKTR", "NKTR", "NKTR", "NKTR", "SNRPA1", "SNRPA1", "SNRPA1", "SNRPA1"};
	  while (toi.hasNext()) {
		  SAMRecord r = toi.next();
		  r.getReadName();
		  String cellName = r.getStringAttribute("ZC");
		  String geneName = r.getStringAttribute("GE");
		  System.out.println("Gene [" + geneName +"] Cell [" + cellName +"]");
		  
		  Assert.assertEquals(cellName, cellOrder[counter]);
		  Assert.assertEquals(geneName, geneOrder[counter]);
		  counter++;
	  }
  }
  
  @Test(enabled=true)
  public void testCellGeneMolecularSorting() {
      final Iterator<SAMRecord> toi = getTagOrderIterator(IN_FILE, "ZC", "GE", "XM");
	  int counter=0;
	  
	  String [] cellOrder={"ATCAGGGACAGA", "ATCAGGGACAGA", "ATCAGGGACAGA", "ATCAGGGACAGA", "ATCAGGGACAGA", "ATCAGGGACAGA", "TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT"};
	  String [] geneOrder={"CHUK", "CHUK", "NKTR", "NKTR", "SNRPA1", "SNRPA1","CHUK", "CHUK", "NKTR", "NKTR", "SNRPA1", "SNRPA1"};
	  String [] molBCOrder={"CGAATTTT", "CGATTTTT", "GACGAGGG", "GCCGGAGC", "CGGCCTGG", "TTGGGGCA", "GGGTGAGA", "TGGAGTTT", "GCAGGGCG", "GGCGGGTT", "AAGGTTGG", "TTCTTCTA"};
	  
		while (toi.hasNext()) {
			SAMRecord r = toi.next();
			r.getReadName();
			String cellName = r.getStringAttribute("ZC");
			String geneName = r.getStringAttribute("GE");
			String molBC = r.getStringAttribute("XM");

			System.out.println("Cell [" + cellName + "] Gene [" + geneName
					+ "] MolBC [" + molBC + "]");

			Assert.assertEquals(cellName, cellOrder[counter]);
			Assert.assertEquals(geneName, geneOrder[counter]);
			Assert.assertEquals(molBC, molBCOrder[counter]);
			counter++;
		}
  }
  
  
  
  
  
}
