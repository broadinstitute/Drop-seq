package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.testng.annotations.Test;

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
	  ReadProcessorCollection filters = new ReadProcessorCollection();
	  List<String> sortingTags = new ArrayList<String>();
	  //sortingTags.add("ZC");
	  sortingTags.add("GE");
	  //sortingTags.add("XM");
	  
	  TagOrderIterator toi = new TagOrderIterator(IN_FILE, sortingTags, filters, true);
	  int counter=0;
	  
	  String [] geneOrder={"CHUK", "CHUK", "CHUK", "CHUK", "NKTR", "NKTR", "NKTR", "NKTR", "SNRPA1", "SNRPA1", "SNRPA1", "SNRPA1"};
	  
	  while (toi.hasNext()) {
		  SAMRecord r = toi.next();
		  String readName = r.getReadName();
		  String geneName = r.getStringAttribute("GE");
		  String expectedName = geneOrder[counter];
		  Assert.assertEquals(geneName, expectedName);
		  counter++;
	  }
  }
  
  @Test(enabled=true)
  public void testCellSorting() {
	  ReadProcessorCollection filters = new ReadProcessorCollection();
	  List<String> sortingTags = new ArrayList<String>();
	  
	  sortingTags.add("ZC");
	  
	  TagOrderIterator toi = new TagOrderIterator(IN_FILE, sortingTags, filters, true);
	  int counter=0;
	  
	  String [] cellOrder={"ATCAGGGACAGA","ATCAGGGACAGA","ATCAGGGACAGA","ATCAGGGACAGA","ATCAGGGACAGA","ATCAGGGACAGA","TGGCGAAGAGAT", "TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT"};
	  
	  
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
  public void testGeneCellSorting() {
	  ReadProcessorCollection filters = new ReadProcessorCollection();
	  List<String> sortingTags = new ArrayList<String>();
	  sortingTags.add("GE");
	  sortingTags.add("ZC");
	  
	  TagOrderIterator toi = new TagOrderIterator(IN_FILE, sortingTags, filters, true);
	  int counter=0;
	  
	  String [] cellOrder={"ATCAGGGACAGA", "ATCAGGGACAGA", "TGGCGAAGAGAT", "TGGCGAAGAGAT","ATCAGGGACAGA", "ATCAGGGACAGA", "TGGCGAAGAGAT", "TGGCGAAGAGAT","ATCAGGGACAGA", "ATCAGGGACAGA", "TGGCGAAGAGAT", "TGGCGAAGAGAT"};
	  String [] geneOrder={"CHUK", "CHUK", "CHUK", "CHUK", "NKTR", "NKTR", "NKTR", "NKTR", "SNRPA1", "SNRPA1", "SNRPA1", "SNRPA1"};
	  while (toi.hasNext()) {
		  SAMRecord r = toi.next();
		  String readName = r.getReadName();
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
	  ReadProcessorCollection filters = new ReadProcessorCollection();
	  List<String> sortingTags = new ArrayList<String>();	  
	  sortingTags.add("ZC");
	  sortingTags.add("GE");
	  sortingTags.add("XM");
	  
	  TagOrderIterator toi = new TagOrderIterator(IN_FILE, sortingTags, filters, true);
	  int counter=0;
	  
	  String [] cellOrder={"ATCAGGGACAGA", "ATCAGGGACAGA", "ATCAGGGACAGA", "ATCAGGGACAGA", "ATCAGGGACAGA", "ATCAGGGACAGA", "TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT","TGGCGAAGAGAT"};
	  String [] geneOrder={"CHUK", "CHUK", "NKTR", "NKTR", "SNRPA1", "SNRPA1","CHUK", "CHUK", "NKTR", "NKTR", "SNRPA1", "SNRPA1"};
	  String [] molBCOrder={"CGAATTTT", "CGATTTTT", "GACGAGGG", "GCCGGAGC", "CGGCCTGG", "TTGGGGCA", "GGGTGAGA", "TGGAGTTT", "GCAGGGCG", "GGCGGGTT", "AAGGTTGG", "TTCTTCTA"};
	  
		while (toi.hasNext()) {
			SAMRecord r = toi.next();
			String readName = r.getReadName();
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
