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
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.MultiComparator;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

/**
 * Note that class AggregatedTagOrderIterator no longer exists, but this class tests the various classes
 * that replace the deleted class.
 *
 * @author nemesh
 *
 */
public class AggregatedTagOrderIteratorTest {
	// See TagOrderIteratorTest for more info about the data in the test BAM
	File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/testTagSorting.bam");

    private Iterator<List<SAMRecord>> filterSortAndGroupByTags(final File bamFile, final String...tags) {
        final SamReader reader = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(bamFile);
        final Iterator<SAMRecord> filteringIterator = new MissingTagFilteringIterator(reader.iterator(), tags);
        final List<Comparator<SAMRecord>> comparators = new ArrayList<>(tags.length);
        for (final String tag: tags) {
            comparators.add(new StringTagComparator(tag));
        }
        final MultiComparator<SAMRecord> comparator = new MultiComparator<>(comparators);
        final CloseableIterator<SAMRecord> sortingIterator =
                SamRecordSortingIteratorFactory.create(reader.getFileHeader(), filteringIterator, comparator, null);
        return new GroupingIterator<>(sortingIterator, comparator);
    }

    private Iterator<List<SAMRecord>> filterSortAndGroupByTagsAndQuality(final File bamFile, final String...tags) {
        final SamReader reader = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(bamFile);
        // Stack the two filters
        final Iterator<SAMRecord> filteringIterator = new MapQualityFilteredIterator(new MissingTagFilteringIterator(reader.iterator(), tags), 10, true);
        final List<Comparator<SAMRecord>> comparators = new ArrayList<>(tags.length);
        for (final String tag: tags) {
            comparators.add(new StringTagComparator(tag));
        }
        final MultiComparator<SAMRecord> comparator = new MultiComparator<>(comparators);
        final CloseableIterator<SAMRecord> sortingIterator =
                SamRecordSortingIteratorFactory.create(reader.getFileHeader(), filteringIterator, comparator, null);
        return new GroupingIterator<>(sortingIterator, comparator);
    }

    @Test(enabled=true)
	public void testGeneSorting() {
        final Iterator<List<SAMRecord>> groupingIterator = filterSortAndGroupByTags(IN_FILE, "GE");
		int counter=0;
		
		String [] geneOrder={"CHUK", "NKTR", "SNRPA1"};
		Map<String, Integer> geneCounts = new HashMap<String, Integer>();
		geneCounts.put("CHUK", 4);
		geneCounts.put("NKTR", 4);
		geneCounts.put("SNRPA1", 4);
		
		 while (groupingIterator.hasNext()) {
			  Collection<SAMRecord> r = groupingIterator.next();
			  int size = r.size();
			  SAMRecord rec = r.iterator().next();
			  
			  String gene = rec.getStringAttribute("GE");
			  System.out.println("Gene [" + gene +"] size [" + size +"]");
			  
			  String expectedGene = geneOrder[counter];
			  int expectedSize = geneCounts.get(gene);
			  
			  Assert.assertEquals(gene, expectedGene);
			  Assert.assertEquals(size, expectedSize);
			  counter++;
		 }
	}
	
	@Test(enabled=true)
	public void testCellSorting() {
        final Iterator<List<SAMRecord>> groupingIterator = filterSortAndGroupByTags(IN_FILE, "ZC");
		int counter=0;
		
		String [] cellOrder={"ATCAGGGACAGA", "TGGCGAAGAGAT"};
		Map<String, Integer> cellCounts = new HashMap<String, Integer>();
		cellCounts.put("ATCAGGGACAGA", 6);
		cellCounts.put("TGGCGAAGAGAT", 6);
		
		
		 while (groupingIterator.hasNext()) {
			  Collection<SAMRecord> r = groupingIterator.next();
			  int size = r.size();
			  SAMRecord rec = r.iterator().next();
			  
			  String cell = rec.getStringAttribute("ZC");
			  System.out.println("Cell [" + cell +"] size [" + size +"]");
			  
			  String expectedGene = cellOrder[counter];
			  int expectedSize = cellCounts.get(cell);
			  
			  Assert.assertEquals(cell, expectedGene);
			  Assert.assertEquals(size, expectedSize);
			  counter++;
		 }
	}
	
	@Test(enabled=true)
	public void testGeneCellSorting() {
        final Iterator<List<SAMRecord>> groupingIterator = filterSortAndGroupByTags(IN_FILE, "GE", "ZC");
		int counter=0;
		  
		
		String [] geneOrder={"CHUK", "CHUK", "NKTR", "NKTR", "SNRPA1", "SNRPA1"};
		String [] cellOrder={"ATCAGGGACAGA", "TGGCGAAGAGAT", "ATCAGGGACAGA", "TGGCGAAGAGAT", "ATCAGGGACAGA", "TGGCGAAGAGAT"};
		int [] expectedSize = {2,2,2,2,2,2};
		
		while (groupingIterator.hasNext()) {
			Collection<SAMRecord> r = groupingIterator.next();
			int size = r.size();
			SAMRecord rec = r.iterator().next();
			
			String readName = rec.getReadName();
			String cellName = rec.getStringAttribute("ZC");
			String geneName = rec.getStringAttribute("GE");
			
			System.out.println("Cell [" + cellName +"] geneName [" + geneName +"] size [" + size +"]");

			Assert.assertEquals(cellName, cellOrder[counter]);
			Assert.assertEquals(geneName, geneOrder[counter]);
			Assert.assertEquals(size, expectedSize[counter]);
			counter++;
		}
	}
	
	
	
	/****
	 * SOME TESTS BY #READS ON A LARGER BAM
	 */
	
	@Test (enabled=true)
	public void testGetCellGeneBatch() {
		File f = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.bam");
        String cellTag = "ZC";
        String geneExonTag = "GE";
        final Iterator<List<SAMRecord>> iter = filterSortAndGroupByTagsAndQuality(f, geneExonTag, cellTag);

		List<String>sortingTags = new ArrayList<String>();

		sortingTags.add(geneExonTag);
		sortingTags.add(cellTag);


		while (iter.hasNext()) {
			Collection<SAMRecord> recs=iter.next();
			SAMRecord r = recs.iterator().next();
			int setSize = recs.size();
			String cell = r.getStringAttribute(cellTag);
			String geneExon = r.getStringAttribute(geneExonTag); 
			int expectedSize = getCellGeneBatchExpectedSize(geneExon, cell);
			Assert.assertTrue(testAllRecordsSamTags(recs, sortingTags, expectedSize, setSize));
			Assert.assertEquals(expectedSize, setSize);
		}
		CloserUtil.close(iter);
	}
	
	@Test (enabled=true)
	public void testGetCellBatch() {
		File f = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.bam");
		String cellTag = "ZC";
        final Iterator<List<SAMRecord>> iter = filterSortAndGroupByTagsAndQuality(f, cellTag);

		List<String>sortingTags = new ArrayList<String>();
		
		sortingTags.add(cellTag);

		while (iter.hasNext()) {
			List<SAMRecord> recs=iter.next();
			SAMRecord r = recs.iterator().next();
			int setSize = recs.size();
			String cell = r.getStringAttribute(cellTag);
			int expectedSize = getCellBatchExpectedSize(cell);
			
			Assert.assertTrue(testAllRecordsSamTags(recs, sortingTags, expectedSize, setSize));
			
			Assert.assertEquals(expectedSize, setSize);
		}
        CloserUtil.close(iter);
	}
	
	
	private boolean testAllRecordsSamTags(Collection<SAMRecord> recs, List<String> sortingTags, int expectedSize, int setSize) {
		
		SAMRecord firstRec = recs.iterator().next();
		List<String> values = getValuesForTags(sortingTags, firstRec);
		System.out.println(buildLogString(values, sortingTags)+ "expected [" + expectedSize +"] observed [" + setSize + "]");
		for (SAMRecord r: recs) {
			List<String> newVals = getValuesForTags(sortingTags, firstRec);
			boolean equalFlag = testTagsEqual(values, newVals);
			if (equalFlag) return (true);
		}
		return true;		
	}
	
	private List<String> getValuesForTags(List<String>tags, SAMRecord r) {
		List<String> currentValues = new ArrayList<String>();
		for (String t: tags) {
			currentValues.add(r.getStringAttribute(t));
		}
		return (currentValues);
	}
	
	private boolean testTagsEqual (List<String> original, List<String> next) {
		for (int i=0; i<original.size(); i++) {
			String s1 = original.get(i);
			String s2 = next.get(i);
			if (!s1.equals(s2)) return (false);
		}
		return (true);
	}
	
	private String buildLogString (List<String> values, List<String> sortingTags) {
		Assert.assertEquals(values.size(), sortingTags.size());
		
		StringBuilder b = new StringBuilder();
		for (int i=0; i<values.size(); i++) {
			b.append(sortingTags.get(i));
			b.append(" [" + values.get(i) +"] ");
		}
		return (b.toString());
	}
	
	// samtools view -q 10 5cell3gene.bam |grep ZC:Z:ATCAGGGACAGA |grep GE:Z:HUMAN_3:42642106-42690227:NKTR |wc -l
	// for all 5 cells / 3 genes.
	// cells=( ZC:Z:ATCAGGGACAGA AGGGAAAATTGA TTGCCTTACGCG TGGCGAAGAGAT TACAATTAAGGC ); genes=( HUMAN_3:42642106-42690227:NKTR HUMAN_10:101948055-101989376:CHUK HUMAN_15:101821715-101835487:SNRPA1 ); for c in "${cells[@]}"; do for g in "${genes[@]}"; do r=$(samtools view -q 10 -F 260 5cell3gene.bam |grep ${c} |grep ${g} |wc -l); echo ${c} ${g} ${r}; done; done
	// assumes filtering to primary reads map quality >=10. 
	private int getCellGeneBatchExpectedSize (String gene, String cell) {
		
		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 845;
		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 200;		
		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 225;
		
		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 428;
		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 212;
		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 238;
		
		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 473;
		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 62;
		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 581;
		
		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 612;
		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 12;
		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 166;
		
		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 385;
		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 160;
		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 294;
		return -1;
		
	}
	
	//cells=( ZC:Z:ATCAGGGACAGA AGGGAAAATTGA TTGCCTTACGCG TGGCGAAGAGAT TACAATTAAGGC ); for c in "${cells[@]}"; do r=$(samtools view -q 10 -F 260 5cell3gene.bam |grep ${c} |wc -l); echo ${c} ${r}; done
	private int getCellBatchExpectedSize (String cell) {
		
		if (cell.equals("ATCAGGGACAGA")) return 1270;
		if (cell.equals("AGGGAAAATTGA")) return 878;
		if (cell.equals("TTGCCTTACGCG")) return 1116;
		if (cell.equals("TGGCGAAGAGAT")) return 790;
		if (cell.equals("TACAATTAAGGC")) return 839;		
		return -1;
		
	}
}
