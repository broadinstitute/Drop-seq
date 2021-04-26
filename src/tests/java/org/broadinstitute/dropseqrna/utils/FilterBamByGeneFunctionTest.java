package org.broadinstitute.dropseqrna.utils;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;

import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import junit.framework.Assert;
import picard.annotation.LocusFunction;

public class FilterBamByGeneFunctionTest {

	/**
	 * Data generated from:
	 * samtools view -H testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_5_cell_2_snp_testdata_retagged.bam > testdata/org/broadinstitute/dropseq/utils/FilterBamByGeneFunction.sam
     * samtools view testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_5_cell_2_snp_testdata_retagged.bam |grep INTRONIC >> testdata/org/broadinstitute/dropseq/utils/FilterBamByGeneFunction.sam
	 */
	private File INPUT= new File ("testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_5_cell_2_snp_testdata_retagged.bam");
	
	private File SMALL_UNFILTERED_RESULT=new File ("testdata/org/broadinstitute/dropseq/utils/SmallUnfilteredResult.bam");
	private File SMALL_UNFILTERED_EMPTY_RESULT=new File ("testdata/org/broadinstitute/dropseq/utils/SmallUnfilteredEmptyResult.bam");
	
	// 2 reads that are both INTRONIC and CODING.
	private File SMALL=new File ("testdata/org/broadinstitute/dropseq/utils/FilterBamByGeneFunction.sam");
	
	// Small tests of 2 reads, make sure reads tagging two functions are retained when filtering to either.
	
	@Test
	public void filterToCodingSmall() {
		// retain 2 reads
		FilterBamByGeneFunction f = new FilterBamByGeneFunction();
		f.INPUT=Collections.singletonList(SMALL);
		f.OUTPUT= TestUtils.getTempReportFile("small", ".bam");
		f.LOCUS_FUNCTION_LIST=Arrays.asList(LocusFunction.CODING);
		int result = f.doWork();
		Assert.assertEquals(0, result);
		TestUtils.assertSamFilesSame(f.OUTPUT, SMALL_UNFILTERED_RESULT);
		
	}
	
	@Test (enabled=true)
	public void filterToUTRSmall() {
		FilterBamByGeneFunction f = new FilterBamByGeneFunction();
		f.INPUT=Collections.singletonList(SMALL);
		f.OUTPUT= TestUtils.getTempReportFile("small", ".bam");
		f.LOCUS_FUNCTION_LIST=Arrays.asList(LocusFunction.UTR);
		int result = f.doWork();
		Assert.assertEquals(0, result);
		TestUtils.assertSamFilesSame(f.OUTPUT, SMALL_UNFILTERED_EMPTY_RESULT);
				
	}
	
	@Test (enabled=true)
	public void filterToIntronicSmall () {
		FilterBamByGeneFunction f = new FilterBamByGeneFunction();
		f.INPUT=Collections.singletonList(SMALL);
		f.OUTPUT= TestUtils.getTempReportFile("small", ".bam");
		f.LOCUS_FUNCTION_LIST=Arrays.asList(LocusFunction.INTRONIC);
		int result = f.doWork();
		Assert.assertEquals(0, result);
		TestUtils.assertSamFilesSame(f.OUTPUT, SMALL_UNFILTERED_RESULT);				
	}
	
	@Test (enabled=true)
	public void filterToUTR() {
		FilterBamByGeneFunction f = new FilterBamByGeneFunction();
		f.INPUT=Collections.singletonList(INPUT);
		f.OUTPUT= TestUtils.getTempReportFile("small", ".bam");
		f.LOCUS_FUNCTION_LIST=Arrays.asList(LocusFunction.UTR);
		int result = f.doWork();
		Assert.assertEquals(0, result);
		int numReads = getNumReads(f.OUTPUT);
		Assert.assertEquals(4, numReads);		
	}
	
	@Test (enabled=true)
	public void filterToIntronic() {
		FilterBamByGeneFunction f = new FilterBamByGeneFunction();
		f.INPUT=Collections.singletonList(INPUT);
		f.OUTPUT= TestUtils.getTempReportFile("small", ".bam");
		f.LOCUS_FUNCTION_LIST=Arrays.asList(LocusFunction.INTRONIC);
		int result = f.doWork();
		Assert.assertEquals(0, result);
		int numReads = getNumReads(f.OUTPUT);
		Assert.assertEquals(2, numReads);		
	}
	
	@Test (enabled=true)
	public void filterToCoding() {
		FilterBamByGeneFunction f = new FilterBamByGeneFunction();
		f.INPUT=Collections.singletonList(INPUT);
		f.OUTPUT= TestUtils.getTempReportFile("small", ".bam");
		f.LOCUS_FUNCTION_LIST=Arrays.asList(LocusFunction.CODING);
		int result = f.doWork();
		Assert.assertEquals(0, result);
		int numReads = getNumReads(f.OUTPUT);
		
		// 18 reads have the CODING annotation but are  on the wrong strand!
		// 765 (total reads) -18 = 747
		
		Assert.assertEquals(747, numReads);		
	}
	
	
	
	
	
	private int getNumReads (File f) {
		final SamReader actualReader = SamReaderFactory.makeDefault().open(f);
		final SAMRecordIterator actualIterator = actualReader.iterator();						
		int count = (int) actualIterator.stream().count();				
		return count;		
	}
	
	
}
