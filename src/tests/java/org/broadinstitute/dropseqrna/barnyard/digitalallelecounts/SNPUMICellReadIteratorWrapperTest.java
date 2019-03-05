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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMICellReadIteratorWrapper;
import org.broadinstitute.dropseqrna.utils.MultiComparator;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
import org.broadinstitute.dropseqrna.utils.readiterators.*;
import org.junit.Assert;
import org.testng.annotations.Test;
import picard.annotation.LocusFunction;

import java.io.File;
import java.util.*;

public class SNPUMICellReadIteratorWrapperTest {

	private final File bamFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_5_cell_2_snp_testdata_retagged.bam");
	private final File cellBCFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_cells_cell_barcodes.txt");
	private final File snpIntervals = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_cells_2snps.intervals");
	private final String cellBarcodeTag = "XC";
	private final String molBCTag = "XM";
	private final String snpTag = "YS";
	private final int readMQ = 10;

	private String GENE_NAME_TAG="gn";
	private String GENE_STRAND_TAG="gs";
	private String GENE_FUNCTION_TAG="gf";
	private StrandStrategy STRAND_STRATEGY = StrandStrategy.SENSE;
	private List<LocusFunction> LOCUS_FUNCTION_LIST=new ArrayList<LocusFunction>(Arrays.asList(LocusFunction.CODING, LocusFunction.UTR));

	private final File smallBAMFile=new File ("testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/smallTest_retagged.sam");

	/**
	 * Integration style test to see if I can throw null pointers, etc with a bit of real data.
	 */
	@Test(enabled=true)
	public void processReads() {
		List<String> cellBarcodeList = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
		IntervalList loci = IntervalList.fromFile(snpIntervals);

        SamReader reader = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(bamFile);
        // Filter records before sorting, to reduce I/O
        final MissingTagFilteringIterator filteringIterator =
                new MissingTagFilteringIterator(reader.iterator(), GENE_NAME_TAG, cellBarcodeTag, molBCTag);

        MapQualityFilteredIterator filteringIterator2 = new MapQualityFilteredIterator(filteringIterator, readMQ, true);

        GeneFunctionIteratorWrapper gfteratorWrapper = new GeneFunctionIteratorWrapper(filteringIterator2, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG, false, STRAND_STRATEGY, LOCUS_FUNCTION_LIST);

        SNPUMICellReadIteratorWrapper snpumiCellReadIterator = new SNPUMICellReadIteratorWrapper(gfteratorWrapper, loci, cellBarcodeTag, cellBarcodeList, GENE_NAME_TAG, snpTag, readMQ);

        // create comparators in the order the data should be sorted
        final MultiComparator<SAMRecord> multiComparator = new MultiComparator<>(
                new StringTagComparator(cellBarcodeTag),
                new StringTagComparator(GENE_NAME_TAG),
                new StringTagComparator(molBCTag));

        final CloseableIterator<SAMRecord> sortingIterator =
                SamRecordSortingIteratorFactory.create(reader.getFileHeader(), snpumiCellReadIterator, multiComparator, null);

        Assert.assertTrue(sortingIterator.hasNext());
        SAMRecord nextRead = sortingIterator.next();
        List<SAMTagAndValue> tagValues= nextRead.getAttributes();
        String snpTagValue = nextRead.getStringAttribute(this.snpTag);
		CloserUtil.close(snpumiCellReadIterator);
	}

	@Test(enabled=true)
	public void testGeneSNPSplitting () {
		List<String> cellBarcodeList = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
		IntervalList loci = IntervalList.fromFile(snpIntervals);
		// make a new SNP that is 10 bases before the first SNP.  Reads in the test data will hit both.
		Interval i1 = loci.getIntervals().iterator().next();
		Interval i2 = new Interval (i1.getContig(), i1.getStart()-10, i1.getStart()-10, true, "testSNP1");
		loci.add(i2);



		//A read that hits only the default SNP: start read at 76227020
		SAMRecord r1 = getReadByName("NS500217:67:H14GMBGXX:1:11308:22039:11268", smallBAMFile);
		Collection<SAMRecord> tempList = processOneRead(r1, loci, cellBarcodeList);
		Assert.assertTrue(tempList.size()==1);
		Assert.assertEquals(1, getNumGenes(tempList));
		Assert.assertEquals(1, getNumSNPs(tempList));

		//A read that hits both the default SNP and the new SNP: start read at 76227000
		r1 = getReadByName("NS500217:67:H14GMBGXX:1:13202:10555:15929", smallBAMFile);
		tempList = processOneRead(r1, loci, cellBarcodeList);
		Assert.assertTrue(tempList.size()==2);
		Assert.assertEquals(1, getNumGenes(tempList));
		Assert.assertEquals(2, getNumSNPs(tempList));


		//  A read that hits only the default SNP: start read at 76227020, has 2 genes.
		r1 = getReadByName("NS500217:67:H14GMBGXX:1:21206:20467:19854", smallBAMFile);
		tempList = processOneRead(r1, loci, cellBarcodeList);
		Assert.assertTrue(tempList.size()==2);
		Assert.assertEquals(2, getNumGenes(tempList));
		Assert.assertEquals(1, getNumSNPs(tempList));


		//A read that hits both the default SNP and the new SNP: start read at 76227000, has 2 genes.
		r1 = getReadByName("NS500217:67:H14GMBGXX:1:22302:3826:3320", smallBAMFile);
		tempList = processOneRead(r1, loci, cellBarcodeList);
		Assert.assertTrue(tempList.size()==4);
		Assert.assertEquals(2, getNumGenes(tempList));
		Assert.assertEquals(2, getNumSNPs(tempList));


	}

	private SAMRecord getReadByName (final String name, final File bamFile) {
		SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
		Iterator<SAMRecord> iter = reader.iterator();
		while (iter.hasNext()) {
			SAMRecord r = iter.next();
			if (r.getReadName().equals(name)) return (r);
		}
		return null;
	}

    private List<SAMRecord> processOneRead(final SAMRecord rec, final IntervalList loci, final List<String> cellBarcodeList) {

        SNPUMICellReadIteratorWrapper it = new SNPUMICellReadIteratorWrapper(Collections.singletonList(rec).iterator(), loci, cellBarcodeTag, cellBarcodeList, GENE_NAME_TAG, snpTag, readMQ);
        final ArrayList<SAMRecord> ret = new ArrayList<>();
        while (it.hasNext())
			ret.add(it.next());
        CloserUtil.close(it);
        return ret;
    }

	private int getNumGenes (final Collection<SAMRecord> recs) {
		Set<String> genes = new HashSet<String>();
		for (SAMRecord r: recs)
			genes.add(r.getStringAttribute(this.GENE_NAME_TAG));
		return genes.size();
	}

	private int getNumSNPs (final Collection<SAMRecord> recs) {
		Set<String> genes = new HashSet<String>();
		for (SAMRecord r: recs)
			genes.add(r.getStringAttribute(this.snpTag));
		return genes.size();
	}

	private void validateGeneSNPSplittingResults (final SAMRecord original, final Collection<SAMRecord> result) {
		//A read that hits only the default SNP: start read at 76227020
		if (original.getReadName().equals("NS500217:67:H14GMBGXX:1:11308:22039:11268"))
			Assert.assertTrue(result.size()==1);
		//A read that hits both the default SNP and the new SNP: start read at 76227000
		if (original.getReadName().equals("NS500217:67:H14GMBGXX:1:13202:10555:15929"))
			Assert.assertTrue(result.size()==2);
		//  A read that hits only the default SNP: start read at 76227020, has 2 genes.
		if (original.getReadName().equals("NS500217:67:H14GMBGXX:1:21206:20467:19854"))
			Assert.assertTrue(result.size()==2);
		//A read that hits both the default SNP and the new SNP: start read at 76227000, has 2 genes.
		if (original.getReadName().equals("NS500217:67:H14GMBGXX:1:22302:3826:3320"))
			Assert.assertTrue(result.size()==4);

	}


}
