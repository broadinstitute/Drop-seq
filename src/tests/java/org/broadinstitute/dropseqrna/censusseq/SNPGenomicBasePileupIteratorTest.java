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
package org.broadinstitute.dropseqrna.censusseq;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.dropseqrna.censusseq.SNPGenomicBasePileUp;
import org.broadinstitute.dropseqrna.censusseq.SNPGenomicBasePileupIterator;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

public class SNPGenomicBasePileupIteratorTest {

	/**
	 * Read pair 1 will be mismatched at the SNP with A/T
	 * Read pair 2 will be matched at the SNP with the best base quality being 30
	 * Read 3 will be unpaired T
	 * Read 4 will be unpaired A
	 * Read 5 will not overlap the SNP.
	 */
	private final File smallBAMFile = new File(
			"testdata/org/broadinstitute/dropseq/censusseq/genomic_pileup_test.sam");

	@Test(enabled=true)
	public void testAllReadsPileup() {
		SamReader reader = SamReaderFactory.makeDefault().open(smallBAMFile);
		int snpPos=23816120;
		Interval snpInterval = new Interval("1", snpPos, snpPos, true, "test");
		IntervalList intervalList = new IntervalList(reader.getFileHeader());
		intervalList.add(snpInterval);
		SNPGenomicBasePileupIterator snpIter = new SNPGenomicBasePileupIterator(reader, intervalList, "ZS", 10, null, null);

		// there's just 1 pileup.
		SNPGenomicBasePileUp p = snpIter.next();

		// test bases.
		List<Character> bases = p.getBasesAsCharacters();
		List<Character> expectedBases = new ArrayList<>();
		expectedBases.add(new Character ('T'));
		expectedBases.add(new Character ('T'));
		expectedBases.add(new Character ('A'));
		Assert.assertEquals(bases, expectedBases);

		// test qualities
		List<Byte> quals = p.getQualities();
		List<Byte> expectedQuals = new ArrayList<>();
		expectedQuals.add(new Byte ((byte)30));
		expectedQuals.add(new Byte ((byte)28));
		expectedQuals.add(new Byte ((byte)28));
		Assert.assertEquals(quals, expectedQuals);

		// test checking for next object (there isn't any) then grabbing it and receiving a null return.
		Assert.assertFalse(snpIter.hasNext());
		Assert.assertTrue(snpIter.next()==null);

		snpIter.close();
	}

	@Test (expectedExceptions=UnsupportedOperationException.class)
	public void testRemove () {
		SamReader reader = SamReaderFactory.makeDefault().open(smallBAMFile);
		int snpPos=23816120;
		Interval snpInterval = new Interval("1", snpPos, snpPos, true, "test");
		IntervalList intervalList = new IntervalList(reader.getFileHeader());
		intervalList.add(snpInterval);
		SNPGenomicBasePileupIterator snpIter = new SNPGenomicBasePileupIterator(reader, intervalList, "ZS", 10, null, null);

		snpIter.remove();
	}
}
