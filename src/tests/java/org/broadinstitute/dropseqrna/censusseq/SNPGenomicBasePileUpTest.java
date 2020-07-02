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
import java.util.Iterator;
import java.util.List;

import org.broadinstitute.dropseqrna.censusseq.SNPGenomicBasePileUp;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;
import junit.framework.Assert;

public class SNPGenomicBasePileUpTest {


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
	public void testReadsMismatchedBase() {
		// read pair 1 test.  Expected no reads in final pileup.
		int snpPos=23816120;
		Interval snpInterval = new Interval("1", snpPos, snpPos, true, "test");
		SNPGenomicBasePileUp p = new SNPGenomicBasePileUp(snpInterval);

		List<SAMRecord> recs = getReadsByName("H3FF2CCXX150427:1:2107:16631:1", this.smallBAMFile);
		p.addReads(recs);

		p.finalizePileup();
		List<Character> bases = p.getBasesAsCharacters();
		List<Character> expectedBases = new ArrayList<>();
		Assert.assertEquals(expectedBases, bases);
	}

	@Test
	public void testToString() {
		// read pair 1 test.  Expected no reads in final pileup.
		int snpPos=23816120;
		Interval snpInterval = new Interval("1", snpPos, snpPos, true, "test");
		SNPGenomicBasePileUp p = new SNPGenomicBasePileUp(snpInterval);

		List<SAMRecord> recs = getReadsByName("H3FF2CCXX150427:1:2107:16631:1", this.smallBAMFile);
		p.addReads(recs);

		p.finalizePileup();
		Assert.assertNotNull(p.toString());
	}


	@Test
	public void testReadsMatchedBase() {
		// read pair 2 test.  Expected 1 base in pileup is A, at quality 30.
		int snpPos=23816120;
		Interval snpInterval = new Interval("1", snpPos, snpPos, true, "test");
		SNPGenomicBasePileUp p = new SNPGenomicBasePileUp(snpInterval);

		List<SAMRecord> recs = getReadsByName("H3FF2CCXX150427:1:2107:16631:2", this.smallBAMFile);
		p.addReads(recs);

		p.finalizePileup();
		List<Character> bases = p.getBasesAsCharacters();
		List<Character> expectedBases = new ArrayList<>();
		expectedBases.add(new Character ('T'));
		Assert.assertEquals(expectedBases, bases);

		List<Byte> quals = p.getQualities();
		List<Byte> expectedQuals = new ArrayList<>();
		expectedQuals.add(new Byte ((byte)30));
		Assert.assertEquals(expectedQuals, quals);

	}

	@Test
	public void testReadsUnpaired() {
		// read pair 2 test.  Expected 1 base in pileup is A, at quality 30.
		int snpPos=23816120;
		Interval snpInterval = new Interval("1", snpPos, snpPos, true, "test");
		SNPGenomicBasePileUp p = new SNPGenomicBasePileUp(snpInterval);

		List<SAMRecord> recs = getReadsByName("H3FF2CCXX150427:1:2107:16631:3", this.smallBAMFile);
		p.addReads(recs);
		recs = getReadsByName("H3FF2CCXX150427:1:2107:16631:4", this.smallBAMFile);
		p.addReads(recs);

		p.finalizePileup();
		List<Character> bases = p.getBasesAsCharacters();
		List<Character> expectedBases = new ArrayList<>();
		expectedBases.add(new Character ('T'));
		expectedBases.add(new Character ('A'));
		Assert.assertEquals(expectedBases, bases);

		List<Byte> quals = p.getQualities();
		List<Byte> expectedQuals = new ArrayList<>();
		expectedQuals.add(new Byte ((byte)28));
		expectedQuals.add(new Byte ((byte)28));
		Assert.assertEquals(expectedQuals, quals);

	}

	@Test
	public void testReadNotOverlappingSNP() {
		// read 5 test.  Doesn't overlap SNP, so no reads added to pileup.
		int snpPos=23816120;
		Interval snpInterval = new Interval("1", snpPos, snpPos, true, "test");
		SNPGenomicBasePileUp p = new SNPGenomicBasePileUp(snpInterval);

		List<SAMRecord> recs = getReadsByName("H3FF2CCXX150427:1:2107:16631:5", this.smallBAMFile);
		p.addReads(recs);

		p.finalizePileup();
		List<Character> bases = p.getBasesAsCharacters();
		List<Character> expectedBases = new ArrayList<>();
		Assert.assertEquals(expectedBases, bases);

		List<Byte> quals = p.getQualities();
		List<Byte> expectedQuals = new ArrayList<>();
		Assert.assertEquals(expectedQuals, quals);

	}

	@Test
	public void testGetFilteredPileupByBaseQuality () {
		int snpPos=23816120;
		Interval snpInterval = new Interval("1", snpPos, snpPos, true, "test");
		SNPGenomicBasePileUp p = new SNPGenomicBasePileUp(snpInterval);

		SamReader reader = SamReaderFactory.makeDefault().open(this.smallBAMFile);
		Iterator<SAMRecord> iter = reader.iterator();
		while (iter.hasNext())
			p.addRead(iter.next());

		p.finalizePileup();
		p=p.getFilteredPileupByBaseQuality((byte)30);
		int numBases = p.getNumBases();
		Assert.assertEquals(1, numBases);
		Assert.assertTrue('T'==p.getBasesAsCharacters().get(0));

	}

	private List<SAMRecord> getReadsByName (final String readName, final File bamFile) {
		List<SAMRecord> result = new ArrayList<>();
		SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
		Iterator<SAMRecord> iter = reader.iterator();
		while (iter.hasNext()) {
			SAMRecord r1 = iter.next();
			if (r1.getReadName().equals(readName))
				result.add(r1);
		}
		return result;
	}
}
