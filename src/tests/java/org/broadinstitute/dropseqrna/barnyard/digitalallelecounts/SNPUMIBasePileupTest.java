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

import java.io.File;
import java.util.Iterator;
import java.util.List;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileup;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;

public class SNPUMIBasePileupTest {

	// see SNPUMICellReadProcessorTest for info on the hek_5_cell_2_snp_testdata BAM.
	// selected reads by:
	// samtools view -H hek_5_cell_2_snp_testdata.bam > smallTest_snpUMIPileUp.sam
	// samtools view -q 10 hek_5_cell_2_snp_testdata.bam HUMAN_1:76227022-76227022 |grep -v 50M |head -n 10 >> smallTest_snpUMIPileUp.sam
	// samtools view -q 10 hek_5_cell_2_snp_testdata.bam HUMAN_1:76227022-76227022 |grep 50M |head -n 5 >> smallTest_snpUMIPileUp.sam
	// java -jar ~/picard/dist/picard.jar SortSam I=smallTest_snpUMIPileUp.sam SO=coordinate O=smallTest_snpUMIPileUp.bam
	// IGV reports 6A, 7G.
	// G37, G37, G37, G37, G37, G32, G27, A37, A32, A27, A13, A8, A8
	// See smallTest_snpUMIPileUp.summary.txt for a listing of the bases/qualities/cell/molecular barcodes for all 15 reads

	private final File smallBAMFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/smallTest_snpUMIPileUp.sam");

	@Test(enabled=true)
	public void testAllReadsPileupError() {

	}
	@Test(enabled=true)
	public void testAllReadsPileup() {
		SamReader reader = SamReaderFactory.makeDefault().open(smallBAMFile);
		Iterator<SAMRecord> iter = reader.iterator();
		int snpPos=76227022;
		Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
		SNPUMIBasePileup p = new SNPUMIBasePileup(snpInterval, "ACADM", "fake_cell", "AAAAA");
		while (iter.hasNext())
			p.addRead(iter.next());

		List<Character> bases= p.getBasesAsCharacters();
		ObjectCounter<Character> baseCounter = new ObjectCounter<>();
		for (Character c: bases)
			baseCounter.increment(c);
		Assert.assertEquals(6, baseCounter.getCountForKey('A'));
		Assert.assertEquals(7, baseCounter.getCountForKey('G'));

		Assert.assertNotNull(p.toString());

	}


	@Test(enabled=true)
	public void singleReadTest1() {

		int snpPos=76227022;
		Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
		SNPUMIBasePileup p = new SNPUMIBasePileup(snpInterval, "ACADM", "fake_cell", "AAAAA");

		SAMRecord r1 = getReadByName("NS500217:67:H14GMBGXX:4:13611:8735:3829", this.smallBAMFile);
		p.addRead(r1);
		List<Character> bases = p.getBasesAsCharacters();
		Assert.assertTrue('A'==bases.get(0));
		Assert.assertTrue(8==p.getQualities().get(0));
	}

	@Test(enabled=true)
	public void singleReadTest2() {
		int snpPos=76227022;
		Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
		SNPUMIBasePileup p = new SNPUMIBasePileup(snpInterval, "ACADM", "fake_cell", "AAAAA");
		SAMRecord r1 = getReadByName("NS500217:67:H14GMBGXX:1:13306:23964:12839", this.smallBAMFile);

		p.addRead(r1);
		List<Character> bases = p.getBasesAsCharacters();
		Assert.assertTrue('G'==bases.get(0));
		Assert.assertTrue(32==p.getQualities().get(0));
	}

	private SAMRecord getReadByName (final String readName, final File bamFile) {
		SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
		Iterator<SAMRecord> iter = reader.iterator();
		while (iter.hasNext()) {
			SAMRecord r1 = iter.next();
			if (r1.getReadName().equals(readName))
				return (r1);
		}
		return null;
	}

	@Test(enabled=true)
	/**
	 * Add an uneven number of bases and quals to trip the exception throw.
	 */
	public void testAddBaseQualsError () {
		int snpPos=76227022;
		Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
		SNPUMIBasePileup p = new SNPUMIBasePileup(snpInterval, "ACADM", "fake_cell", "AAAAA");
		char [] bases = {'A', 'A'};
		byte [] quals = {27,17,55};
		byte [] bases2 = new byte [bases.length];
		StringUtil.charsToBytes(bases, 0, bases.length, bases2, 0);
		boolean passes=false;
		try {
			p.setBasesAndQualities(bases2, quals);
		} catch (IllegalArgumentException e) {
			Assert.assertNotNull(e);
			passes=true;
		}
		Assert.assertTrue(passes);
	}

	@Test(enabled=true)
	public void testAddBaseQuals () {
		int snpPos=76227022;
		Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
		SNPUMIBasePileup p = new SNPUMIBasePileup(snpInterval, "ACADM", "fake_cell", "AAAAA");
		char [] bases = {'A', 'A'};
		byte [] quals = {27,55};
		byte [] bases2 = new byte [bases.length];
		StringUtil.charsToBytes(bases, 0, bases.length, bases2, 0);
		boolean passes=true;
		try {
			p.setBasesAndQualities(bases2, quals);
		} catch (IllegalArgumentException e) {
			passes=false;
		}
		Assert.assertTrue(passes);
	}
}
