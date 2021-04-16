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

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.DigitalAlleleCounts;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileup;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.statistics.BinomialStatistics;
import org.junit.Assert;
import org.testng.annotations.Test;

import java.util.Collection;

public class DigitalAlleleCountsTest {

	private int snpPos=76227022;
	private Interval snpInterval = new Interval("HUMAN_1", snpPos, snpPos, true, "test");
	private String gene="ACADM";
	private String cell="fake_cell";

	private DigitalAlleleCounts getTestDataSet1 (char refAllele) {
		DigitalAlleleCounts dac = new DigitalAlleleCounts(snpInterval, gene, cell, 10, refAllele, 'N');
		SNPUMIBasePileup p1 = buildPileUp(dac, "AAAAA", new char [] {'A', 'A'}, new byte [] {27,17});
		SNPUMIBasePileup p2 = buildPileUp(dac, "TTTTT", new char [] {'T', 'T'}, new byte [] {27,17});
		SNPUMIBasePileup p3 = buildPileUp(dac, "AATTT", new char [] {'T', 'T', 'T'}, new byte [] {27,17,11});
		SNPUMIBasePileup p4 = buildPileUp(dac, "GGGGG", new char [] {'A', 'T', 'T'}, new byte [] {27,17,11});
		SNPUMIBasePileup p5 = buildPileUp(dac, "CCCCC", new char [] {'A', 'A', 'T'}, new byte [] {27,17,4});
		dac.addPileup(p1);
		dac.addPileup(p2);
		dac.addPileup(p3);
		dac.addPileup(p4);
		dac.addPileup(p5);
		return (dac);
	}

	private void testUMIPurityFiltering () {

	}

	@Test(enabled=true)
	public void testUMIMeanPurity() {
		// AAAAA collapses with the other two.
		DigitalAlleleCounts dac = new DigitalAlleleCounts(snpInterval, gene, cell, 10, 'N','N');
		SNPUMIBasePileup p1 = buildPileUp(dac, "AAAAA", new char [] {'A', 'A','A','A', 'A'}, new byte [] {27,17,18,25,17});
		SNPUMIBasePileup p2 = buildPileUp(dac, "AAAAT", new char [] {'T', 'T','T', 'T'}, new byte [] {27,17,20,21});
		SNPUMIBasePileup p3 = buildPileUp(dac, "TAAAA", new char [] {'A','A','A'}, new byte [] {27,17,11});
		dac.addPileup(p1);
		dac.addPileup(p2);
		dac.addPileup(p3);

		double purityNoCollapse = dac.getMeanUMIPurity();
		Assert.assertEquals(1, purityNoCollapse, 0.001);
		DigitalAlleleCounts dacCollapsed=dac.collapseUMIs(1);
		double purityWithCollapse = dacCollapsed.getMeanUMIPurity();
		Assert.assertEquals(0.6666, purityWithCollapse, 0.001);


	}

	@Test(enabled=true)
	public void testUMICollapse () {
		// AAAAA absorbs AAAAT before AAATTT can get to it.
		DigitalAlleleCounts dac = new DigitalAlleleCounts(snpInterval, gene, cell, 10, 'N', 'N');
		SNPUMIBasePileup p1 = buildPileUp(dac, "AAAAA", new char [] {'A', 'A','A','A', 'A'}, new byte [] {27,17,18,25,17});
		SNPUMIBasePileup p2 = buildPileUp(dac, "AAATT", new char [] {'T', 'T','T', 'T'}, new byte [] {27,17,20,21});
		SNPUMIBasePileup p3 = buildPileUp(dac, "AAAAT", new char [] {'A','A','A'}, new byte [] {27,17,11});
		dac.addPileup(p1);
		dac.addPileup(p2);
		dac.addPileup(p3);

		dac=dac.collapseUMIs(1);
		Collection<String> umiList = dac.umis();
		for (String umi: umiList) {
			ObjectCounter<Character> result = dac.getReadsPerUMI(umi);
			if (umi.equals("AAAAA"))
				Assert.assertEquals(8, result.getCountForKey('A'));
			if (umi.equals("AAATT"))
				Assert.assertEquals(4, result.getCountForKey('T'));
		}
		// this UMI no longer exists because of collapse!
		ObjectCounter<Character> t = dac.getReadsPerUMI("AAAAT");
		Assert.assertNull("", t);
	}

	private SNPUMIBasePileup buildPileUp (final DigitalAlleleCounts dac, final String molecularBarcode, final char [] bases, final byte [] qualities) {
		SNPUMIBasePileup p = new SNPUMIBasePileup(dac.getSnpInterval(), dac.getGene(), dac.getCell(), molecularBarcode);
		byte [] bases2 = new byte [bases.length];
		StringUtil.charsToBytes(bases, 0, bases.length, bases2, 0);
		p.setBasesAndQualities(bases2, qualities);
		return (p);
	}


	@Test(enabled=true)
	public void getReadCounts() {
		DigitalAlleleCounts dac = getTestDataSet1 ('N');
		ObjectCounter<Character> result = dac.getReadCounts();
		Assert.assertEquals(5, result.getCountForKey('A'));
		Assert.assertEquals(7, result.getCountForKey('T'));
	}

	@Test(enabled=true)
	public void getReadsPerUMI() {
		DigitalAlleleCounts dac = getTestDataSet1 ('N');
		Collection<String> umiList = dac.umis();
		for (String umi: umiList) {
			ObjectCounter<Character> result = dac.getReadsPerUMI(umi);
			if (umi.equals("AAAAA"))
				Assert.assertEquals(2, result.getCountForKey('A'));
			if (umi.equals("TTTTT"))
				Assert.assertEquals(2, result.getCountForKey('T'));
			if (umi.equals("AATTT"))
				Assert.assertEquals(3, result.getCountForKey('T'));
			if (umi.equals("GGGGG")) {
				Assert.assertEquals(2, result.getCountForKey('T'));
				Assert.assertEquals(1, result.getCountForKey('A'));
			}
			if (umi.equals("CCCCC"))
				Assert.assertEquals(2, result.getCountForKey('A'));
				// no T, it's filtered out
		}

	}

	@Test(enabled=true)
	public void getUMIAlleleCount() {
		DigitalAlleleCounts dac = getTestDataSet1 ('N');
		ObjectCounter<Character> result = dac.getUMIAlleleCount();
		Assert.assertEquals(2, result.getCountForKey('A'));
		Assert.assertEquals(3, result.getCountForKey('T'));


	}

	@Test(enabled=true)
	public void getUMIPurity() {
		DigitalAlleleCounts dac = getTestDataSet1 ('N');
		Collection<String> umiList = dac.umis();
		for (String umi: umiList) {
			double purity = dac.getUMIPurity(umi);
			if (umi.equals("AAAAA"))
				Assert.assertEquals(1, purity, 0.01);
			if (umi.equals("TTTTT"))
				Assert.assertEquals(1, purity, 0.01);
			if (umi.equals("AATTT"))
				Assert.assertEquals(1, purity, 0.01);
			if (umi.equals("GGGGG"))
				Assert.assertEquals(0.66, purity, 0.01);
			if (umi.equals("CCCCC"))
				Assert.assertEquals(1, purity, 0.01);
				// no T, it's filtered out
		}

	}

	@Test
	public void umis() {
		DigitalAlleleCounts dac = getTestDataSet1 ('N');
		Collection<String> umiList = dac.umis();

		String [] umisExpected = new String [] {"AAAAA", "TTTTT", "AATTT", "GGGGG", "CCCCC"};
		int counter=0;
		for (String s: umisExpected)
			if (umiList.contains(s))
				counter++;
		Assert.assertTrue(counter==umisExpected.length);
	}

	
	@Test
	public void testMostCommonCounts () {
		DigitalAlleleCounts dac = getTestDataSet1 ('A');
		
		BinomialStatistics readsA = dac.getReadConfidenceInterval(0.95);
		Assert.assertEquals(0.416, readsA.getRatio(), 0.001);
		
		dac = getTestDataSet1 ('T');
		BinomialStatistics readsT = dac.getReadConfidenceInterval(0.95);
		Assert.assertEquals(0.583, readsT.getRatio(), 0.001);
	}
	

}
