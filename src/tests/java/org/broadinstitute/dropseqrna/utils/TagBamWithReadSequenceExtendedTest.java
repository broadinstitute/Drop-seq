/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.List;

import junit.framework.Assert;

import org.testng.annotations.Test;

public class TagBamWithReadSequenceExtendedTest {

	@Test (enabled=true)
	public void testHardClip () {
		SAMRecord r = new SAMRecord(null);
		r.setReadString("AAAAAAAAAATTTTTT"); //10A6T
		byte [] quals = {10,10,10,0,0,0,10,10,10,10,0,0,0,0,10,10}; // 0 qual at 4-6, 11-14
		r.setBaseQualities(quals);
		List<BaseRange> ranges = new ArrayList<BaseRange>();
		ranges.add(new BaseRange (4,6));
		ranges.add(new BaseRange (11,14));

		r=TagBamWithReadSequenceExtended.hardClipBasesFromRead(r, ranges);

		String newSeq = r.getReadString();
		byte [] newQuals = r.getBaseQualities();
		String expectedSeq ="AAAAAAATT";
		byte [] expectedQuals = {10,10,10,10,10,10,10,10,10};

		Assert.assertEquals(expectedSeq.length(), newSeq.length());
		Assert.assertEquals(expectedSeq, newSeq);
		for (int i=0; i<expectedQuals.length;i++)
			Assert.assertEquals(expectedQuals[i], newQuals[i]);
	}
}
