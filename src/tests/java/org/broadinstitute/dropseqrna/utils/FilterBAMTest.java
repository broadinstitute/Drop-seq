package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.dropseqrna.utils.FilterBAM;
import org.junit.Assert;
import org.testng.annotations.Test;

public class FilterBAMTest {

	@Test(enabled=true, groups = { "dropseq","transcriptome" })
	public void testNonPrimary() {
		FilterBAM b = new FilterBAM();
		b.RETAIN_ONLY_PRIMARY_READS=true;
		
		SAMRecord r = new SAMRecord(null);
		r.setNotPrimaryAlignmentFlag(true);
		boolean result = b.filterRead(r);
		Assert.assertTrue(result);
		
		r = new SAMRecord(null);
		r.setNotPrimaryAlignmentFlag(false);
		result = b.rejectNonPrimaryReads(r);
		Assert.assertFalse(result);
		
		b.RETAIN_ONLY_PRIMARY_READS=false;
		r.setNotPrimaryAlignmentFlag(true);
		result = b.filterRead(r);
		Assert.assertFalse(result);
	}
	
	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testMapQuality () {
		FilterBAM b = new FilterBAM();
		b.MINIMUM_MAPPING_QUALITY=10;
		SAMRecord r = new SAMRecord(null);
		r.setMappingQuality(12);
		Assert.assertFalse(b.filterRead(r));
		
		r.setMappingQuality(9);
		Assert.assertTrue(b.filterRead(r));
		
		b.MINIMUM_MAPPING_QUALITY=null;
		Assert.assertFalse(b.filterRead(r));
	}
	
	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testSoftMatch () {
		FilterBAM b = new FilterBAM();
		List<String> retained = new ArrayList<String>();
		retained.add("HUMAN");
		b.REF_SOFT_MATCHED_RETAINED=retained;
		b.buildPatterns();
		SAMRecord r = new SAMRecord(null);
		r.setReferenceName("HUMAN_CHR1");
		Assert.assertFalse(b.filterRead(r));
		r.setReferenceName("CHIMP_CHR1");
		Assert.assertTrue(b.filterRead(r));
		
		retained.add("CHIMP");
		b.REF_SOFT_MATCHED_RETAINED=retained;
		b.buildPatterns();
		Assert.assertFalse(b.filterRead(r));
		
		
	}
	
	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testExactMatch () {
		FilterBAM b = new FilterBAM();
		List<String> retained = new ArrayList<String>();
		retained.add("HUMAN_CHR1");
		b.REF_HARD_MATCHED_RETAINED=retained;
		b.buildPatterns();
		SAMRecord r = new SAMRecord(null);
		r.setReferenceName("HUMAN_CHR1");
		Assert.assertFalse(b.filterRead(r));
		r.setReferenceName("CHIMP_CHR1");
		Assert.assertTrue(b.filterRead(r));
		
		
	}
	
	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testRejectOnTags () {
		FilterBAM b = new FilterBAM();
		b.TAG_REJECT_COMBINE_FLAG="UNION";
		List<String> tags = new ArrayList<String>();
		tags.add("XE");
		tags.add("XZ");
		b.TAG_REJECT=tags;
		SAMRecord r = new SAMRecord(null);
		r.setAttribute("XE", "Foo");
		boolean result = b.filterRead(r);
		Assert.assertTrue(result);
		// if you need to see 2 tags and only see 1, then you don't reject the read.
		b.TAG_REJECT_COMBINE_FLAG="INTERSECT";
		result = b.rejectOnTags(tags, r);
		Assert.assertFalse(result);
		
	}
	
}
