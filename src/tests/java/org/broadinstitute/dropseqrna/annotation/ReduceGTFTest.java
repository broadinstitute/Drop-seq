package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.OverlapDetector;

import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.annotation.GTFRecord;
import org.broadinstitute.dropseqrna.annotation.GeneFromGTF;
import org.broadinstitute.dropseqrna.annotation.ReduceGTF;
import org.testng.Assert;
import org.testng.annotations.Test;


public class ReduceGTFTest {

	File GTF_FILE1 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_ISG15.gtf.gz");
	File GTF_FILE2 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_ISG15_FAM41C.gtf.gz");
	File GTF_FILE3 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_AL592188.5.gtf.gz");
	// contains just the APITD1 gene
	File GTF_FILE4 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_APITD1.gtf.gz");
	// contains both the APITD1 gene and the APITD1-CORT genes, that occupy the same bounds, but have different sets of exons.
	File GTF_FILE5 = new File("testdata/org/broadinstitute/transcriptome/annotation/human_APITD1_both.gtf.gz");
	
	File SD = new File("testdata/org/broadinstitute/transcriptome/annotation/human_g1k_v37_decoy_50.dict");
	
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void test1() {
		ReduceGTF r = new ReduceGTF();
		Set<String> featureTypes = new HashSet<>(r.FEATURE_TYPE);
		Set<String> ignoredFunctionalTypes = new HashSet<>(r.IGNORE_FUNC_TYPE);
		
		List<GTFRecord> records = r.parseGTF (GTF_FILE1, SD, featureTypes, ignoredFunctionalTypes, true);
		Assert.assertNotNull(records);
		
	}
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void test2() {
		ReduceGTF r = new ReduceGTF();
        Set<String> featureTypes = new HashSet<>(r.FEATURE_TYPE);
        Set<String> ignoredFunctionalTypes = new HashSet<>(r.IGNORE_FUNC_TYPE);
		
		List<GTFRecord> records = r.parseGTF (GTF_FILE3, SD, featureTypes, ignoredFunctionalTypes, true);
		Assert.assertNotNull(records);
		EnhanceGTFRecords e = new EnhanceGTFRecords();
		records = e.enhanceGTFRecords(records);
		Assert.assertNotNull(records);
		
		for (GTFRecord a: records) {
			Assert.assertEquals(a.getStart(), 152252);
			Assert.assertEquals(a.getEnd(), 152312);
			Assert.assertEquals(a.getStrandAsString(), "+");
			
			if (a.getFeatureType().equals("gene")) {
				Assert.assertEquals(a.getTranscriptID(), null);
				Assert.assertEquals(a.getTranscriptName(), null);
				Assert.assertEquals(a.getTranscriptType(), "miRNA");
			} else {
				Assert.assertEquals(a.getTranscriptID(), "ENST00000577630");
				Assert.assertEquals(a.getTranscriptName(), "AL592188.5-201");
			}
			
		}
		
	}
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testAPITD1() {
		ReduceGTF r = new ReduceGTF();
        Set<String> featureTypes = new HashSet<>(r.FEATURE_TYPE);
        Set<String> ignoredFunctionalTypes = new HashSet<>(r.IGNORE_FUNC_TYPE);
		
		List<GTFRecord> records = r.parseGTF (GTF_FILE4, SD, featureTypes, ignoredFunctionalTypes, true);
		Assert.assertNotNull(records);
		// gunzip -c human_APITD1.gtf.gz | grep -v CDS |grep -v start_codon |grep -v stop_codon |wc -l
		Assert.assertEquals(records.size(),26);
		
		EnhanceGTFRecords e = new EnhanceGTFRecords();
		OverlapDetector<GeneFromGTF> od = e.getOverlapDetector(records);
		Assert.assertEquals(od.getAll().size(),1);
		
		records = e.enhanceGTFRecords(records);
		Assert.assertNotNull(records);
		
		
		
	}
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testAPITD1Complex() {
		ReduceGTF r = new ReduceGTF();
        Set<String> featureTypes = new HashSet<>(r.FEATURE_TYPE);
        Set<String> ignoredFunctionalTypes = new HashSet<>(r.IGNORE_FUNC_TYPE);
		
		List<GTFRecord> records = r.parseGTF (GTF_FILE5, SD, featureTypes, ignoredFunctionalTypes, true);
		Assert.assertNotNull(records);
		// gunzip -c human_APITD1.gtf.gz | grep -v CDS |grep -v start_codon |grep -v stop_codon |wc -l
		Assert.assertEquals(records.size(),42);
		
		EnhanceGTFRecords e = new EnhanceGTFRecords();
		OverlapDetector<GeneFromGTF> od = e.getOverlapDetector(records);
		// this fails because the two genes share the same bounds, and one replaces the other in the OverlapDetector class.
		Assert.assertEquals(od.getAll().size(),2);
		
		records = e.enhanceGTFRecords(records);
		Assert.assertNotNull(records);
		
		
	}
}
