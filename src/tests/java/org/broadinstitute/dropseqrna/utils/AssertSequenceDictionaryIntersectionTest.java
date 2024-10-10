package org.broadinstitute.dropseqrna.utils;

import java.io.File;
import java.util.List;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import org.testng.annotations.Test;

public class AssertSequenceDictionaryIntersectionTest {
	private static final Log log = Log.getInstance(AssertSequenceDictionaryIntersectionTest.class);
    private final File TESTDATA_DIR =
            new File("testdata/org/broadinstitute/dropseq/utils/SequenceDictionaryIntersectionTest");
    private final File BAM_FILE = new File(TESTDATA_DIR, "chr.sam");
    private final File VCF_FILE = new File(TESTDATA_DIR, "chr.vcf");
    private final File INTERVAL_LIST_FILE = new File(TESTDATA_DIR, "chr.interval_list");

    @Test
    public void testAssertIntersectionVcfBam() {
        AssertSequenceDictionaryIntersection.assertIntersectionVcfBam(VCF_FILE, BAM_FILE, log);
        AssertSequenceDictionaryIntersection.assertIntersectionVcfBam(VCF_FILE.toPath(), BAM_FILE.toPath(), null);
    }

    @Test
    public void testAssertIntersectionObjectBam() {
        final IntervalList intervalList = IntervalList.fromFile(INTERVAL_LIST_FILE);
        final String intervalListDescription = INTERVAL_LIST_FILE.getName();
        AssertSequenceDictionaryIntersection.assertIntersectionObjectBam(
                intervalList, intervalListDescription, BAM_FILE, log);
        AssertSequenceDictionaryIntersection.assertIntersectionObjectBam(
                intervalList, intervalListDescription, BAM_FILE.toPath(), null);
    }

    @Test
    public void testAssertIntersectionObjectVcf() {
        final IntervalList intervalList = IntervalList.fromFile(INTERVAL_LIST_FILE);
        final String intervalListDescription = INTERVAL_LIST_FILE.getName();
        AssertSequenceDictionaryIntersection.assertIntersectionObjectVcf(
                intervalList, intervalListDescription, VCF_FILE, log);
        AssertSequenceDictionaryIntersection.assertIntersectionObjectVcf(
                intervalList, intervalListDescription, VCF_FILE.toPath(), null);
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testAssertNoIntersection() {
        final SAMSequenceRecord samSequenceRecord1 = new SAMSequenceRecord("chr1", 1000);
        final SAMSequenceRecord samSequenceRecord2 = new SAMSequenceRecord("chr2", 2000);
        final SAMSequenceDictionary samSequenceDictionary1 = new SAMSequenceDictionary(List.of(samSequenceRecord1));
        final SAMSequenceDictionary samSequenceDictionary2 = new SAMSequenceDictionary(List.of(samSequenceRecord2));
        AssertSequenceDictionaryIntersection.assertIntersection(
                samSequenceDictionary1, "dict1", samSequenceDictionary2, "dict2", null);
    }
}
