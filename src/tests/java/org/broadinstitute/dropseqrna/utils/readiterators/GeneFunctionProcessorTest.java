package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorStrategy;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.annotation.LocusFunction;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.*;

public class GeneFunctionProcessorTest {

    private static final String BAM = "testdata/org/broadinstitute/dropseq/utils/readiterators/SPANXB1.bam";

    @Test
    /**
     * A specific test for when a metagene is comprised of two genes on opposite strands and the pair of reads
     * correctly match the strands.
     * This was being filtered out due to a bug, this test is to make sure that bug is not reintroduced.
     */
    public void testMetaGeneTwoStrands() {
        GeneFunctionProcessor gfp = new GeneFunctionProcessor("mn", "ms", "mf", false, StrandStrategy.SENSE,
                List.of(LocusFunction.CODING), FunctionalDataProcessorStrategy.DROPSEQ);

        final SamReader in = SamReaderFactory.makeDefault().open(new File(BAM));
        List<SAMRecord> samRecords = in.iterator().stream().toList();
        SAMRecord result = gfp.processReads(samRecords);
        Assert.assertNotNull(result);
        Assert.assertEquals(result.getReadName(),"LH00118:69:22CNFNLT3:3:2141:14906:8787");

    }
}