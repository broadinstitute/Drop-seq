package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;
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

    /**
     * getSingleFunction must agree with processRead for every read, and must not modify the read.
     */
    @Test
    public void testGetSingleFunctionMatchesProcessRead() {
        for (boolean allGenes : new boolean[]{false, true}) {
            GeneFunctionProcessor gfp = new GeneFunctionProcessor("mn", "ms", "mf", allGenes, StrandStrategy.SENSE,
                    List.of(LocusFunction.CODING), FunctionalDataProcessorStrategy.DROPSEQ);
            final SamReader in = SamReaderFactory.makeDefault().open(new File(BAM));
            for (SAMRecord r : in.iterator().stream().toList()) {
                final String before = r.getSAMString();
                final FunctionalData fd = gfp.getSingleFunction(r);
                Assert.assertEquals(r.getSAMString(), before, "getSingleFunction must not modify the read");
                final List<SAMRecord> processed = gfp.processRead(r.deepCopy());
                if (processed.size() == 1) {
                    Assert.assertNotNull(fd);
                    Assert.assertEquals(fd.getGene(), processed.get(0).getStringAttribute("mn"));
                    Assert.assertEquals(fd.getGeneStrand(), processed.get(0).getStringAttribute("ms"));
                } else {
                    Assert.assertNull(fd);
                }
            }
        }
    }
}
