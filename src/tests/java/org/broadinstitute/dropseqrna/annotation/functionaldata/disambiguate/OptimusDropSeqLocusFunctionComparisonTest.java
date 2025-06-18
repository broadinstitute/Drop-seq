package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Collections;

import static org.testng.Assert.*;

public class OptimusDropSeqLocusFunctionComparisonTest {

    private final String ROOT_DIR="testdata/org/broadinstitute/dropseq/annotation/functionaldata/disambiguate/";

    private final String bamfile = ROOT_DIR+"testdata.bam";

    private final File EXPECTED_OUTPUT = new File(ROOT_DIR, "testdata.function_comparison_ambiguous_umis.txt");
    private final File EXPECTED_SUMMARY = new File(ROOT_DIR, "testdata.function_comparison_summary.txt");
    private final File EXPECTED_CONFUSION_MATRIX =  new File(ROOT_DIR, "testdata.function_comparison_confusion_matrix.txt");

    private final File EXPECTED_OUT_BAM = new File(ROOT_DIR, "testdata.function_comparison_error_reads.bam");

    @Test
    /**
     * Writing this mostly to satisfy code coverage, but this does cover at least a little bit of most cases / functional types
     * so if I break something later this may find it - but I trust the other more explicit unit tests written to test
     * the functional types and validation more.
     */
    public void testDoWork()  {
        OptimusDropSeqLocusFunctionComparison cmd = new OptimusDropSeqLocusFunctionComparison();
        cmd.INPUT= Collections.singletonList(new File(this.bamfile));
        cmd.OUTPUT=TestUtils.getTempReportFile("OptimusDropSeqLocusFunctionComparisonTest.", ".function_comparison_ambiguous_umis.txt");
        cmd.SUMMARY=TestUtils.getTempReportFile("OptimusDropSeqLocusFunctionComparisonTest.", ".function_comparison_summary.txt");
        cmd.OUT_CONFUSION_MATRIX=TestUtils.getTempReportFile("OptimusDropSeqLocusFunctionComparisonTest.", ".function_comparison_confusion_matrix.txt");
        cmd.OUT_BAM=TestUtils.getTempReportFile("OptimusDropSeqLocusFunctionComparisonTest.", ".function_comparison_error_reads.bam");

        int result = cmd.doWork();
        Assert.assertEquals(result, 0);

        Assert.assertTrue(TestUtils.testFilesSame(this.EXPECTED_OUTPUT, cmd.OUTPUT));
        Assert.assertTrue(TestUtils.testFilesSame(this.EXPECTED_CONFUSION_MATRIX, cmd.OUT_CONFUSION_MATRIX, false));
        Assert.assertTrue(TestUtils.testFilesSame(this.EXPECTED_SUMMARY, cmd.SUMMARY));
        TestUtils.assertSamFilesSame(this.EXPECTED_OUT_BAM, cmd.OUT_BAM, false);


//        cmd.OUTPUT=/downloads/DGE/2023-08-25_v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1.retagged_pct50.function_comparison_ambiguous_umis.txt
//        SUMMARY=/downloads/DGE/2023-08-25_v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1.retagged_pct50.function_comparison_summary.txt
//        OUT_CONFUSION_MATRIX=/downloads/DGE/2023-08-25_v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1.retagged_pct50.function_comparison_confusion_matrix.txt
//        OUT_BAM=/downloads/DGE/2023-08-25_v3_Bamboo_18d_10X_RNAseq_Optiprep8000_CaH_17k_rxn1.retagged_pct50.function_comparison_error_reads.bam
//
            //TestUtils.assertSamRecordsSame();
    }

    @Test
    public void testCustomCommandLineValidation() {
    }
}