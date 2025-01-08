package org.broadinstitute.dropseqrna.utils;


import htsjdk.samtools.util.Log;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;

public class CompareBAMTagValuesTest {

    private static final File TESTDATA_DIR = new File("testdata/org/broadinstitute/dropseq/utils");

    private static final PicardHtsPath STARsolo = new PicardHtsPath(TESTDATA_DIR + "/compare_tags_STARsolo.sam");
    private static final PicardHtsPath CELLRANGER = new PicardHtsPath(TESTDATA_DIR + "/compare_tags_CellRanger.sam");
    private static final File test1Expected = new File(TESTDATA_DIR + "/testCompareBAMTagValuesUnpairedReads_report.txt");

    // private static final PicardHtsPath STARsolo = new PicardHtsPath(TESTDATA_DIR+"/test1.sam");
    // private static final PicardHtsPath CELLRANGER = new PicardHtsPath(TESTDATA_DIR+"/test2.sam");

    private static final PicardHtsPath DROPSEQ = new PicardHtsPath(TESTDATA_DIR + "/compare_tags_DropSeq.sam");
    private final Log log = Log.getInstance(CompareBAMTagValuesTest.class);

    @Test
    // TEST 1.
    public void testCompareBAMTagValuesUnpairedReads() {
        CompareBAMTagValues c = new CompareBAMTagValues();
        c.INPUT_1 = Collections.singletonList(STARsolo);
        c.INPUT_2 = Collections.singletonList(CELLRANGER);
        c.TAGS_1 = Arrays.asList("CB", "CR");
        c.TAGS_2 = Arrays.asList("CB", "CR");
        c.BAM_OUTPUT_1 = TestUtils.getTempReportFile("compare_starsolo_cellranger_1.", ".bam");
        c.BAM_OUTPUT_2 = TestUtils.getTempReportFile("compare_starsolo_cellranger_2.", ".bam");
        c.TAG_VALUES_OUTPUT = TestUtils.getTempReportFile("compare_starsolo_cellranger", ".report");
        c.READ_COUNT_OUTPUT = TestUtils.getTempReportFile("compare_starsolo_cellranger", ".read_count");
        c.useFixedLengthBaseCompression = true;
        c.STRICT = false;
        c.BAM_OUTPUT_1.deleteOnExit();
        c.BAM_OUTPUT_2.deleteOnExit();
        c.TAG_VALUES_OUTPUT.deleteOnExit();
        c.READ_COUNT_OUTPUT.deleteOnExit();
        int result = c.doWork();
        log.info("Result file: " + c.TAG_VALUES_OUTPUT.getAbsolutePath());
        Assert.assertEquals(result, 0);
        boolean outputSame = TestUtils.testFilesSame(test1Expected, c.TAG_VALUES_OUTPUT);
        Assert.assertTrue(outputSame);
    }

    @Test
    // do not use fixed length base compression for the tags.
    public void testCompareBAMTagValuesUnpairedReads2() {
        CompareBAMTagValues c = new CompareBAMTagValues();
        c.INPUT_1 = Collections.singletonList(STARsolo);
        c.INPUT_2 = Collections.singletonList(CELLRANGER);
        c.TAGS_1 = Arrays.asList("CB", "CR");
        c.TAGS_2 = Arrays.asList("CB", "CR");
        c.BAM_OUTPUT_1 = TestUtils.getTempReportFile("compare_starsolo_cellranger_1.", ".bam");
        c.BAM_OUTPUT_2 = TestUtils.getTempReportFile("compare_starsolo_cellranger_2.", ".bam");
        c.TAG_VALUES_OUTPUT = TestUtils.getTempReportFile("compare_starsolo_cellranger", ".report");
        c.READ_COUNT_OUTPUT = TestUtils.getTempReportFile("compare_starsolo_cellranger", ".read_count");
        c.useFixedLengthBaseCompression = false;
        c.STRICT = false;
        c.BAM_OUTPUT_1.deleteOnExit();
        c.BAM_OUTPUT_2.deleteOnExit();
        c.TAG_VALUES_OUTPUT.deleteOnExit();
        c.READ_COUNT_OUTPUT.deleteOnExit();
        int result = c.doWork();
        log.info("Result file: " + c.TAG_VALUES_OUTPUT.getAbsolutePath());
        Assert.assertEquals(result, 0);
        boolean outputSame = TestUtils.testFilesSame(test1Expected, c.TAG_VALUES_OUTPUT);
        Assert.assertTrue(outputSame);
    }

    @Test
    public void testCompareBAMTagValuesPairedReads() {
        // test 2 - DROPSEQ is a paired end read version of STARsolo.
        CompareBAMTagValues c = new CompareBAMTagValues();
        c.INPUT_1 = Collections.singletonList(DROPSEQ);
        c.INPUT_2 = Collections.singletonList(CELLRANGER);
        c.TAGS_1 = Arrays.asList("CB", "CR");
        c.TAGS_2 = Arrays.asList("CB", "CR");
        c.BAM_OUTPUT_1 = TestUtils.getTempReportFile("compare_starsolo_dropseq_1.", ".bam");
        c.BAM_OUTPUT_2 = TestUtils.getTempReportFile("compare_starsolo_dropseq_2.", ".bam");
        c.TAG_VALUES_OUTPUT = TestUtils.getTempReportFile("compare_starsolo_dropseq", ".report");
        c.READ_COUNT_OUTPUT = TestUtils.getTempReportFile("compare_starsolo_dropseq", ".read_count");
        c.useFixedLengthBaseCompression = true;
        c.STRICT = false;
        c.BAM_OUTPUT_1.deleteOnExit();
        c.BAM_OUTPUT_2.deleteOnExit();
        c.TAG_VALUES_OUTPUT.deleteOnExit();
        c.READ_COUNT_OUTPUT.deleteOnExit();
        int result = c.doWork();
        log.info("Result file: " + c.TAG_VALUES_OUTPUT.getAbsolutePath());
        Assert.assertEquals(result, 0);
        boolean outputSame = TestUtils.testFilesSame(test1Expected, c.TAG_VALUES_OUTPUT);
        Assert.assertTrue(outputSame);

    }
}
