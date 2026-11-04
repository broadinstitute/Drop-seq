package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.dropseqrna.utils.CompareBAMTagValues;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import picard.nio.PicardHtsPath;


public class TagReadWithGeneFunctionTest {
	File testBAMFile= new File ("testdata/org/broadinstitute/dropseq/metrics/NucBYReg4Reg.MOUSE.GCTAAGTAAGAT.Elp2.fixed.bam");
	File annotationsFile=new File ("testdata/org/broadinstitute/dropseq/metrics/mm10_Elp2.gtf");
	File OUT_BAM = new File ("testdata/org/broadinstitute/dropseq/metrics/NucBYReg4Reg.MOUSE.GCTAAGTAAGAT.Elp2.gene_function_tagged.bam");

	@Test
	public void testDoWork() throws IOException {
		TagReadWithGeneFunction t = new TagReadWithGeneFunction();
		File tempBAM = File.createTempFile("TagReadWithGeneFunctionTest", ".bam");
		tempBAM.deleteOnExit();
		File tempSummary=File.createTempFile("TagReadWithGeneFunctionTest", ".summary");
		tempSummary.deleteOnExit();

		t.INPUT=new PicardHtsPath(testBAMFile);
		t.OUTPUT=tempBAM;
		t.ANNOTATIONS_FILE=annotationsFile;
		t.SUMMARY=tempSummary;
		int returnVal = t.doWork();
		Assert.assertTrue(returnVal==0);

		// test output BAM
		List<String> tags = new ArrayList<>(Arrays.asList("XC", "gn", "gs", "gf", "XF"));
		compareBAMTagValues(OUT_BAM, tempBAM, tags, 0);


	}

	private int writeRepeatedBAM (final File out, final int repeats) {
		SamReader inputSam = SamReaderFactory.makeDefault().open(testBAMFile);
		List<SAMRecord> reads = new ArrayList<>();
		for (SAMRecord r: inputSam)
			reads.add(r);
		// the repeated reads are not in coordinate order, so let the writer sort them to match the header.
		SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(inputSam.getFileHeader(), false, out);
		int count=0;
		for (int i=0; i<repeats; i++)
			for (SAMRecord r: reads) {
				SAMRecord copy = r.deepCopy();
				copy.setReadName(r.getReadName() + "." + i);
				// every 7th read is unmapped, and every 5th is on the other strand.
				if (count % 7 == 0) {
					copy.setReadUnmappedFlag(true);
					copy.setMappingQuality(0);
					copy.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
				} else if (count % 5 == 0)
					copy.setReadNegativeStrandFlag(!copy.getReadNegativeStrandFlag());
				writer.addAlignment(copy);
				count++;
			}
		writer.close();
		return count;
	}

	private List<String> readAll (final File bam) {
		List<String> result = new ArrayList<>();
		for (SAMRecord r: SamReaderFactory.makeDefault().open(bam))
			result.add(r.getSAMString());
		return result;
	}

	private void runTagger (final File input, final File output, final File summary, final int numThreads, final int batchSize) {
		TagReadWithGeneFunction t = new TagReadWithGeneFunction();
		t.INPUT=new PicardHtsPath(input);
		t.OUTPUT=output;
		t.ANNOTATIONS_FILE=annotationsFile;
		t.SUMMARY=summary;
		t.NUM_THREADS=numThreads;
		t.BATCH_SIZE=batchSize;
		Assert.assertNull(t.customCommandLineValidation());
		Assert.assertEquals(t.doWork(), 0);
	}

	@Test
	// the output (every read, in the same order) and the summary must not depend on the number of threads.
	public void testNumThreadsDoesNotChangeOutput () throws IOException {
		File input = File.createTempFile("TagReadWithGeneFunctionTest.input", ".bam");
		input.deleteOnExit();
		int numReads = writeRepeatedBAM(input, 50);

		File expectedBAM = File.createTempFile("TagReadWithGeneFunctionTest.expected", ".bam");
		expectedBAM.deleteOnExit();
		File expectedSummary = File.createTempFile("TagReadWithGeneFunctionTest.expected", ".summary");
		expectedSummary.deleteOnExit();
		runTagger(input, expectedBAM, expectedSummary, 1, 1000);
		List<String> expected = readAll(expectedBAM);
		Assert.assertEquals(expected.size(), numReads);
		byte[] expectedSummaryBytes = Files.readAllBytes(expectedSummary.toPath());

		// a batch size of 7 gives hundreds of batches, none of which line up with the unmapped/strand pattern above.
		for (int numThreads: new int [] {2, 4}) {
			for (int batchSize: new int[] {1, 7, 1000}) {
				File bam = File.createTempFile("TagReadWithGeneFunctionTest.threaded", ".bam");
				bam.deleteOnExit();
				File summary = File.createTempFile("TagReadWithGeneFunctionTest.threaded", ".summary");
				summary.deleteOnExit();
				runTagger(input, bam, summary, numThreads, batchSize);
				Assert.assertEquals(readAll(bam), expected, "threads=" + numThreads + " batchSize=" + batchSize);
				Assert.assertEquals(Files.readAllBytes(summary.toPath()), expectedSummaryBytes, "summary threads=" + numThreads + " batchSize=" + batchSize);
			}
		}
	}

	@Test
	public void testInvalidThreadOptions () {
		TagReadWithGeneFunction t = new TagReadWithGeneFunction();
		Assert.assertNull(t.customCommandLineValidation());
		t.NUM_THREADS=0;
		Assert.assertNotNull(t.customCommandLineValidation());
		t.NUM_THREADS=2;
		t.BATCH_SIZE=0;
		Assert.assertNotNull(t.customCommandLineValidation());
	}

	private void compareBAMTagValues(File input1, File input2, List<String> tags, int expectedProgramValue) {
		CompareBAMTagValues cbtv = new CompareBAMTagValues();
		cbtv.INPUT_1 = Collections.singletonList(new PicardHtsPath(input1));
		cbtv.INPUT_2 = Collections.singletonList(new PicardHtsPath(input2));
		cbtv.TAGS_1 = tags;
		cbtv.TAGS_2 = tags;
		cbtv.STRICT = true;
		int result = cbtv.doWork();
		Assert.assertTrue(result == expectedProgramValue);
	}
}
