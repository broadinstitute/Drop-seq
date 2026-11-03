package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.dropseqrna.annotation.GeneAnnotationReader;
import org.broadinstitute.dropseqrna.annotation.GeneFromGTF;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.OverlapDetector;
import picard.annotation.Gene;
import picard.annotation.LocusFunction;

/**
 * Tests the per-read tagging logic on its own, without the TagReadWithGeneFunction command.
 * TagReadWithGeneFunctionTest has the end to end tests.
 */
public class ReadGeneFunctionTaggerTest {

	File testBAMFile= new File ("testdata/org/broadinstitute/dropseq/metrics/NucBYReg4Reg.MOUSE.GCTAAGTAAGAT.Elp2.fixed.bam");
	File annotationsFile=new File ("testdata/org/broadinstitute/dropseq/metrics/mm10_Elp2.gtf");

	/** A tagger with the command's default tag names and settings. */
	private ReadGeneFunctionTagger newTagger (final OverlapDetector<Gene> geneOverlapDetector) {
		return new ReadGeneFunctionTagger(geneOverlapDetector, "gn", "gs", "gf", "XF", true, 0);
	}

	/**
	 * Gene "A" on the given strand, one transcript, exons 1-10 (coding) and 91-100 (UTR).
	 */
	private OverlapDetector<Gene> detectorWithGeneA(final boolean negativeStrand) {
		GeneFromGTF gene = new GeneFromGTF("18", 1, 100, negativeStrand, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> detector = new OverlapDetector<>(0, 0);
		detector.addLhs(gene, gene);
		return detector;
	}

	private ReadGeneFunctionTagger tagger(final OverlapDetector<Gene> detector, final boolean useStrandInfo) {
		return new ReadGeneFunctionTagger(detector, "gn", "gs", "gf", "XF", useStrandInfo, 0);
	}

	private SAMRecord read(final int start, final int end, final boolean negativeStrand) {
		SamReader inputSam = SamReaderFactory.makeDefault().open(testBAMFile);
		SAMRecord r = inputSam.iterator().next();
		r.setReferenceName("18");
		r.setAlignmentStart(start);
		int length = (end - start) + 1;
		r.setCigarString(length + "M");
		r.setReadNegativeStrandFlag(negativeStrand);
		r.setReadBases(Arrays.copyOf(r.getReadBases(), length));
		r.setBaseQualities(Arrays.copyOf(r.getBaseQualities(), length));
		return r;
	}

	@Test
	public void testCustomTagNames() {
		ReadGeneFunctionTagger t = new ReadGeneFunctionTagger(detectorWithGeneA(false), "N1", "S1", "F1", "R1", true, 0);
		SAMRecord r = t.tag(read(1, 10, false));
		Assert.assertEquals(r.getStringAttribute("N1"), "A");
		Assert.assertEquals(r.getStringAttribute("S1"), "+");
		Assert.assertEquals(r.getStringAttribute("F1"), "CODING");
		Assert.assertEquals(r.getStringAttribute("R1"), "CODING");
		Assert.assertNull(r.getStringAttribute("gn"));
	}

	@Test
	public void testStrandMetrics() {
		ReadGeneFunctionTagger t = tagger(detectorWithGeneA(false), true);
		// right strand, wrong strand, and a read that overlaps no gene (counted as right strand: nothing to be wrong about).
		t.tag(read(1, 10, false));
		SAMRecord wrongStrand = t.tag(read(1, 10, true));
		t.tag(read(200, 210, false));

		TagReadWithGeneFunction.ReadTaggingMetric m = t.getMetrics();
		Assert.assertEquals(m.TOTAL_READS, 3);
		Assert.assertEquals(m.READS_RIGHT_STRAND, 2);
		Assert.assertEquals(m.READS_WRONG_STRAND, 1);
		// the gene tags are strand independent, but the summary function only uses genes on the read's strand.
		Assert.assertEquals(wrongStrand.getStringAttribute("gn"), "A");
		Assert.assertEquals(wrongStrand.getStringAttribute("XF"), LocusFunction.INTERGENIC.name());
	}

	@Test
	public void testNoStrandMetricsWithoutStrandInfo() {
		ReadGeneFunctionTagger t = tagger(detectorWithGeneA(false), false);
		SAMRecord r = t.tag(read(1, 10, true));
		Assert.assertEquals(t.getMetrics().TOTAL_READS, 0);
		// without strand info the wrong strand gene still counts.
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.CODING.name());
	}

	@Test
	public void testTaggersDoNotShareMetrics() {
		OverlapDetector<Gene> detector = detectorWithGeneA(false);
		ReadGeneFunctionTagger t1 = tagger(detector, true);
		ReadGeneFunctionTagger t2 = tagger(detector, true);
		t1.tag(read(1, 10, false));
		t1.tag(read(1, 10, true));
		t2.tag(read(1, 10, false));
		Assert.assertEquals(t1.getMetrics().TOTAL_READS, 2);
		Assert.assertEquals(t2.getMetrics().TOTAL_READS, 1);
	}

	@Test
	public void testMergeMetrics() {
		OverlapDetector<Gene> detector = detectorWithGeneA(false);
		ReadGeneFunctionTagger t1 = tagger(detector, true);
		ReadGeneFunctionTagger t2 = tagger(detector, true);
		t1.tag(read(1, 10, false));
		t1.tag(read(1, 10, true));
		t2.tag(read(1, 10, true));
		t2.tag(read(200, 210, false));

		TagReadWithGeneFunction.ReadTaggingMetric total = new TagReadWithGeneFunction.ReadTaggingMetric();
		total.merge(t1.getMetrics());
		total.merge(t2.getMetrics());
		Assert.assertEquals(total.TOTAL_READS, 4);
		Assert.assertEquals(total.READS_RIGHT_STRAND, 2);
		Assert.assertEquals(total.READS_WRONG_STRAND, 2);
		Assert.assertEquals(total.READ_AMBIGUOUS_GENE_FIXED, 0);
		Assert.assertEquals(total.AMBIGUOUS_READS_REJECTED, 0);
	}

	@Test
	public void testCompoundStrings() {
		GeneFromGTF a = new GeneFromGTF("18", 1, 100, false, "A", "coding", "A", "coding", 1);
		GeneFromGTF b = new GeneFromGTF("18", 1, 100, true, "B", "coding", "B", "coding", 1);
		Assert.assertEquals(ReadGeneFunctionTagger.getCompoundName(Arrays.asList(a, b)), "A,B");
		Assert.assertEquals(ReadGeneFunctionTagger.getCompoundStrand(Arrays.asList(a, b)), "+,-");
		Assert.assertEquals(ReadGeneFunctionTagger.getCompoundFunctionName(
				Arrays.asList(LocusFunction.CODING, LocusFunction.INTRONIC)), "CODING,INTRONIC");
		Assert.assertNull(ReadGeneFunctionTagger.getCompoundName(Arrays.<Gene>asList()));
		Assert.assertNull(ReadGeneFunctionTagger.getCompoundStrand(Arrays.<Gene>asList()));
		Assert.assertNull(ReadGeneFunctionTagger.getCompoundFunctionName(Arrays.<LocusFunction>asList()));
	}

	@Test
	public void testTagIntronRead () {
		SamReader inputSam = SamReaderFactory.makeDefault().open(testBAMFile);
		SAMFileHeader header = inputSam.getFileHeader();
		SAMSequenceDictionary bamDict = header.getSequenceDictionary();
		final OverlapDetector<Gene> geneOverlapDetector = GeneAnnotationReader.loadAnnotationsFile(annotationsFile, bamDict);

		String recName="H575CBGXY:4:23606:1714:16102";
		SAMRecord r = getRecord(testBAMFile, recName);

		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		r=tagger.tag(r);

		String geneName = r.getStringAttribute("gn");
		String geneStrand = r.getStringAttribute("gs");
		String function = r.getStringAttribute("gf");

		Assert.assertEquals(geneName, "Elp2");
		Assert.assertEquals(geneStrand, "+");
		Assert.assertEquals(function, "INTRONIC");
		Assert.assertEquals(r.getStringAttribute("XF"), "INTRONIC");

	}

	@Test
	public void testTagCodingRead () {
		SamReader inputSam = SamReaderFactory.makeDefault().open(testBAMFile);
		SAMFileHeader header = inputSam.getFileHeader();
		SAMSequenceDictionary bamDict = header.getSequenceDictionary();
		final OverlapDetector<Gene> geneOverlapDetector = GeneAnnotationReader.loadAnnotationsFile(annotationsFile, bamDict);

		String recName="H575CBGXY:3:13502:3588:9959";
		SAMRecord r = getRecord(testBAMFile, recName);

		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		r=tagger.tag(r);

		Assert.assertEquals(r.getStringAttribute("gn"), "Elp2");
		Assert.assertEquals(r.getStringAttribute("gs"), "+");
		Assert.assertEquals(r.getStringAttribute("gf"), "CODING");
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.CODING.name());
	}

	//H575CBGXY:2:21206:5655:8103
	@Test
	public void testTagCodingUTR () {
		SamReader inputSam = SamReaderFactory.makeDefault().open(testBAMFile);
		SAMFileHeader header = inputSam.getFileHeader();
		SAMSequenceDictionary bamDict = header.getSequenceDictionary();
		final OverlapDetector<Gene> geneOverlapDetector = GeneAnnotationReader.loadAnnotationsFile(annotationsFile, bamDict);

		String recName="H575CBGXY:2:21206:5655:8103";
		SAMRecord r = getRecord(testBAMFile, recName);

		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		r=tagger.tag(r);

		String geneName = r.getStringAttribute("gn");
		String geneStrand = r.getStringAttribute("gs");
		String function = r.getStringAttribute("gf");

		Assert.assertEquals(geneName, "Elp2");
		Assert.assertEquals(geneStrand, "+");
		Assert.assertEquals(function, "UTR");
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.UTR.name());
	}

	@Test
	// One intergenic read (no other gene models overlapping)
	public void testIntergeicRead () {
		SAMRecord r = getFakeRecord(testBAMFile, 200, 210, false);
		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, true, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);
		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		r=tagger.tag(r);
		Assert.assertEquals(r.getStringAttribute("gn"), null);
		Assert.assertEquals(r.getStringAttribute("gs"), null);
		Assert.assertEquals(r.getStringAttribute("gf"), null);
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.INTERGENIC.name());
	}

	@Test
	// A read that is partly UTR exon and partly past the end of the gene: the gene keeps both UTR and INTERGENIC annotations.
	public void testPartiallyIntergenicRead () {
		SAMRecord r = getFakeRecord(testBAMFile, 95, 105, false);
		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);
		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		r=tagger.tag(r);
		Assert.assertEquals(r.getStringAttribute("gn"), "A,A");
		Assert.assertEquals(r.getStringAttribute("gs"), "+,+");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.UTR.name() + "," + LocusFunction.INTERGENIC.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.UTR.name());
	}

	@Test
	// A read that is partly intronic and partly before the start of the gene.
	public void testPartiallyIntergenicIntronicRead () {
		SAMRecord r = getFakeRecord(testBAMFile, 20, 40, false);
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 31, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 31, 100, 31, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(31, 35);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);
		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		r=tagger.tag(r);
		// bases 20-30 intergenic, 31-35 coding, 36-40 intronic.
		Assert.assertEquals(r.getStringAttribute("gn"), "A,A,A");
		Assert.assertEquals(r.getStringAttribute("gf"), "CODING,INTRONIC,INTERGENIC");
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.CODING.name());
	}

	@Test
	// A split read where one alignment block is in the gene and the other is intergenic: no gene is common to both blocks.
	public void testSplitReadOneBlockIntergenic () {
		SAMRecord r = getFakeSplitRecord(testBAMFile, 5, 10, 200, 210, false);
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);
		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		r=tagger.tag(r);
		Assert.assertEquals(r.getStringAttribute("gn"), null);
		Assert.assertEquals(r.getStringAttribute("gf"), null);
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.INTERGENIC.name());
	}

	@Test
	// 	2) One intronic read (no other gene models overlapping)
	public void testIntronicRead () {
		SAMRecord r = getFakeRecord(testBAMFile, 50, 75, false);
		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);
		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		r=tagger.tag(r);
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.INTRONIC.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.INTRONIC.name());
	}

	@Test
	// 	One UTR read (correct strand, no other gene models overlapping)
	public void testUTRRead () {
		SAMRecord r = getFakeRecord(testBAMFile, 91, 95, false);
		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);
		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		r=tagger.tag(r);
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.UTR.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.UTR.name());
	}

	@Test
	// 	One UTR read (wrong strand, no other gene models overlapping)
	public void testUTRReadWrongStrand () {
		SAMRecord r = getFakeRecord(testBAMFile, 91, 95, true);
		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);
		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		r=tagger.tag(r);
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.UTR.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.INTERGENIC.name());
	}


	@Test
	// One CODING read (correct strand, no other gene models overlapping)
	public void testCodingRead () {
		SAMRecord r = getFakeRecord(testBAMFile, 2, 8, false);
		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);
		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		r=tagger.tag(r);
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.CODING.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.CODING.name());
	}

	@Test
	// One CODING read (wrong strand, no other gene models overlapping)
	public void testCodingReadWrongStrand () {
		SAMRecord r = getFakeRecord(testBAMFile, 2, 8, true);
		boolean negStrandFlag = r.getReadNegativeStrandFlag();
		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);
		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		r=tagger.tag(r);
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.CODING.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.INTERGENIC.name());
	}

	// One intronic read (correct strand) that's UTR (wrong strand) on another gene
	@Test
	public void testIntronicCorrectUTRWrong () {
		// read on the negative strand
		SAMRecord r = getFakeRecord(testBAMFile, 91, 100, true);

		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);

		// gene with 2 exons, 1 coding from 50-60, 1 coding from 150-160. Negative strand gene.
		GeneFromGTF gene2 = new GeneFromGTF(r.getContig(), 50, 160, true, "B", "coding", "B", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx2 = gene2.addTranscript("trans2", 50, 160, 50, 150, 2, "trans2", "trans2", "coding");
		tx2.addExon(50, 60);
		tx2.addExon(150, 160);
		geneOverlapDetector.addLhs(gene2, gene2);

		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, ReadGeneFunctionTagger.GENE_NAME_COMPARATOR);

		r=tagger.tag(r);
		String gn = r.getStringAttribute("gn");
		String gs = r.getStringAttribute("gs");
		String gf = r.getStringAttribute("gf");

		// names always come out alphabetically sorted.
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName()+","+gene2.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+,-");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.UTR.name() + "," + LocusFunction.INTRONIC.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.INTRONIC.name());
	}


	// One intronic read (correct strand) that's CODING (wrong strand) on another gene
	@Test
	public void testIntronicCorrectCodingWrong () {
		// read on the negative strand
		SAMRecord r = getFakeRecord(testBAMFile, 91, 100, true);

		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 100, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);

		// gene with 2 exons, 1 coding from 50-60, 1 coding from 150-160. Negative strand gene.
		GeneFromGTF gene2 = new GeneFromGTF(r.getContig(), 50, 160, true, "B", "coding", "B", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx2 = gene2.addTranscript("trans2", 50, 160, 50, 150, 2, "trans2", "trans2", "coding");
		tx2.addExon(50, 60);
		tx2.addExon(150, 160);
		geneOverlapDetector.addLhs(gene2, gene2);

		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, ReadGeneFunctionTagger.GENE_NAME_COMPARATOR);

		r=tagger.tag(r);
		String gn = r.getStringAttribute("gn");
		String gs = r.getStringAttribute("gs");
		String gf = r.getStringAttribute("gf");

		// names always come out alphabetically sorted.
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName()+","+gene2.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+,-");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.CODING.name() + "," + LocusFunction.INTRONIC.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.INTRONIC.name());
	}

	// One UTR read (correct strand) that's intronic (wrong strand) on another gene
	@Test
	public void testUTRCorrectIntronicWrong () {
		// read on the negative strand
		SAMRecord r = getFakeRecord(testBAMFile, 91, 100, false);

		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);

		// gene with 2 exons, 1 coding from 50-60, 1 coding from 150-160. Negative strand gene.
		GeneFromGTF gene2 = new GeneFromGTF(r.getContig(), 50, 160, true, "B", "coding", "B", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx2 = gene2.addTranscript("trans2", 50, 160, 50, 150, 2, "trans2", "trans2", "coding");
		tx2.addExon(50, 60);
		tx2.addExon(150, 160);
		geneOverlapDetector.addLhs(gene2, gene2);

		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, ReadGeneFunctionTagger.GENE_NAME_COMPARATOR);

		r=tagger.tag(r);
		String gn = r.getStringAttribute("gn");
		String gs = r.getStringAttribute("gs");
		String gf = r.getStringAttribute("gf");

		// names always come out alphabetically sorted.
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName()+","+gene2.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+,-");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.UTR.name() + "," + LocusFunction.INTRONIC.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.UTR.name());
	}

	//One CODING read (correct strand) that's intronic (wrong strand) on another gene
	@Test
	public void testCodingCorrectIntronicWrong () {
		// read on the negative strand
		SAMRecord r = getFakeRecord(testBAMFile, 91, 100, false);

		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 100, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);

		// gene with 2 exons, 1 coding from 50-60, 1 coding from 150-160. Negative strand gene.
		GeneFromGTF gene2 = new GeneFromGTF(r.getContig(), 50, 160, true, "B", "coding", "B", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx2 = gene2.addTranscript("trans2", 50, 160, 50, 150, 2, "trans2", "trans2", "coding");
		tx2.addExon(50, 60);
		tx2.addExon(150, 160);
		geneOverlapDetector.addLhs(gene2, gene2);

		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, ReadGeneFunctionTagger.GENE_NAME_COMPARATOR);

		r=tagger.tag(r);
		String gn = r.getStringAttribute("gn");
		String gs = r.getStringAttribute("gs");
		String gf = r.getStringAttribute("gf");

		// names always come out alphabetically sorted.
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName()+","+gene2.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+,-");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.CODING.name() + "," + LocusFunction.INTRONIC.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.CODING.name());
	}

	// One intergenic read, with an intronic read on the wrong strand
	@Test
	public void testIntergenicCorrectIntronicWrong () {
		SAMRecord r = getFakeRecord(testBAMFile, 50, 60, false);
		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, true, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);
		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		r=tagger.tag(r);

		String gn = r.getStringAttribute("gn");
		String gs = r.getStringAttribute("gs");
		String gf = r.getStringAttribute("gf");

		Assert.assertEquals(gn, gene.getName());
		Assert.assertEquals(gs, "-");
		Assert.assertEquals(gf, LocusFunction.INTRONIC.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.INTERGENIC.name());
	}


	@Test
	// Read overlaps exon on one gene, and intron on another, all on the same strand.
	public void testCodingAndIntronic () {
		// read on the negative strand
		SAMRecord r = getFakeRecord(testBAMFile, 91, 100, false);

		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 100, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);

		// gene with 2 exons, 1 coding from 50-60, 1 coding from 150-160. Negative strand gene.
		GeneFromGTF gene2 = new GeneFromGTF(r.getContig(), 50, 160, false, "B", "coding", "B", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx2 = gene2.addTranscript("trans2", 50, 160, 50, 150, 2, "trans2", "trans2", "coding");
		tx2.addExon(50, 60);
		tx2.addExon(150, 160);
		geneOverlapDetector.addLhs(gene2, gene2);

		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, ReadGeneFunctionTagger.GENE_NAME_COMPARATOR);

		r=tagger.tag(r);

		// names always come out alphabetically sorted.
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName()+","+gene2.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+,+");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.CODING.name() + "," + LocusFunction.INTRONIC.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.CODING.name());
	}

	@Test(enabled=false)
	public void testAmbiguousRead() {
		SAMRecord r = getFakeRecord(testBAMFile, 50, 100, false);

		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 100, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);

		// gene with 2 exons, 1 coding from 50-60, 1 coding from 150-160. Positive strand gene.
		GeneFromGTF gene2 = new GeneFromGTF(r.getContig(), 50, 160, false, "B", "coding", "B", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx2 = gene2.addTranscript("trans2", 50, 160, 50, 150, 2, "trans2", "trans2", "coding");
		tx2.addExon(50, 60);
		tx2.addExon(150, 160);
		geneOverlapDetector.addLhs(gene2, gene2);

		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, ReadGeneFunctionTagger.GENE_NAME_COMPARATOR);

		r=tagger.tag(r);

		// names always come out alphabetically sorted.
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName()+","+gene2.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+,+");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.CODING.name() + "," + LocusFunction.INTRONIC.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.CODING.name());
	}

	@Test(enabled=true)
	// a read with 2 alignment blocks that both overlap a single gene, once in an intron and once an exon.
	public void testSplitReadExonicIntronicSameGene () {
		SAMRecord r = getFakeSplitRecord(testBAMFile, 91, 100, 131,140, false);

		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 200, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 200, 1, 200, 3, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		tx.addExon(150, 160);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);

		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, ReadGeneFunctionTagger.GENE_NAME_COMPARATOR);

		r=tagger.tag(r);
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName()+","+gene.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+,+");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.CODING.name() + "," + LocusFunction.INTRONIC.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.CODING.name());
	}

	@Test(enabled=true)
	// a read with 2 alignment blocks that both overlap a single gene, both in an exon.
	public void testSplitReadExonExonSameGene () {
		SAMRecord r = getFakeSplitRecord(testBAMFile, 91, 100, 151,160, false);

		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 200, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 200, 1, 200, 3, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		tx.addExon(150, 160);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);

		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, ReadGeneFunctionTagger.GENE_NAME_COMPARATOR);

		r=tagger.tag(r);
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.CODING.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.CODING.name());
	}

	// split read that overlaps 2 genes.  The alignment blocks both overlap gene A, but only one block overlaps gene B.
	// So, only A is a consistent gene.  B is removed.
	@Test
	public void testSplitReadExonicIntronicDifferentGenes() {
		SAMRecord r = getFakeSplitRecord(testBAMFile, 1, 10, 91,100, false);

		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 100, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);

		// gene with 2 exons, 1 coding from 50-60, 1 coding from 150-160. Positive strand gene.
		GeneFromGTF gene2 = new GeneFromGTF(r.getContig(), 50, 160, false, "B", "coding", "B", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx2 = gene2.addTranscript("trans2", 50, 160, 50, 150, 2, "trans2", "trans2", "coding");
		tx2.addExon(50, 60);
		tx2.addExon(150, 160);
		geneOverlapDetector.addLhs(gene2, gene2);

		ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, ReadGeneFunctionTagger.GENE_NAME_COMPARATOR);

		r=tagger.tag(r);
		Assert.assertEquals(r.getStringAttribute("gn"), gene.getName());
		Assert.assertEquals(r.getStringAttribute("gs"), "+");
		Assert.assertEquals(r.getStringAttribute("gf"), LocusFunction.CODING.name());
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.CODING.name());

	}

	/**
	 * Write a BAM of the test reads repeated many times, with some unmapped reads mixed in.  Big enough for many batches.
	 * Returns the number of reads written.
	 */
	private SAMRecord getFakeRecord (final File testBAMFile, final int start, final int end, final boolean negativeStrand) {
		SamReader inputSam = SamReaderFactory.makeDefault().open(testBAMFile);
		SAMRecord r = inputSam.iterator().next();
		r.setAlignmentStart(start);
		int length = (end -start) +1;
		r.setCigarString(length+"M");
		r.setMappingQuality(255);
		r.setReadNegativeStrandFlag(negativeStrand);
		r.setReadBases(Arrays.copyOf(r.getReadBases(), length));
		r.setBaseQualities(Arrays.copyOf(r.getBaseQualities(), length));
		return (r);
	}

	private SAMRecord getFakeSplitRecord (final File testBAMFile, final int start1, final int end1, final int start2, final int end2, final boolean negativeStrand) {
		SamReader inputSam = SamReaderFactory.makeDefault().open(testBAMFile);
		SAMRecord r = inputSam.iterator().next();
		r.setAlignmentStart(start1);
		int length = (end1 -start1) +1;
		int gap = ((start2-1)-(end1+1)) +1;
		int length2 = (end2 -start2) +1;
		r.setCigarString(length+"M"+gap+"N"+length2+"M");
		r.setMappingQuality(255);
		r.setReadNegativeStrandFlag(negativeStrand);
		r.setReadBases(Arrays.copyOf(r.getReadBases(), length+length2));
		r.setBaseQualities(Arrays.copyOf(r.getBaseQualities(), length+length2));
		return (r);
	}




	private SAMRecord getRecord (final File bamFile, final String recName) {
		SamReader inputSam = SamReaderFactory.makeDefault().open(testBAMFile);
		for (SAMRecord r: inputSam)
			if (r.getReadName().equals(recName))
				return (r);
		return null;
	}
}
