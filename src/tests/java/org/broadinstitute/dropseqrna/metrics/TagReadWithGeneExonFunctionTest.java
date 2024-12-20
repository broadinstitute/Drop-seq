package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.dropseqrna.annotation.GeneAnnotationReader;
import org.broadinstitute.dropseqrna.annotation.GeneFromGTF;
import org.broadinstitute.dropseqrna.utils.CompareBAMTagValues;
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
import picard.nio.PicardHtsPath;

public class TagReadWithGeneExonFunctionTest {

	File testBAMFile= new File ("testdata/org/broadinstitute/dropseq/metrics/NucBYReg4Reg.MOUSE.GCTAAGTAAGAT.Elp2.fixed.bam");
	File annotationsFile=new File ("testdata/org/broadinstitute/dropseq/metrics/mm10_Elp2.gtf");


	@Test
	public void testDoWork() {
		TagReadWithGeneExonFunction t = new TagReadWithGeneExonFunction();
		File tempBAM = getTempReportFile("TagReadWithGeneExonFunctionTest", ".bam");
		File tempSummary=getTempReportFile("TagReadWithGeneExonFunctionTest", ".summary");
		tempBAM.deleteOnExit();
		tempSummary.deleteOnExit();

		t.INPUT=testBAMFile;
		t.OUTPUT=tempBAM;
		t.ANNOTATIONS_FILE=annotationsFile;
		t.SUMMARY=tempSummary;
		int returnVal = t.doWork();
		Assert.assertTrue(returnVal==0);

		// test output BAM
		List<String> tags = new ArrayList<>(Arrays.asList("XC", "GE", "GS", "XF"));
		compareBAMTagValues(testBAMFile, tempBAM, tags, 0);

	}

	@Test
	public void testTagIntronRead () {
		SamReader inputSam = SamReaderFactory.makeDefault().open(testBAMFile);
		SAMFileHeader header = inputSam.getFileHeader();
		SAMSequenceDictionary bamDict = header.getSequenceDictionary();
		final OverlapDetector<Gene> geneOverlapDetector = GeneAnnotationReader.loadAnnotationsFile(annotationsFile, bamDict);

		String recName="H575CBGXY:4:23606:1714:16102";
		SAMRecord r = getRecord(testBAMFile, recName);

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		r=tagger.setAnnotations(r, geneOverlapDetector);

		String geneName = r.getStringAttribute("GE");
		String geneStrand = r.getStringAttribute("GS");
		String function = r.getStringAttribute("XF");

		Assert.assertNull(geneName);
		Assert.assertNull(geneStrand);
		Assert.assertEquals(function, "INTRONIC");

	}

	@Test
	public void testTagCodingRead () {
		SamReader inputSam = SamReaderFactory.makeDefault().open(testBAMFile);
		SAMFileHeader header = inputSam.getFileHeader();
		SAMSequenceDictionary bamDict = header.getSequenceDictionary();
		final OverlapDetector<Gene> geneOverlapDetector = GeneAnnotationReader.loadAnnotationsFile(annotationsFile, bamDict);

		String recName="H575CBGXY:3:13502:3588:9959";
		SAMRecord r = getRecord(testBAMFile, recName);

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		r=tagger.setAnnotations(r, geneOverlapDetector);

		Assert.assertEquals(r.getStringAttribute("GE"), "Elp2");
		Assert.assertEquals(r.getStringAttribute("GS"), "+");
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

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		r=tagger.setAnnotations(r, geneOverlapDetector);

		String geneName = r.getStringAttribute("GE");
		String geneStrand = r.getStringAttribute("GS");

		Assert.assertEquals(geneName, "Elp2");
		Assert.assertEquals(geneStrand, "+");
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
		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		r=tagger.setAnnotations(r, geneOverlapDetector);

		Assert.assertEquals(r.getStringAttribute("GE"), null);
		Assert.assertEquals(r.getStringAttribute("GS"), null);
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
		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		r=tagger.setAnnotations(r, geneOverlapDetector);

		String geneName = r.getStringAttribute("GE");
		String geneStrand = r.getStringAttribute("GS");
		String XF = r.getStringAttribute("XF");

		Assert.assertNull(geneName);
		Assert.assertNull(geneStrand);

		Assert.assertEquals(XF, LocusFunction.INTRONIC.name());
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
		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		r=tagger.setAnnotations(r, geneOverlapDetector);
		Assert.assertEquals(r.getStringAttribute("GE"), gene.getName());
		Assert.assertEquals(r.getStringAttribute("GS"), "+");
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
		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		r=tagger.setAnnotations(r, geneOverlapDetector);

		String GE=r.getStringAttribute("GE");
		String GS=r.getStringAttribute("GS");
		String XF = r.getStringAttribute("XF");

		Assert.assertNull(GE);
		Assert.assertNull(GS);
		Assert.assertEquals(XF, LocusFunction.INTERGENIC.name());
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
		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		r=tagger.setAnnotations(r, geneOverlapDetector);


		Assert.assertEquals(r.getStringAttribute("GE"), gene.getName());
		Assert.assertEquals(r.getStringAttribute("GS"), "+");
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.CODING.name());
	}

	@Test
	// One CODING read (correct strand, 2nd gene has overlapping coordinates with read)==ambiguous.
	public void testCodingReadAmbiguous () {
		SAMRecord r = getFakeRecord(testBAMFile, 2, 8, false);
		// gene with 2 exons, 1 coding from 1-10, one UTR from 91-100.  Positive strand gene.
		GeneFromGTF gene = new GeneFromGTF(r.getContig(), 1, 100, false, "A", "coding", "A", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx.addExon(1, 10);
		tx.addExon(91, 100);
		GeneFromGTF gene2 = new GeneFromGTF(r.getContig(), 1, 100, false, "B", "coding", "B", "coding", 1);
		final GeneFromGTF.TranscriptFromGTF tx2 = gene2.addTranscript("trans1", 1, 100, 1, 90, 2, "trans1", "trans1", "coding");
		tx2.addExon(1, 10);
		tx2.addExon(91, 100);
		OverlapDetector<Gene> geneOverlapDetector = new OverlapDetector<>(0, 0);
		geneOverlapDetector.addLhs(gene, gene);
		geneOverlapDetector.addLhs(gene2, gene2);
		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		r=tagger.setAnnotations(r, geneOverlapDetector);


		Assert.assertNull(r.getStringAttribute("GE"));
		Assert.assertNull(r.getStringAttribute("GS"));
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
		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		r=tagger.setAnnotations(r, geneOverlapDetector);

		String geneTagged=r.getStringAttribute("GE");
		String strandTagged=r.getStringAttribute("GS");

		Assert.assertNull(geneTagged);
		Assert.assertNull(strandTagged);
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

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, TagReadWithGeneFunction.GENE_NAME_COMPARATOR);

		r=tagger.setAnnotations(r, geneOverlapDetector);
		String GE = r.getStringAttribute("GE");
		String GS = r.getStringAttribute("GS");
		String XF = r.getStringAttribute("XF");

		// names always come out alphabetically sorted.
		Assert.assertNull(GE);
		Assert.assertNull(GS);

		Assert.assertEquals(XF, LocusFunction.INTRONIC.name());
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

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, TagReadWithGeneFunction.GENE_NAME_COMPARATOR);

		r=tagger.setAnnotations(r, geneOverlapDetector);
		String ge = r.getStringAttribute("GE");
		String gs = r.getStringAttribute("GS");
		String xf = r.getStringAttribute("XF");

		Assert.assertNull(ge);
		Assert.assertNull(gs);

		Assert.assertEquals(xf, LocusFunction.INTRONIC.name());
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

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, TagReadWithGeneFunction.GENE_NAME_COMPARATOR);

		r=tagger.setAnnotations(r, geneOverlapDetector);
		String GE = r.getStringAttribute("GE");
		String GS = r.getStringAttribute("GS");
		String XF = r.getStringAttribute("XF");

		// names always come out alphabetically sorted.
		Assert.assertEquals(GE, gene.getName());
		Assert.assertEquals(GS, "+");
		Assert.assertEquals(XF, LocusFunction.UTR.name());
	}

	//One CODING read (correct strand) that's intronic (wrong strand) on another gene
	@Test
	public void testCodingCorrectIntronicWrong () {
		// read on the positive strand
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

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, TagReadWithGeneFunction.GENE_NAME_COMPARATOR);

		r=tagger.setAnnotations(r, geneOverlapDetector);
		String geneTag = r.getStringAttribute("GE");
		String geneStrand = r.getStringAttribute("GS");


		// names always come out alphabetically sorted.
		Assert.assertEquals(geneTag, gene.getName());
		Assert.assertEquals(geneStrand, "+");
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
		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		r=tagger.setAnnotations(r, geneOverlapDetector);

		String geneTag = r.getStringAttribute("GE");
		String geneStrand = r.getStringAttribute("GS");

		Assert.assertNull(geneTag);
		Assert.assertNull(geneStrand);
		// Assert.assertEquals(geneTag, gene.getName());
		// Assert.assertEquals(geneStrand, "-");
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.INTERGENIC.name());
	}


	@Test
	// Read overlaps exon on one gene, and intron on another, all on the same strand.
	public void testCodingAndIntronic () {
		// read on the positive strand
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

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, TagReadWithGeneFunction.GENE_NAME_COMPARATOR);

		r=tagger.setAnnotations(r, geneOverlapDetector);

		String geneTagged = r.getStringAttribute("GE");
		String strandTagged = r.getStringAttribute("GS");

		// names always come out alphabetically sorted.
		Assert.assertEquals(geneTagged, gene.getName());
		Assert.assertEquals(strandTagged, "+");
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

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, TagReadWithGeneFunction.GENE_NAME_COMPARATOR);

		r=tagger.setAnnotations(r, geneOverlapDetector);

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

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, TagReadWithGeneFunction.GENE_NAME_COMPARATOR);

		r=tagger.setAnnotations(r, geneOverlapDetector);
		Assert.assertEquals(r.getStringAttribute("GE"), gene.getName());
		Assert.assertEquals(r.getStringAttribute("GS"), "+");
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

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, TagReadWithGeneFunction.GENE_NAME_COMPARATOR);

		r=tagger.setAnnotations(r, geneOverlapDetector);
		Assert.assertEquals(r.getStringAttribute("GE"), gene.getName());
		Assert.assertEquals(r.getStringAttribute("GS"), "+");
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

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		List <Gene> genes = new ArrayList <> (geneOverlapDetector.getAll());
		Collections.sort(genes, TagReadWithGeneFunction.GENE_NAME_COMPARATOR);

		r=tagger.setAnnotations(r, geneOverlapDetector);
		Assert.assertEquals(r.getStringAttribute("GE"), gene.getName());
		Assert.assertEquals(r.getStringAttribute("GS"), "+");
		Assert.assertEquals(r.getStringAttribute("XF"), LocusFunction.CODING.name());

	}

	@Test
	public void testGetCompoundGeneName () {
		GeneFromGTF gene = new GeneFromGTF("1", 1, 100, true, "A", "coding", "A", "coding", 1);
		GeneFromGTF gene2 = new GeneFromGTF("1", 1, 100, true, "B", "coding", "B", "coding", 1);
		List<Gene> genes = new ArrayList<> (Arrays.asList(gene, gene2));

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		String result = tagger.getCompoundGeneName(genes);
		String expected="A" + tagger.RECORD_SEP+"B";
		Assert.assertEquals(result, expected);
	}

	@Test
	public void testGetCompoundStrand () {
		GeneFromGTF gene = new GeneFromGTF("1", 1, 100, true, "A", "coding", "A", "coding", 1);
		GeneFromGTF gene2 = new GeneFromGTF("1", 1, 100, false, "B", "coding", "B", "coding", 1);
		List<Gene> genes = new ArrayList<> (Arrays.asList(gene, gene2));

		TagReadWithGeneExonFunction tagger = new TagReadWithGeneExonFunction();
		String result = tagger.getCompoundStrand(genes);
		String expected="-" + tagger.RECORD_SEP+"+";
		Assert.assertEquals(result, expected);
	}


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

	private File getTempReportFile (final String prefix, final String suffix) {
		File tempFile=null;

		try {
			tempFile = File.createTempFile(prefix, suffix);
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		return tempFile;
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
