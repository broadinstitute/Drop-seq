package org.broadinstitute.dropseqrna.annotation;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.OverlapDetector;
import picard.annotation.Gene;
import picard.annotation.LocusFunction;

public class AnnotationUtilsTest {

	File GTF = new File ("testdata/org/broadinstitute/dropseq/annotation/test.gtf.gz");
	File BAM = new File ("testdata/org/broadinstitute/dropseq/annotation/test.bam");

	private final static double PCT_REQUIRED_OVERLAP_DROPSEQ=0;
	private final static double PCT_REQUIRED_OVERLAP_STARSOLO=50;
	@Test
	public void testGetFunctionalDataForRead() {
		AnnotationUtils u = AnnotationUtils.getInstance();
		OverlapDetector<Gene> geneOverlapDetector = getGeneOverlapDetector(BAM, GTF);
		Map<String, Gene> map = getGeneMap(geneOverlapDetector);
		// expected UTR and coding.
		SAMRecord rec = getReadByName("000000000-AMY9M:1:2104:20987:12124", BAM);
		Map<Gene, List<LocusFunction>> funcMap = u.getFunctionalDataForRead(rec, geneOverlapDetector, PCT_REQUIRED_OVERLAP_DROPSEQ);
		Set<LocusFunction> actual = new HashSet<>(funcMap.get(map.get("RPL22")));
		Set<LocusFunction> expected = new HashSet<> (Arrays.asList(LocusFunction.CODING, LocusFunction.UTR));
		Assert.assertEquals(actual, expected);

		// expected CODING once for each gene.
		rec = getReadByName("000000000-AMY9M:1:1114:17300:12082", BAM);
		funcMap = u.getFunctionalDataForRead(rec, geneOverlapDetector, PCT_REQUIRED_OVERLAP_DROPSEQ);
		actual = new HashSet<>(funcMap.get(map.get("PLEKHM2")));
		expected = new HashSet<> (Arrays.asList(LocusFunction.CODING));
		Assert.assertEquals(actual, expected);
		actual = new HashSet<>(funcMap.get(map.get("RP11-288I21.1")));
		Assert.assertEquals(actual, expected);

		// expected CODING AND INTRONIC for the same gene.
		rec = getReadByName("000000000-AMY9M:1:2109:19958:23569", BAM);
		funcMap = u.getFunctionalDataForRead(rec, geneOverlapDetector, PCT_REQUIRED_OVERLAP_DROPSEQ);
		actual = new HashSet<>(funcMap.get(map.get("UTP11L")));
		expected = new HashSet<> (Arrays.asList(LocusFunction.CODING, LocusFunction.INTRONIC));
		Assert.assertEquals(actual, expected);

		//  CODING, CODING, INTRONIC.
		rec = getReadByName("000000000-AMY9M:1:2104:15433:7166", BAM);
		funcMap = u.getFunctionalDataForRead(rec, geneOverlapDetector, PCT_REQUIRED_OVERLAP_DROPSEQ);
		actual = new HashSet<>(funcMap.get(map.get("RP11-288I21.1")));
		expected = new HashSet<> (Arrays.asList(LocusFunction.CODING, LocusFunction.INTRONIC));
		Assert.assertEquals(actual, expected);



	}

	/**
	 * Test functional annotations where >50% of the read must match the locus function.
	 */
	@Test
	public void testGetFunctionalDataForReadStarSolo() {
		AnnotationUtils u = AnnotationUtils.getInstance();
		OverlapDetector<Gene> geneOverlapDetector = getGeneOverlapDetector(BAM, GTF);
		Map<String, Gene> map = getGeneMap(geneOverlapDetector);
		// expected UTR and coding.
		SAMRecord rec = getReadByName("000000000-AMY9M:1:2104:20987:12124", BAM);
		Map<Gene, List<LocusFunction>> funcMap = u.getFunctionalDataForRead(rec, geneOverlapDetector, PCT_REQUIRED_OVERLAP_STARSOLO);
		Set<LocusFunction> actual = new HashSet<>(funcMap.get(map.get("RPL22")));
		// UTR has the majority of bases.
		Set<LocusFunction> expected = new HashSet<> (Arrays.asList(LocusFunction.UTR, LocusFunction.CODING));
		Assert.assertEquals(actual, expected);

		// expected CODING once for each gene.
		rec = getReadByName("000000000-AMY9M:1:1114:17300:12082", BAM);
		funcMap = u.getFunctionalDataForRead(rec, geneOverlapDetector, PCT_REQUIRED_OVERLAP_STARSOLO);
		actual = new HashSet<>(funcMap.get(map.get("PLEKHM2")));
		expected = new HashSet<> (Arrays.asList(LocusFunction.CODING));
		Assert.assertEquals(actual, expected);
		actual = new HashSet<>(funcMap.get(map.get("RP11-288I21.1")));
		Assert.assertEquals(actual, expected);

		// expected CODING AND INTRONIC for the same gene.
		rec = getReadByName("000000000-AMY9M:1:2109:19958:23569", BAM);
		funcMap = u.getFunctionalDataForRead(rec, geneOverlapDetector, PCT_REQUIRED_OVERLAP_STARSOLO);
		actual = new HashSet<>(funcMap.get(map.get("UTP11L")));
		// CODING 38, INTRONIC 2 - coding is the annotation > 50%
		expected = new HashSet<> (Arrays.asList(LocusFunction.CODING));
		Assert.assertEquals(actual, expected);

		//  CODING, INTRONIC.
		rec = getReadByName("000000000-AMY9M:1:2104:15433:7166", BAM);
		funcMap = u.getFunctionalDataForRead(rec, geneOverlapDetector, PCT_REQUIRED_OVERLAP_STARSOLO);
		// {CODING=6, INTRONIC=44}
		actual = new HashSet<>(funcMap.get(map.get("RP11-288I21.1")));
		expected = new HashSet<> (Arrays.asList(LocusFunction.INTRONIC));
		Assert.assertEquals(actual, expected);
	}
	@Test
	public void testGetLocusFunctionForReadByGene () {
		AnnotationUtils u = AnnotationUtils.getInstance();
		OverlapDetector<Gene> geneOverlapDetector = getGeneOverlapDetector(BAM, GTF);
		Map<String, Gene> map = getGeneMap(geneOverlapDetector);
		// expected UTR and coding.
		SAMRecord rec = getReadByName("000000000-AMY9M:1:2104:20987:12124", BAM);
		Map<Gene, LocusFunction> funcMap = u.getLocusFunctionForReadByGene(rec, geneOverlapDetector, PCT_REQUIRED_OVERLAP_DROPSEQ);
		LocusFunction actual = funcMap.get(map.get("RPL22"));
		LocusFunction expected = LocusFunction.CODING;
		Assert.assertEquals(actual, expected);

		// expected CODING,INTRONIC->CODING.
		rec = getReadByName("000000000-AMY9M:1:2109:19958:23569", BAM);
		funcMap = u.getLocusFunctionForReadByGene(rec, geneOverlapDetector, PCT_REQUIRED_OVERLAP_DROPSEQ);
		actual = funcMap.get(map.get("UTP11L"));
		expected = LocusFunction.CODING;
		Assert.assertEquals(actual, expected);

		//
		rec = getReadByName("000000000-AMY9M:1:2104:15433:7166", BAM);
		funcMap = u.getLocusFunctionForReadByGene(rec, geneOverlapDetector, PCT_REQUIRED_OVERLAP_DROPSEQ);
		actual = funcMap.get(map.get("PLEKHM2"));
		expected = LocusFunction.CODING;
		Assert.assertEquals(actual, expected);
		// RP11-288I21.1 is both CODING and INTRONIC, CODING wins.
		actual = funcMap.get(map.get("RP11-288I21.1"));
		expected = LocusFunction.CODING;
		Assert.assertEquals(actual, expected);

	}

	@Test
	public void testGetConsistentExons () {
		// need to build a new read that falls into 2 genes to check if it's inconsistent.
		AnnotationUtils u = AnnotationUtils.getInstance();
		OverlapDetector<Gene> geneOverlapDetector = getGeneOverlapDetector(BAM, GTF);
		Map<String, Gene> map = getGeneMap(geneOverlapDetector);
		// expected UTR and coding.
		SAMRecord rec = getReadByName("000000000-AMY9M:1:2104:20987:12124", BAM);
		// change the cigar so the read maps to two genes exons.
		// rec.setCigar(cigar);
		rec.setAlignmentStart(6269430);
		Cigar c = new Cigar ();
		c.add(new CigarElement(20, CigarOperator.M));
		c.add(new CigarElement(9741377, CigarOperator.N));
		c.add(new CigarElement(30, CigarOperator.M));
		rec.setCigar(c);
		Set<Gene> genes = new HashSet<>();
		genes.add(map.get("RPL22"));
		genes.add(map.get("PLEKHM2"));

		// This read spans 2 genes, so it's not consistent under the strict test.
		Set<Gene> actual = u.getConsistentExons(rec, genes, false);
		Set<Gene> expected = new HashSet<>();
		Assert.assertEquals(actual, expected);

		// under the less strict option, we allow the read to span 2 different genes.
		actual = u.getConsistentExons(rec, genes, true);
		expected = genes;
		Assert.assertEquals(actual, expected);

	}

	@Test
	public void testGetLocusFunction () {
		AnnotationUtils u = AnnotationUtils.getInstance();

		// Strict means all functions must be the same.
		Set <LocusFunction> c = new HashSet<>(Arrays.asList(LocusFunction.CODING, LocusFunction.CODING));
		LocusFunction result = u.getLocusFunction(c, true);
		Assert.assertEquals(result, LocusFunction.CODING);

		// functions don't match, conservative = null
		c = new HashSet<>(Arrays.asList(LocusFunction.CODING, LocusFunction.INTRONIC));
		result = u.getLocusFunction(c, true);
		Assert.assertNull (result);

		// non conservative collapses to the "best" [closest to coding] locus function.
		result = u.getLocusFunction(c, false);
		Assert.assertEquals(result, LocusFunction.CODING);

		c = new HashSet<>(Arrays.asList(LocusFunction.RIBOSOMAL, LocusFunction.INTRONIC));
		result = u.getLocusFunction(c, true);
		Assert.assertNull (result);
		result = u.getLocusFunction(c, false);
		Assert.assertEquals(result, LocusFunction.INTRONIC);

	}

	@Test
	public void testGetLocusFunctionsByBlock() {
		AnnotationUtils u = AnnotationUtils.getInstance();
		OverlapDetector<Gene> geneOverlapDetector = getGeneOverlapDetector(BAM, GTF);
		Map<String, Gene> map = getGeneMap(geneOverlapDetector);

		// this interval is either UTR or CODING at every base.
		Interval i = new Interval ("1", 16010827 , 16011113);
		// Set<Gene> overlappingGenes = new HashSet<>(Arrays.asList(map.get("PLEKHM2")));
		Set<Gene> overlappingGenes = geneOverlapDetector.getOverlaps(i);
		LocusFunction [] result = u.getLocusFunctionsByBlock(i, overlappingGenes);
		for (LocusFunction r: result) {
			boolean test = r==LocusFunction.CODING || r==LocusFunction.UTR;
			Assert.assertTrue(test);
		}

		// this interval is intergenic
		i = new Interval ("1", 16000000 , 16000050);
		overlappingGenes = geneOverlapDetector.getOverlaps(i);
		result = u.getLocusFunctionsByBlock(i, overlappingGenes);
		for (LocusFunction r: result) {
			boolean test = r==LocusFunction.INTERGENIC;
			Assert.assertTrue(test);
		}

		// this interval is intronic.
		i = new Interval ("1", 16020000 , 16020050);
		overlappingGenes = geneOverlapDetector.getOverlaps(i);
		result = u.getLocusFunctionsByBlock(i, overlappingGenes);
		for (LocusFunction r: result) {
			boolean test = r==LocusFunction.INTRONIC;
			Assert.assertTrue(test);
		}



	}

	@Test
	// A block that is a mix of INTERGENIC, INTRONIC, UTR and CODING bases must keep the per-base mix.
	public void testGetLocusFunctionsByIntervalMixed() {
		AnnotationUtils u = AnnotationUtils.getInstance();
		GeneFromGTF gene = new GeneFromGTF("1", 11, 100, false, "A", "coding", "A", "coding", 1);
		// exons 11-20 (UTR 11-14, coding 15-20), 41-50 (coding), 81-100 (coding 81-90, UTR 91-100)
		final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("trans1", 11, 100, 15, 90, 3, "trans1", "trans1", "coding");
		tx.addExon(11, 20);
		tx.addExon(41, 50);
		tx.addExon(81, 100);
		LocusFunction[] actual = u.getLocusFunctionsByInterval(new Interval("1", 8, 105), gene);
		Assert.assertEquals(actual.length, 98);
		for (int pos=8; pos<=105; pos++) {
			LocusFunction expected;
			if (pos<11 || pos>100) expected=LocusFunction.INTERGENIC;
			else if (pos<=14) expected=LocusFunction.UTR;
			else if (pos<=20) expected=LocusFunction.CODING;
			else if (pos<41) expected=LocusFunction.INTRONIC;
			else if (pos<=50) expected=LocusFunction.CODING;
			else if (pos<81) expected=LocusFunction.INTRONIC;
			else if (pos<=90) expected=LocusFunction.CODING;
			else expected=LocusFunction.UTR;
			Assert.assertEquals(actual[pos-8], expected, "position " + pos);
		}
	}

	@Test
	// The range based implementation must give the same per-base answer as picard's Transcript.assignLocusFunctionForRange.
	public void testGetLocusFunctionsByIntervalMatchesPicard() {
		AnnotationUtils u = AnnotationUtils.getInstance();
		Random random = new Random(42);
		for (int g=0; g<300; g++) {
			GeneFromGTF gene = new GeneFromGTF("1", 100, 400, false, "G"+g, "coding", "G"+g, "coding", 1);
			int numTranscripts = 1 + random.nextInt(3);
			for (int t=0; t<numTranscripts; t++) {
				int txStart = 100 + random.nextInt(100);
				int txEnd = 300 + random.nextInt(100);
				// exons sorted by start, may overlap each other.
				int numExons = 1 + random.nextInt(5);
				int[] starts = new int[numExons];
				for (int e=0; e<numExons; e++) starts[e]=txStart + random.nextInt(txEnd-txStart);
				Arrays.sort(starts);
				// coding region is sometimes degenerate (non-coding transcript).
				int codingStart = txStart + random.nextInt(txEnd-txStart);
				int codingEnd = random.nextInt(4)==0 ? codingStart-1 : Math.min(txEnd, codingStart + random.nextInt(150));
				final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript("t"+t, txStart, txEnd, codingStart, codingEnd, numExons, "t"+t, "t"+t, "coding");
				for (int e=0; e<numExons; e++)
					tx.addExon(starts[e], Math.min(txEnd, starts[e] + random.nextInt(60)));
			}
			for (int k=0; k<30; k++) {
				int start = 50 + random.nextInt(400);
				int end = start + random.nextInt(120);
				LocusFunction[] expected = new LocusFunction[end-start+1];
				Arrays.fill(expected, LocusFunction.INTERGENIC);
				for (Gene.Transcript t: gene)
					t.assignLocusFunctionForRange(start, expected);
				LocusFunction[] actual = u.getLocusFunctionsByInterval(new Interval("1", start, end), gene);
				Assert.assertEquals(actual, expected, "gene " + g + " interval " + start + "-" + end);
			}
		}
	}

	/**
	 * Only suitable for small BAMs.
	 * @param name
	 * @param bamFile
	 * @return
	 */
	private SAMRecord getReadByName (final String name, final File bamFile) {
		SamReader inputSam = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(bamFile);
		Iterator<SAMRecord> iter = inputSam.iterator();
		boolean found = false;
		SAMRecord result=null;
		while (iter.hasNext() && !found) {
			SAMRecord r = iter.next();
			if (r.getReadName().equals(name)) {
				found=true;
				result=r;
			}
		}
		CloserUtil.close(inputSam);
		return result;
	}

	private OverlapDetector<Gene> getGeneOverlapDetector (final File bam, final File gtf) {
		SamReader inputSam = SamReaderFactory.makeDefault().open(bam);
		SAMFileHeader header = inputSam.getFileHeader();
		SAMSequenceDictionary bamDict = header.getSequenceDictionary();
		return GeneAnnotationReader.loadAnnotationsFile(gtf, bamDict);
	}

	private Map<String, Gene> getGeneMap (final OverlapDetector<Gene> genes) {
		Map<String, Gene> result = new HashMap<>();
		for (Gene g: genes.getAll())
			result.put(g.getName(), g);
		return (result);
	}


}
