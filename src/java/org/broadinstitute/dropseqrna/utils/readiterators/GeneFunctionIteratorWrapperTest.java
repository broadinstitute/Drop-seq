package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import junit.framework.Assert;
import org.testng.annotations.Test;
import picard.annotation.LocusFunction;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

public class GeneFunctionIteratorWrapperTest {

	File testBAMFile= new File ("testdata/org/broadinstitute/dropseq/private/metrics/NucBYReg4Reg.MOUSE.GCTAAGTAAGAT.Elp2.fixed.bam");


	@Test
	// don't need to test FunctionalDataProcessorTest functions, but do need to test assignReadsToAllGenes and decoding LocusFunction from string attributes.
	public void testAssignReadsToAllGenesTrue() {
		String geneTag="gn";
		String strandTag="gs";
		String functionTag="gf";
		boolean assignReadsToAllGenes=true;
		StrandStrategy strandFilterStrategy = StrandStrategy.SENSE;
		LocusFunction [] alf = {LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC};
		Collection<LocusFunction> acceptedLociFunctions = Arrays.asList(alf);

		// set up test data.  for assignReadsToAllGenes=T, this returns two SAMRecords.
		String recName="H575CBGXY:4:23606:1714:16102";
		SAMRecord r1 = getRecord(testBAMFile, recName);
		r1.setAttribute(geneTag, "A,B");
		r1.setAttribute(strandTag, "+,+");
		r1.setAttribute(functionTag, "CODING,CODING");

		List<SAMRecord> list = new ArrayList<SAMRecord>();
		list.add(r1);

		GeneFunctionIteratorWrapper gfiw = new GeneFunctionIteratorWrapper(list.iterator(), geneTag, strandTag, functionTag, assignReadsToAllGenes, strandFilterStrategy, acceptedLociFunctions);
		int counter=0;
		while (gfiw.hasNext()) {
			SAMRecord result = gfiw.next();
			if (result.getStringAttribute(geneTag).equals("A")) {
				Assert.assertEquals(result.getStringAttribute(geneTag), "A");
				Assert.assertEquals(result.getStringAttribute(strandTag), "+");
				Assert.assertEquals(result.getStringAttribute(functionTag), "CODING");
			}
			if (result.getStringAttribute(geneTag).equals("B")) {
				Assert.assertEquals(result.getStringAttribute(geneTag), "B");
				Assert.assertEquals(result.getStringAttribute(strandTag), "+");
				Assert.assertEquals(result.getStringAttribute(functionTag), "CODING");
			}
			counter++;
		}

		Assert.assertEquals(2, counter);
	}

	@Test
	// don't need to test FunctionalDataProcessorTest functions, but do need to test assignReadsToAllGenes and decoding LocusFunction from string attributes.
	public void testAssignReadsToAllGenesFalse() {
		String geneTag="gn";
		String strandTag="gs";
		String functionTag="gf";
		boolean assignReadsToAllGenes=false;
		StrandStrategy strandFilterStrategy = StrandStrategy.SENSE;
		LocusFunction [] alf = {LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC};
		Collection<LocusFunction> acceptedLociFunctions = Arrays.asList(alf);

		// set up test data.
		String recName="H575CBGXY:4:23606:1714:16102";
		SAMRecord r1 = getRecord(testBAMFile, recName);
		r1.setAttribute(geneTag, "A,B,C");
		r1.setAttribute(strandTag, "+,+,+");
		r1.setAttribute(functionTag, "CODING,INTRONIC,CODING");

		List<SAMRecord> list = new ArrayList<SAMRecord>();
		list.add(r1);

		GeneFunctionIteratorWrapper gfiw = new GeneFunctionIteratorWrapper(list.iterator(), geneTag, strandTag, functionTag, assignReadsToAllGenes, strandFilterStrategy, acceptedLociFunctions);
		int counter=0;
		while (gfiw.hasNext()) {
			SAMRecord result = gfiw.next();
			if (result.getStringAttribute(geneTag).equals("A")) {
				Assert.assertEquals(result.getStringAttribute(geneTag), "A");
				Assert.assertEquals(result.getStringAttribute(strandTag), "+");
				Assert.assertEquals(result.getStringAttribute(functionTag), "CODING");
			}
			if (result.getStringAttribute(geneTag).equals("C")) {
				Assert.assertEquals(result.getStringAttribute(geneTag), "B");
				Assert.assertEquals(result.getStringAttribute(strandTag), "+");
				Assert.assertEquals(result.getStringAttribute(functionTag), "CODING");
			}
			counter++;
		}
		// No reads emerge.
		Assert.assertEquals(0, counter);
	}

	@Test (enabled=true)
	// don't need to test FunctionalDataProcessorTest functions, but do need to test assignReadsToAllGenes and decoding LocusFunction from string attributes.
	public void testAmbiguousGenes() {
		String geneTag="gn";
		String strandTag="gs";
		String functionTag="gf";

		StrandStrategy strandFilterStrategy = StrandStrategy.SENSE;
		LocusFunction [] alf = {LocusFunction.CODING, LocusFunction.UTR};
		Collection<LocusFunction> acceptedLociFunctions = Arrays.asList(alf);

		// set up test data.  for assignReadsToAllGenes=T, this returns two SAMRecords.
		String recName="H575CBGXY:4:23606:1714:16102";
		SAMRecord r1 = getRecord(testBAMFile, recName);
		r1.setAttribute(geneTag, "A,B");
		r1.setAttribute(strandTag, "+,+");
		r1.setAttribute(functionTag, "CODING,CODING");

		List<SAMRecord> list = new ArrayList<SAMRecord>();
		list.add(r1);
		// should remove read because it's ambiguous.
		boolean assignReadsToAllGenes=false;
		GeneFunctionIteratorWrapper gfiw = new GeneFunctionIteratorWrapper(list.iterator(), geneTag, strandTag, functionTag, assignReadsToAllGenes, strandFilterStrategy, acceptedLociFunctions);
		int counter=0;
		while (gfiw.hasNext()) {

		}
		// neither of the reads will emerge because they are ambiguously assiged.
		Assert.assertEquals(0, counter);

		// should return 2 reads because it's ambiguous and we're ok with that.
		assignReadsToAllGenes=true;
		gfiw = new GeneFunctionIteratorWrapper(list.iterator(), geneTag, strandTag, functionTag, assignReadsToAllGenes, strandFilterStrategy, acceptedLociFunctions);
		counter=0;
		int countA=0;
		int countB=0;
		while (gfiw.hasNext()) {
			SAMRecord result = gfiw.next();
			if (result.getStringAttribute(geneTag).equals("A")) {
				countA++;
				Assert.assertEquals(result.getStringAttribute(geneTag), "A");
				Assert.assertEquals(result.getStringAttribute(strandTag), "+");
				Assert.assertEquals(result.getStringAttribute(functionTag), "CODING");
			}
			if (result.getStringAttribute(geneTag).equals("B")) {
				countB++;
				Assert.assertEquals(result.getStringAttribute(geneTag), "B");
				Assert.assertEquals(result.getStringAttribute(strandTag), "+");
				Assert.assertEquals(result.getStringAttribute(functionTag), "CODING");
			}
			counter++;
		}
		// Both reads will emerge.
		Assert.assertEquals(2, counter);
		Assert.assertEquals(1, countA);
		Assert.assertEquals(1, countB);
	}

	private SAMRecord getRecord (final File bamFile, final String recName) {
		SamReader inputSam = SamReaderFactory.makeDefault().open(testBAMFile);
		for (SAMRecord r: inputSam)
			if (r.getReadName().equals(recName))
				return (r);
		return null;
	}
}
