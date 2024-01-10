package org.broadinstitute.dropseqrna.utils.readiterators;

import org.broadinstitute.dropseqrna.annotation.functionaldata.DropSeqFunctionalDataProcessor;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;
import org.broadinstitute.dropseqrna.annotation.functionaldata.AbstractFunctionalDataProcessor;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorI;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.annotation.LocusFunction;

import java.util.List;

public class DropSeqFunctionalDataProcessorTest {


	@Test
	public void testSingleFD () {
		LocusFunction [] acceptedFunctions = {LocusFunction.CODING, LocusFunction.UTR};
		FunctionalDataProcessorI fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.SENSE, acceptedFunctions);
		String [] genes = {"A"};
		String [] strands = {"+"};
		LocusFunction [] locusFunctions = {LocusFunction.CODING};
		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, false);
		Assert.assertEquals(fdList.size(), 1);
		FunctionalData fd = fdList.get(0);
		Assert.assertEquals(fd.getGene(), genes[0]);
		Assert.assertEquals(fd.getStrand(), strands[0]);
		Assert.assertEquals(fd.getLocusFunction(), locusFunctions[0]);

		fdList = fdp.filterToPreferredAnnotations(fdList);
		Assert.assertEquals(fdList.size(), 1);

	}

	@Test
	public void testNoResult () {
		LocusFunction [] acceptedFunctions = {LocusFunction.CODING, LocusFunction.UTR};
		FunctionalDataProcessorI fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.SENSE, acceptedFunctions);
		String [] genes = {"A"};
		String [] strands = {"+"};
		LocusFunction [] locusFunctions = {LocusFunction.INTRONIC};
		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, false);
		Assert.assertEquals(fdList.size(), 0);
		fdList = fdp.filterToPreferredAnnotations(fdList);
		Assert.assertEquals(fdList.size(), 0);
	}

	//Desired intronic reads, have intronic and CODING on same gene.  Return nothing, as CODING wasn't specified.
	// Then make both CODING and INTRONIC accepted and return the functional data.  In this case, CODING is preferred over intronic.
	@Test
	public void testFunctionalFiltering () {
		LocusFunction [] acceptedFunctions = {LocusFunction.INTRONIC};
		FunctionalDataProcessorI fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.SENSE, acceptedFunctions);
		String [] genes = {"A", "A"};
		String [] strands = {"+", "+"};
		LocusFunction [] locusFunctions = {LocusFunction.CODING, LocusFunction.INTRONIC};
		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, false);
		Assert.assertEquals(fdList.size(), 0);

		// now accept UTR and INTRONIC.
		LocusFunction [] acceptedFunctions2 = {LocusFunction.INTRONIC, LocusFunction.CODING};
		fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.SENSE, acceptedFunctions2);
		fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, false);
		Assert.assertEquals(fdList.size(), 1);

		FunctionalData fd = fdList.get(0);
		Assert.assertEquals(fd.getGene(), genes[1]);
		Assert.assertEquals(fd.getStrand(), strands[1]);
		Assert.assertEquals(fd.getLocusFunction(), locusFunctions[0]);

		fdList = fdp.filterToPreferredAnnotations(fdList);
		Assert.assertEquals(fdList.size(), 1);

	}

	//Desired CODING reads, have CODING read on SENSE and ANTISENESE.  Test both.
	@Test
	public void testStrandStrategy () {
		LocusFunction [] acceptedFunctions = {LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC};
		AbstractFunctionalDataProcessor fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.SENSE, acceptedFunctions);
		String [] genes = {"A", "B"};
		String [] strands = {"+", "-"};
		LocusFunction [] locusFunctions = {LocusFunction.CODING, LocusFunction.CODING};
		// SENSE STRAND TEST positive strand
		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, false);
		Assert.assertEquals(fdList.size(), 1);
		FunctionalData fd = fdList.get(0);
		Assert.assertEquals(fd.getGene(), genes[0]);
		Assert.assertEquals(fd.getStrand(), strands[0]);
		Assert.assertEquals(fd.getLocusFunction(), locusFunctions[0]);

		// SENSE STRAND TEST negative strand
		fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, true);
		Assert.assertEquals(fdList.size(), 1);
		fd = fdList.get(0);
		Assert.assertEquals(fd.getGene(), genes[1]);
		Assert.assertEquals(fd.getStrand(), strands[1]);
		Assert.assertEquals(fd.getLocusFunction(), locusFunctions[1]);

		// ANTISENSE STRAND TESTs
		fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.ANTISENSE, acceptedFunctions);

		// ANTISENSE STRAND TEST negative strand
		fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, true);
		Assert.assertEquals(fdList.size(), 1);
		fd = fdList.get(0);
		Assert.assertEquals(fd.getGene(), genes[0]);
		Assert.assertEquals(fd.getStrand(), strands[0]);
		Assert.assertEquals(fd.getLocusFunction(), locusFunctions[0]);

		// ANTISENSE SENSE STRAND TEST positive strand
		fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, false);
		Assert.assertEquals(fdList.size(), 1);
		fd = fdList.get(0);
		Assert.assertEquals(fd.getGene(), genes[1]);
		Assert.assertEquals(fd.getStrand(), strands[1]);
		Assert.assertEquals(fd.getLocusFunction(), locusFunctions[1]);

		// BOTH returns both A and B.
		fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.BOTH, acceptedFunctions);
		fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, true);
		Assert.assertEquals(fdList.size(), 2);
		Assert.assertEquals(fdList.get(0).getGene(), genes[0]);
		Assert.assertEquals(fdList.get(1).getGene(), genes[1]);

	}

	// For a gene that has a read that's intronic and coding, should collapse to the preferred (CODING) even though both are on the initial accepted list.
	@Test
	public void testFunctionCollapse () {
		LocusFunction [] acceptedFunctions = {LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC};
		AbstractFunctionalDataProcessor fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.SENSE, acceptedFunctions);
		String [] genes = {"A", "A"};
		String [] strands = {"+", "+"};
		LocusFunction [] locusFunctions = {LocusFunction.CODING, LocusFunction.INTRONIC};

		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, false);
		Assert.assertEquals(fdList.size(), 1);
		FunctionalData fd = fdList.get(0);
		Assert.assertEquals(fd.getGene(), genes[0]);
		Assert.assertEquals(fd.getStrand(), strands[0]);
		Assert.assertEquals(fd.getLocusFunction(), locusFunctions[0]);

		fdList = fdp.filterToPreferredAnnotations(fdList);
		Assert.assertEquals(fdList.size(), 1);
	}

	@Test
	// for two genes, one intronic, one coding.
	// at first two genes should be emitted.
	// after filtering for the preferred type, only the coding read should be emitted.
	public void testPreferredAnnotationTypeCodingIntronic () {
		LocusFunction [] acceptedFunctions = {LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC};
		AbstractFunctionalDataProcessor fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.SENSE, acceptedFunctions);
		String [] genes = {"A", "B"};
		String [] strands = {"+", "+"};
		LocusFunction [] locusFunctions = {LocusFunction.CODING, LocusFunction.INTRONIC};

		// both genes come out
		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, false);
		Assert.assertEquals(fdList.size(), 2);
		FunctionalData fd = fdList.get(0);
		Assert.assertEquals(fd.getGene(), genes[0]);
		Assert.assertEquals(fd.getStrand(), strands[0]);
		Assert.assertEquals(fd.getLocusFunction(), locusFunctions[0]);
		FunctionalData fd2 = fdList.get(1);
		Assert.assertEquals(fd2.getGene(), genes[1]);
		Assert.assertEquals(fd2.getStrand(), strands[1]);
		Assert.assertEquals(fd2.getLocusFunction(), locusFunctions[1]);

		// once you filter to the preferred annotations, coding is "the best", so you only have a single FunctionalData.
		fdList = fdp.filterToPreferredAnnotations(fdList);
		Assert.assertEquals(fdList.size(), 1);
		fd = fdList.get(0);
		Assert.assertEquals(fd.getGene(), genes[0]);
		Assert.assertEquals(fd.getStrand(), strands[0]);
		Assert.assertEquals(fd.getLocusFunction(), locusFunctions[0]);
	}

	@Test
	// for two genes, both coding.
	// at first two genes should be emitted.
	// after filtering for the preferred type, both are emitted.  Sad days, this doesn't resolve ambiguity.
	public void testPreferredAnnotationTypeCodingCoding () {
		LocusFunction [] acceptedFunctions = {LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC};
		AbstractFunctionalDataProcessor fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.SENSE, acceptedFunctions);
		String [] genes = {"A", "B"};
		String [] strands = {"+", "+"};
		LocusFunction [] locusFunctions = {LocusFunction.CODING, LocusFunction.CODING};

		// both genes come out
		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, false);
		Assert.assertEquals(fdList.size(), 2);
		FunctionalData fd = fdList.get(0);
		Assert.assertEquals(fd.getGene(), genes[0]);
		Assert.assertEquals(fd.getStrand(), strands[0]);
		Assert.assertEquals(fd.getLocusFunction(), locusFunctions[0]);
		FunctionalData fd2 = fdList.get(1);
		Assert.assertEquals(fd2.getGene(), genes[1]);
		Assert.assertEquals(fd2.getStrand(), strands[1]);
		Assert.assertEquals(fd2.getLocusFunction(), locusFunctions[1]);

		// once you filter to the preferred annotations, both are still present.
		fdList = fdp.filterToPreferredAnnotations(fdList);
		Assert.assertEquals(fdList.size(), 2);
		fd = fdList.get(0);
		Assert.assertEquals(fd.getGene(), genes[0]);
		Assert.assertEquals(fd.getStrand(), strands[0]);
		Assert.assertEquals(fd.getLocusFunction(), locusFunctions[0]);
		fd2 = fdList.get(1);
		Assert.assertEquals(fd2.getGene(), genes[1]);
		Assert.assertEquals(fd2.getStrand(), strands[1]);
		Assert.assertEquals(fd2.getLocusFunction(), locusFunctions[1]);
	}

	// gf:Z:CODING,INTRONIC,UTR	gn:Z:CHKB-CPT1B,CHKB-CPT1B,CPT1B	gs:Z:-,-,-
	@Test
	public void testFromRealData1() {
		LocusFunction [] acceptedFunctions = {LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC};
		AbstractFunctionalDataProcessor fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.SENSE, acceptedFunctions);
		String [] genes = {"CHKB-CPT1B", "CHKB-CPT1B", "CPT1B"};
		String [] strands = {"-", "-", "-"};
		LocusFunction [] locusFunctions = {LocusFunction.CODING, LocusFunction.INTRONIC, LocusFunction.UTR};

		// both genes come out.
		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, true);
		Assert.assertEquals(fdList.size(), 2);

		// once you filter to the preferred annotations, both are still present.
		fdList = fdp.filterToPreferredAnnotations(fdList);
		Assert.assertEquals(fdList.size(), 2);

	}

	// Demonstrate interaction of a coding and intronic gene on the same read, with coding+intronic interpretation.
	@Test
	public void testExample1() {
		LocusFunction [] acceptedFunctions = {LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC};
		AbstractFunctionalDataProcessor fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.SENSE, acceptedFunctions);
		String [] genes = {"A", "A", "B", "C"};
		String [] strands = {"+", "+", "+", "-"};
		LocusFunction [] locusFunctions = {LocusFunction.CODING, LocusFunction.INTRONIC, LocusFunction.INTRONIC, LocusFunction.INTRONIC};

		// Read on the + strand.
		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, false);
		Assert.assertEquals(fdList.size(), 2);

		// once you filter to the preferred annotations, only the coding gene is retained.
		
		fdList = fdp.filterToPreferredAnnotations(fdList);
		Assert.assertEquals(fdList.size(), 1);

	}

	// Demonstrate interaction of a coding and intronic gene on the same read, with intronic interpretation.
	@Test
	public void testExample2() {
		
		LocusFunction [] acceptedFunctions = {LocusFunction.INTRONIC};
		AbstractFunctionalDataProcessor fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.SENSE, acceptedFunctions);
		String [] genes = {"A", "A", "B", "C"};
		String [] strands = {"+", "+", "+", "-"};
		LocusFunction [] locusFunctions = {LocusFunction.CODING, LocusFunction.INTRONIC, LocusFunction.INTRONIC, LocusFunction.INTRONIC};

		// Read on the + strand.
		// Because the read overlaps Gene A on both the intronic and exonic portions, the gene annotation is discarded.
		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, false);
		Assert.assertEquals(fdList.size(), 1);

		// once you filter to the preferred annotations, only the coding gene is retained.		
		fdList = fdp.filterToPreferredAnnotations(fdList);
		Assert.assertEquals(fdList.size(), 1);

	}

	// Demonstrate interaction of a coding and intronic gene on the same read, with coding interpretation.
	@Test
	public void testExample3() {
		
		LocusFunction [] acceptedFunctions = {LocusFunction.CODING, LocusFunction.UTR};
		AbstractFunctionalDataProcessor fdp = new DropSeqFunctionalDataProcessor(StrandStrategy.SENSE, acceptedFunctions);
		String [] genes = {"A", "A", "B", "C"};
		String [] strands = {"+", "+", "+", "-"};
		LocusFunction [] locusFunctions = {LocusFunction.CODING, LocusFunction.INTRONIC, LocusFunction.INTRONIC, LocusFunction.INTRONIC};

		// Read on the + strand.
		// Because the read overlaps Gene A on both the intronic and exonic portions, the gene annotation is discarded.
		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, false);
		Assert.assertEquals(fdList.size(), 0);

		// once you filter to the preferred annotations, only the coding gene is retained.		
		fdList = fdp.filterToPreferredAnnotations(fdList);
		Assert.assertEquals(fdList.size(), 0);

	}

}
