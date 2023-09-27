package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeIterator.DgeLine;
import org.testng.Assert;
import org.testng.annotations.Test;

import com.google.common.collect.Sets;

public class DgeIteratorTest {

	private final File inFile = new File("testdata/org/broadinstitute/transcriptome/barnyard/digitalexpression/test_with_header.dge.txt.gz");
	@Test
	public void testParse() {
		DgeIterator iter = new DgeIterator(this.inFile);
		DgeHeader h = iter.getDgeHeader();
		Assert.assertTrue(iter.hasNext());
		DgeLine l = iter.next();
		String gene = l.getGene();
		double [] exp = l.getExpression();
		List<String> ids = new ArrayList<>(l.getIdentifiers());

		Assert.assertTrue(h!=null);

		// header: GENE    GATCGATAGAAACCAT        GGATTACTCATTATCC        TCGAGGCTCAGCCTAA        AACCATGCACATCTTT        CTAATGGCAATACGCT        ACTTACTCAGGCGATA        AGGTCCGCATGTAGTC        TTTGCGCAGCAACGGT        ATAACGCCACTTCTGC        CCACCTAGTGTCCTCT
		// first line expected: A4GALT  2       2       1       0       0       3       0       0       0       0

		// validate cell barcode order in header.
		String [] cellBarcodes = {"GATCGATAGAAACCAT","GGATTACTCATTATCC","TCGAGGCTCAGCCTAA","AACCATGCACATCTTT","CTAATGGCAATACGCT","ACTTACTCAGGCGATA","AGGTCCGCATGTAGTC","TTTGCGCAGCAACGGT","ATAACGCCACTTCTGC","CCACCTAGTGTCCTCT"};
		for (int j=0; j<cellBarcodes.length; j++)
			Assert.assertEquals(cellBarcodes[j], ids.get(j));

		// validate first line.
		Assert.assertEquals("A4GALT", gene);
		int [] expectedExpression = {2,2,1,0,0,3,0,0,0,0};
		Assert.assertEquals(expectedExpression.length, exp.length);

		for (int j=0; j<expectedExpression.length; j++)
			Assert.assertEquals(expectedExpression[j], exp[j], 0.0001);

		// validate random access.
		for (int j=0; j<cellBarcodes.length; j++)
			Assert.assertEquals(expectedExpression[j], l.getExpression(cellBarcodes[j]), 0.0001);

		// validate 10 lines in file.
		int counter=1;
		while (iter.hasNext()) {
			iter.next();
			counter++;
		}
		Assert.assertEquals(10,  counter);

		//check close.
		iter.close();

	}

	public void testSubset () {
		Set<String> cellBarcodes = new HashSet<>(Arrays.asList("GGATTACTCATTATCC","TCGAGGCTCAGCCTAA","CTAATGGCAATACGCT","AGGTCCGCATGTAGTC","TTTGCGCAGCAACGGT","CCACCTAGTGTCCTCT"));
		DgeIterator iter = new DgeIterator(this.inFile);

		while (iter.hasNext()) {
			DgeLine l = iter.next();
			DgeLine filtered = l.subset(cellBarcodes);
			// check that the filtered identifiers are a match for the cell barcodes list with no extra values, and that
			// the values in the filtered line match the original line.
			Set<String> intersect = Sets.intersection(cellBarcodes, filtered.getIdentifiers());
			Assert.assertTrue(intersect.size()==cellBarcodes.size());
			for (String k : filtered.getIdentifiers()) {
				double valueOld = l.getExpression(k);
				double valueNew = filtered.getExpression(k);
				Assert.assertEquals(valueOld, valueNew, 0.00001);
			}
		}

		iter.close();


	}
}
