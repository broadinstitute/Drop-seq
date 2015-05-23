package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.PeekableIterator;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.readIterators.DEIteratorUtils;
import org.broadinstitute.dropseqrna.utils.readIterators.UMIIterator;
import org.testng.Assert;
import org.testng.annotations.Test;


public class DigitalExpressionTest {

	
	// generate an input file with 5 barcodes and 3 genes.
	/**
	 * For the integration test, make a relatively small BAM file and read it in.  
	 * Manually test the counts of the reads across UMIs - umi testing to collapse molecular barcodes is in it's own test.
	 * 
	 * Data taken from 9-27-14 100 cells data set.
	 * 
	 * cells
	 * ATCAGGGACAGA
	 * AGGGAAAATTGA
	 * TTGCCTTACGCG
	 * TGGCGAAGAGAT
	 * TACAATTAAGGC
	 * 
	 * genes
	 * HUMAN_10:101948055-101989376:CHUK
	 * HUMAN_15:101821715-101835487:SNRPA1
	 * HUMAN_3:42642106-42690227:NKTR
	 * 
	 * /fg/software/gap/gap_analysis/FilterBAMByTag I=100cells_star_bq10_noPseudoGenes.bam O=test.bam TAG=ZC TAG_VALUES_FILE=bc.txt ACCEPT_TAG=true
	 * samtools view -H test.bam > 5cell3gene.sam
	 * samtools view test.bam |grep -f genes.txt >> 5cell3gene.sam 
	 * samtools view -Sb 5cell3gene.sam > 5cell3gene.bam
	 * 
	 */
		
	File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.bam");
	
	private int MAX_RECORDS_IN_RAM=100000;
	private String GENE_EXON_TAG="GE";
	private String STRAND_TAG="GS";
	private String CELL_BARCODE_TAG="ZC";
	private String MOLECULAR_BARCODE_TAG = "XM";
	private int READ_MQ=10;
	private Boolean USE_STRAND_INFO=true;
	
	
	private UMIIterator getUMIIterator () {
		String [] barcodes ={"ATCAGGGACAGA", "AGGGAAAATTGA", "TTGCCTTACGCG", "TGGCGAAGAGAT", "TACAATTAAGGC"};
		List<String> cellBarcodes = Arrays.asList(barcodes);
		
		UMIIterator iter = new UMIIterator(this.IN_FILE, this.GENE_EXON_TAG, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.STRAND_TAG, this.READ_MQ, 
				true, this.USE_STRAND_INFO, cellBarcodes, this.MAX_RECORDS_IN_RAM);
		
		return (iter);
	}
	
	
	@Test(groups={"dropseq", "transcriptome"})
	public void DGEIntegrationTest () {
		UMIIterator u = getUMIIterator ();
		UMICollection batch;
		while ((batch=u.next())!=null) {
			if (batch==null || batch.isEmpty()){
				continue;	
			}
			String currentGene = batch.getGeneName();
			String currentCell = batch.getCellBarcode();
			int expectedReadcount = getReadCounts(currentGene, currentCell);
			int dgeReadCount = batch.getDigitalExpression(1, 0, true);
			Assert.assertEquals(dgeReadCount,expectedReadcount);
			
			int dgeNoCollapseExpected=getDGEWithoutCollapse(currentGene, currentCell);
			int dgeNoCollapseActual = batch.getDigitalExpression(1, 0, false);
			Assert.assertEquals(dgeNoCollapseActual,dgeNoCollapseExpected);
			
		}
		
	}
	
	//GET COUNTS OF READS
	// ~/samtools view -q 10 5cell3gene.bam |grep ZC:Z:ATCAGGGACAGA |grep HUMAN_3:42642106-42690227:NKTR  | sed -E 's/.*XM
	private int getReadCounts (String gene, String cell) {
	
		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 845;
		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 200;
		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 225;
		
		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 428;
		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 212;
		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 238;
		
		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 473;
		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 62;
		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 581;
		
		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 612;
		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 12;
		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 166;
		
		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 385;
		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 160;
		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 294;
		return -1;
		
	}
	
	// GET DGE (WITHOUT BC COLLAPSE)
	// ~/samtools view -q 10 5cell3gene.bam |grep ZC:Z:ATCAGGGACAGA |grep HUMAN_3:42642106-42690227:NKTR  | sed -E 's/.*XM:Z:([ACGT]+).*/\1/' |sort |uniq -c |wc -l
	private int getDGEWithoutCollapse (String gene, String cell) {
		
		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 84;
		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 19;
		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 34;
		
		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 54;
		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 14;
		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 45;
		
		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 40;
		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 4;
		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 92;
		
		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 59;
		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 2;
		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 32;
		
		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 40;
		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 11;
		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 32;
		return -1;
		
	}
	
	
	
	
}
