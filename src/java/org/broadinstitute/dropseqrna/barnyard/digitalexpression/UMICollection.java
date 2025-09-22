/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import java.io.File;
import java.util.*;

import org.broadinstitute.dropseqrna.barnyard.GatherMolecularBarcodeDistributionByGene;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.editdistance.MapBarcodesByEditDistance;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.PeekableIterator;
import picard.sam.util.Pair;
import picard.util.TabbedTextFileWithHeaderParser;

/**
 * Models a collection of UMIs for a gene and cell barcode.
 * @author nemesh
 *
 */
public class UMICollection {

	private String cellBarcode;
	private String geneName;
	private ObjectCounter<String> molecularBarcodeCounts;	
	private static MapBarcodesByEditDistance mbed =new MapBarcodesByEditDistance(false);
	private Map<String, List<SAMRecord>> reads=null;

	public UMICollection (final String cellBarcode, final String geneName) {
		this.cellBarcode = cellBarcode;
		this.geneName = geneName;
		molecularBarcodeCounts=new ObjectCounter<>();
	}

	public void incrementMolecularBarcodeCount (final String molecularBarcode, final int count) {
		molecularBarcodeCounts.incrementByCount(molecularBarcode, count);
	}

	public void incrementMolecularBarcodeCount (final String molecularBarcode) {
		molecularBarcodeCounts.incrementByCount(molecularBarcode, 1);
	}

	public String getCellBarcode() {
		return cellBarcode;
	}

	public String getGeneName() {
		return geneName;
	}

	public ObjectCounter<String> getMolecularBarcodeCounts() {
		return this.molecularBarcodeCounts;
	}

	public Collection<String> getMolecularBarcodes() {
		return this.molecularBarcodeCounts.getKeys();
	}

	public boolean isEmpty () {
		return this.molecularBarcodeCounts.getSize()==0;
	}

	public ObjectCounter<String> getMolecularBarcodeCountsCollapsed(final int editDistance) {
		ObjectCounter<String> counts = collapseByEditDistance(this.molecularBarcodeCounts, editDistance);
		return counts;
	}
	
	public void addRead (String molecularBarcode, SAMRecord rec) {
		if (reads==null) {
			reads=new HashMap<>();
		}
		List<SAMRecord> l = reads.get(molecularBarcode);
		if (l==null) {
			l = new ArrayList<>();
			reads.put(molecularBarcode, l);
		}
		l.add(rec);
	}
	
	/**
	 * If reads were added to this collection via addRead, they can be retrieved for a given molecular barcode via this method.
	 * @param molecularBarcode The molecular barcode to request reads from
	 * @return The list of reads for this molecular barcode.  This will return null when there were no reads added for a requested barcode.
	 */
	public List<SAMRecord> getReads (String molecularBarcode) {
		if (this.reads==null)
			return null;
		return this.reads.get(molecularBarcode);		
	}

	/**
	 * Take this data set and filter it by the frequency of the UMIs observed across the cell/gene pair.
	 * Calculate the mean #observations across the UMIs, the drop UMIs with less than the mean * the frequency.
	 * @param freq
	 */
	public void filterByUMIFrequency (final double freq) {
		int total = 0;
		double average = this.molecularBarcodeCounts.getTotalCount() / this.molecularBarcodeCounts.getSize();
		int minNumReads = (int) Math.floor(average * freq);
		this.molecularBarcodeCounts.filterByMinCount(minNumReads);
	}


	/**
	 * Generate digital expression output for a cell and gene combo (which this object handily tracks.)
	 * This method allows barcode collapse of molecular barcodes by varying edit distance, and can also count reads instead of UMIs if desired.
	 * @param editDistance How close two molecular barcodes should be in edit distance space to be collapsed
	 * @param minBCReadThreshold The minimum number of reads a molecular barcode should have to be considered.  This is done AFTER edit distance collapse of barcodes.
	 * @param outputReads Instead of counting UMIs, count the number of reads on those UMIs.
	 * @return
	 */
	public int getDigitalExpression (final int minBCReadThreshold, final int editDistance, final boolean outputReads) {
		// because this is confusing, I've hard coded the distance one should naively search for edit distance to 3*edit distance.
		// this number has seemed reasonable for all our practical tests.
		// int threshold = editDistance*3;

		if (outputReads) {
			int count = this.molecularBarcodeCounts.getTotalCount();
			return (count);
		}
		// harder, collapse molecular barcodes and count them.
		ObjectCounter<String> counts = collapseByEditDistance(this.molecularBarcodeCounts, editDistance);
		counts.filterByMinCount(minBCReadThreshold);
		int count = counts.getKeys().size();
		return count;
	}

	/**
	 * For a list of molecular barcodes, collapse them by edit distance.
	 * @param counts
	 * @param editDistance
	 * @return A new object counter with the molecular barcodes collapsed by edit distance.
	 */
	private ObjectCounter<String> collapseByEditDistance (final ObjectCounter<String> counts, final int editDistance) {
		ObjectCounter<String> result = mbed.collapseAndMergeBarcodes(counts, false, editDistance);
		return (result);
	}

	/**
	 * Collapse the molecular barcodes in this collection by edit distance.  This only makes sense if reads != null,
	 * because otherwise, one could just call getMolecularBarcodeCountsCollapsed to get collapsed read counts.
	 * The reads are re-tagged with the collapsed molecular barcode, and re-grouped, and the molecular barcode counts are updated.
	 * @param editDistance
	 * @param molecularBarcodeTag
	 */
	public void collapseThisByEditDistance(final int editDistance, final String molecularBarcodeTag) {
		if (reads == null) {
			throw new IllegalStateException("It doesn't make sense to collapse in place without reads");
		}
		if (editDistance < 1) {
			throw new IllegalArgumentException("Edit distance must be at least 1");
		}
		Map<String, List<String>> collapsedToUncollapsed = mbed.collapseBarcodes(this.molecularBarcodeCounts, false, editDistance);
		for (Map.Entry<String, List<String>> entry: collapsedToUncollapsed.entrySet()) {
			final String goodUmi = entry.getKey();
			final List<String> umisToBeCollapsed = entry.getValue();
			final List<SAMRecord> readsForGoodUmi = reads.get(goodUmi);
			for (final String umiToBeCollapsed : umisToBeCollapsed) {
				final List<SAMRecord> readsToBeCollapsed = reads.get(umiToBeCollapsed);
				readsToBeCollapsed.stream().forEach(r -> r.setAttribute(molecularBarcodeTag, goodUmi));
				readsForGoodUmi.addAll(readsToBeCollapsed);
				reads.remove(umiToBeCollapsed);
			}
		}
		molecularBarcodeCounts = new ObjectCounter<>();
		for (Map.Entry<String, List<SAMRecord>> entry: this.reads.entrySet()) {
			molecularBarcodeCounts.setCount(entry.getKey(), entry.getValue().size());
		}
	}

	public static Collection<UMICollection> parseUMICollectionFile (final File input) {
		IOUtil.assertFileIsReadable(input);

		Collection<UMICollection> result = new ArrayList<>();
		TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(input);

		UMICollection currentUMI = null;

		PeekableIterator<TabbedTextFileWithHeaderParser.Row> parserIter = new PeekableIterator<>(parser.iterator());

		if (parserIter.hasNext()) {
			TabbedTextFileWithHeaderParser.Row row  =parserIter.peek();
			String cell = row.getField(GatherMolecularBarcodeDistributionByGene.CELL_BARCODE_COLUMN);
			String gene = row.getField(GatherMolecularBarcodeDistributionByGene.GENE_COLUMN);
			currentUMI = new UMICollection(cell, gene);
		}

		//when a new gene/cell is seen, make a new object and put records into that, and store the old gene/cell.
		while(parserIter.hasNext()) {
			TabbedTextFileWithHeaderParser.Row row =parserIter.next();
			String cell = row.getField(GatherMolecularBarcodeDistributionByGene.CELL_BARCODE_COLUMN);
			String gene = row.getField(GatherMolecularBarcodeDistributionByGene.GENE_COLUMN);
			String molBC = row.getField(GatherMolecularBarcodeDistributionByGene.MOLECULAR_BARCODE_COLUMN);
			int count = Integer.parseInt(row.getField(GatherMolecularBarcodeDistributionByGene.NUM_OBS_COLUMN));
			// if you're in not in the same cell/gene, add the current UMI to the collection and make a new one.
			if (!cell.equals(currentUMI.cellBarcode) || !gene.equals(currentUMI.geneName)) {
				result.add(currentUMI);
				currentUMI = new UMICollection(cell, gene);
			}
			currentUMI.incrementMolecularBarcodeCount(molBC, count);
		}
		if (currentUMI.molecularBarcodeCounts.getSize()>0)
			result.add(currentUMI);
		CloserUtil.close(parserIter);
		return (result);
	}

	@Override
	public String toString () {
		return "[" + this.cellBarcode + "] + [" + this.geneName + "] " + this.molecularBarcodeCounts.toString();
	}


	/**
	 * Estimate number of reads for each UMI at a particular downsampling rate.
	 * @author dmeyer
	 * @param downsampleRate
	 * @param random
     * @param minDownsampledCount Minimum number of reads a UMI must have after downsampling to be retained.
	 * @return List of pairs of UMI and estimated downsampled read count.
	 */
	public List<Pair<String, Integer>> getDownsampledMolecularBarcodeCounts(double downsampleRate,
                                                                            SplittableRandom random,
                                                                  int minDownsampledCount) {
		assert downsampleRate >= 0 && downsampleRate <= 1;

		// Iterate over umis create a downsampled objectCounter that has umi counts at a given downsampleRate
        List<Pair<String, Integer>> res = new ArrayList<>(this.molecularBarcodeCounts.getSize());
		int count; int downsampledCount;
		for (String umi : this.molecularBarcodeCounts.getKeys()) {
			count = this.getMolecularBarcodeCounts().getCountForKey(umi);
			downsampledCount =  (int)random.doubles(count).filter(r -> r < downsampleRate).count();
            if (downsampledCount >= minDownsampledCount) {
                res.add(new Pair<>(umi, downsampledCount));
            }
		}
		return res;
	}

}
