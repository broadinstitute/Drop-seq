package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTagUtil;
import htsjdk.samtools.util.*;
import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.MultiComparator;

import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

public class UMIIterator implements CloseableIterator<UMICollection>  {

    private static final Log log = Log.getInstance(UMIIterator.class);
    private static final ProgressLogger prog = new ProgressLogger(log);

	private final GroupingIterator<SAMRecord> atoi;
	private final String geneExonTag;
	private final String cellBarcodeTag;
	private final String molecularBarcodeTag;
	
	/**
	 * Construct an object that generates UMI objects from a BAM file
     * @param headerAndIterator The BAM records to extract UMIs from
	 * @param geneExonTag The geneExon tag on BAM records
	 * @param cellBarcodeTag The cell barcode tag on BAM records
	 * @param molecularBarcodeTag The molecular barcode tag on BAM records
	 * @param strandTag The strand tag on BAM records
	 * @param readMQ The minimum map quality of the reads
	 * @param assignReadsToAllGenes Should records tagged with multiple genes be double counted, once for each gene?
	 * @param useStrandInfo should the gene and read strand match for the read to be accepted
	 * @param cellBarcodes The list of cell barcode tag values that match the <cellBarcodeTag> tag on the BAM records.
     *                     Only reads with these values will be used.  If set to null, all cell barcodes are used.
	 */
	public UMIIterator(SamHeaderAndIterator headerAndIterator,
                       String geneExonTag,
                       String cellBarcodeTag,
                       String molecularBarcodeTag,
                       String strandTag,
                       int readMQ,
                       boolean assignReadsToAllGenes,
                       boolean useStrandInfo,
                       List<String> cellBarcodes) {
		this(headerAndIterator, geneExonTag, cellBarcodeTag, molecularBarcodeTag, strandTag, readMQ, assignReadsToAllGenes, useStrandInfo, cellBarcodes, false);
	}
	
	/**
	 * Construct an object that generates UMI objects from a BAM file
	 * @param headerAndIterator The BAM records to extract UMIs from
	 * @param geneExonTag The geneExon tag on BAM records
	 * @param cellBarcodeTag The cell barcode tag on BAM records
	 * @param molecularBarcodeTag The molecular barcode tag on BAM records
	 * @param strandTag The strand tag on BAM records
	 * @param readMQ The minimum map quality of the reads
	 * @param assignReadsToAllGenes Should records tagged with multiple genes be double counted, once for each gene?
	 * @param useStrandInfo should the gene and read strand match for the read to be accepted
	 * @param cellBarcodes The list of cell barcode tag values that match the <cellBarcodeTag> tag on the BAM records.
     *                     Only reads with these values will be used.  If set to null, all cell barcodes are used.
	 * @param cellFirstSort if true, then cell barcodes are sorted first, followed by gene/exon tags.
     *                      If false, then gene/exon tags are sorted first, followed by cells.  false is the default and used in the other constructor.
	 */
	public UMIIterator(SamHeaderAndIterator headerAndIterator,
                       String geneExonTag,
                       String cellBarcodeTag,
                       String molecularBarcodeTag,
                       String strandTag,
                       int readMQ,
                       boolean assignReadsToAllGenes,
                       boolean useStrandInfo,
                       List<String> cellBarcodes,
                       boolean cellFirstSort) {

        this.geneExonTag=geneExonTag;
        this.cellBarcodeTag=cellBarcodeTag;
        this.molecularBarcodeTag=molecularBarcodeTag;

        final StringTagComparator cellBarcodeTagComparator = new StringTagComparator(cellBarcodeTag);
        final StringTagComparator geneExonTagComparator = new StringTagComparator(geneExonTag);
        final MultiComparator<SAMRecord> multiComparator;
		if (cellFirstSort) {
            multiComparator = new MultiComparator<>(cellBarcodeTagComparator, geneExonTagComparator);
		} else {
            multiComparator = new MultiComparator<>(geneExonTagComparator, cellBarcodeTagComparator);
		}
        // Filter records before sorting, to reduce I/O
		final MissingTagFilteringIterator filteringIterator =
                new MissingTagFilteringIterator(headerAndIterator.iterator, cellBarcodeTag, geneExonTag);

        SortingCollection<SAMRecord> alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
                new BAMRecordCodec(headerAndIterator.header), multiComparator,
                SAMFileWriterImpl.getDefaultMaxRecordsInRam());

        log.info("Reading in records for TAG name sorting");
        while (filteringIterator.hasNext()) {
            final SAMRecord rec = filteringIterator.next();
            prog.record(rec);
            alignmentSorter.add(rec);
        }
        filteringIterator.close();


        // assign the tags in the order you want data sorted.
		
		UmiIteratorWrapper umiIteratorWrapper = new UmiIteratorWrapper(alignmentSorter.iterator(), cellBarcodeTag,
                cellBarcodes, geneExonTag, strandTag, readMQ, assignReadsToAllGenes, useStrandInfo);
        // Not really -- merge sort is ongoing.
        log.info("Sorting finished.");

		this.atoi = new GroupingIterator<>(umiIteratorWrapper, multiComparator);
	}
	
	private static class StringTagComparator implements Comparator<SAMRecord> {
        private final short tag;

        StringTagComparator(String tag) {
            this.tag = SAMTagUtil.getSingleton().makeBinaryTag(tag);
        }

        @Override
        public int compare(SAMRecord rec1, SAMRecord rec2) {
            final String s1 = (String)rec1.getAttribute(tag);
            final String s2 = (String)rec2.getAttribute(tag);

            if (s1 != null) {
                if (s2 == null)
                    return 1;
                else {
                    return s1.compareTo(s2);
                }
            } else if (s2 != null) {
                return -1;
            } else {
                return 0;
            }
        }
    }

    private static class MissingTagFilteringIterator extends FilteredIterator<SAMRecord> {
        final short[] requiredTags;

        private MissingTagFilteringIterator(final Iterator<SAMRecord> underlyingIterator, final String...requiredTags) {
            super(underlyingIterator);
            this.requiredTags = new short[requiredTags.length];
            for (int i = 0; i < requiredTags.length; ++i) {
                this.requiredTags[i] = SAMTagUtil.getSingleton().makeBinaryTag(requiredTags[i]);
            }
        }

        @Override
        protected boolean filterOut(SAMRecord rec) {
            for (final short tag : requiredTags) {
                if (rec.getAttribute(tag) == null) {
                    return true;
                }
            }
            return false;
        }
    }
	
	/**
	 * Gets the next UMI Collection - all the molecular barcodes and reads on those for a particular gene/cell.
	 * @return Null if there are no reads left in the iterator.  Otherwise, returns a UMICollection.
	 */
	@Override
	public UMICollection next () {
		if (!this.atoi.hasNext()) {
			return null;
		}
		
		Collection<SAMRecord> records = this.atoi.next();
		PeekableIterator<SAMRecord> recordCollectionIter = new PeekableIterator<SAMRecord>(records.iterator());
		
		// the next record is the first of the "batch"
		SAMRecord r = recordCollectionIter.peek();
		String currentGene = r.getStringAttribute(geneExonTag);
		String currentCell = Utils.getCellBC(r, cellBarcodeTag);
		UMICollection umi = new UMICollection(currentCell, currentGene);		
		
		while (recordCollectionIter.hasNext()) {
			// if there's a next read, set up the current gene/cell variables.
			r=recordCollectionIter.next();
			String molecularBarcode = r.getStringAttribute(this.molecularBarcodeTag);
			umi.incrementMolecularBarcodeCount(molecularBarcode);
						
		}
		// I ran out of reads
		recordCollectionIter.close();
		return (umi);
	}
	
	@Override
	public void remove() {
		this.atoi.remove();
	}

	@Override
	public void close() {
        CloserUtil.close(this.atoi);
	}

	@Override
	public boolean hasNext() {
		return this.atoi.hasNext();
	}
	
	
}
