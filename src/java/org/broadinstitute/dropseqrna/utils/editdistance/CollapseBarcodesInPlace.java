package org.broadinstitute.dropseqrna.utils.editdistance;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.*;
import org.broadinstitute.dropseqrna.barnyard.BarcodeListRetrieval;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.metrics.BAMTagHistogram;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Fold down barcodes, possibly in the context of another barcode (that has been folded down already.) 
 * @author nemesh
 *
 */
@CommandLineProgramProperties(usage = "Fold down barcodes, possibly in the context of another barcode (that has been folded down already.)",
        usageShort = "Fold down barcodes, possibly in the context of another barcode (that has been folded down already.)",
        programGroup = DropSeq.class)
public class CollapseBarcodesInPlace extends CommandLineProgram {


	private final Log log = Log.getInstance(CollapseBarcodesInPlace.class);
	private ProgressLogger pl = new ProgressLogger(this.log);
	
	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.  Must be coordinate sorted. ",
            minElements = 1)
	public List<File> INPUT;
	
	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BAM file with an extra tag.")
	public File OUTPUT;
	
	@Option(doc="Barcode to collapse")
	public String PRIMARY_BARCODE;
	
	@Option(doc="The edit distance to collapse barcodes")
	public Integer EDIT_DISTANCE=1;
	
	@Option(doc = "Should indels be considered in edit distance calculations?  Doing this correctly is far slower than a simple edit distance test, but gives a more complete result.")
	public boolean FIND_INDELS=true;
	
	@Option(doc="The output barcode tag for the newly collapsed barcodes")
	public String OUT_BARCODE;
	
	@Option(doc="Read quality filter.  Filters all reads lower than this mapping quality.  Defaults to 10.  Set to 0 to not filter reads by map quality.")
	public Integer READ_QUALITY=10;
	
	@Option(doc="Number of reads a barcode would need to have in order to have other barcodes get merged into it.  All barcodes are candidates to be merged into another barcode." +
			"For cell barcodes you probably want to set this to a relatively high number like 100, since we expect cells to have thousands or more reads, and this signficantly speeds up analysis.  " +
			"For molecular barcodes, you probably want to set this to 1, as you want to include all molecular barcodes, unless you have very high sequencing depth.", optional=true)
	public Integer MIN_NUM_READS_CORE=null;
	
	@Option(doc="Number of cells that you think are in the library.  This accomplishes the same goals as the MIN_NUM_READS_CORE argument, but instead of defining barcodes as important based on the number of reads, it picks the top <X> barcodes as core.", optional=true)
	public Integer NUM_CORE_BARCODES=null;
	
	@Option(doc="The number of reads a non-core barcode must have to be merged with a core barcode.", optional=true)
	public Integer MIN_NUM_READS_NONCORE=1;
	
	@Option(doc="Filter PCR Duplicates.  Defaults to false")
	public boolean FILTER_PCR_DUPLICATES=false;
	
	@Option(doc="Number of threads to use.  Defaults to 1.")
	public int NUM_THREADS=1;
	
	
	private int REPORT_PROGRESS_INTERVAL=100;
	private CollapseBarcodeThreaded cbt=null;
	private int threadedBlockSize=20000;
	
	@Override
	protected int doWork() {
		log.info("Number of cores selected [" + Integer.toString(this.NUM_THREADS) + "]"); 
		
		if (this.NUM_THREADS>1) cbt= new CollapseBarcodeThreaded(this.threadedBlockSize, this.NUM_THREADS);
		
		IOUtil.assertFileIsWritable(OUTPUT);
		for (final File inputFile: INPUT) {
            IOUtil.assertFileIsReadable(inputFile);
        }
		
        processOnlyPrimary();
		
		return 0;
	}
	
	
	public void processOnlyPrimary () {
        final SamHeaderAndIterator inputs = openInputs();
		CloseableIterator<SAMRecord> inputSam = inputs.iterator;
		SAMFileHeader header = inputs.header;
		header.addComment("Edit distance collapsed tag " +  this.PRIMARY_BARCODE + " to new tag " + this.OUT_BARCODE+ " with edit distance "+ this.EDIT_DISTANCE);
        SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, this.OUTPUT);
        
		// gather up the barcodes that exist in the BAM
        final SamHeaderAndIterator inputs2 = openInputs();
		ObjectCounter<String> barcodes = new BAMTagHistogram().getBamTagCounts(inputs2.iterator, this.PRIMARY_BARCODE,this.READ_QUALITY, this.FILTER_PCR_DUPLICATES);
        CloserUtil.close(inputs2.iterator);
        
		// filter barcodes by #reds in each barcode.
		barcodes=filterBarcodesByNumReads(barcodes, this.MIN_NUM_READS_NONCORE);
		
		// collapse them
		Map<String, String> childParentBarcodes=collapseBarcodes(this.MIN_NUM_READS_CORE, this.NUM_CORE_BARCODES, barcodes, this.FIND_INDELS, this.EDIT_DISTANCE);
		// iterate through the reads and retag with the proper reads.
		// log.info("STUFF");
		retagReads(inputSam, writer, childParentBarcodes, this.PRIMARY_BARCODE, this.OUT_BARCODE);
		// collapsed.size();
		
		CloserUtil.close(inputSam);
		writer.close();
	}

    private SamHeaderAndIterator openInputs() {
        final SamHeaderAndIterator ret = SamFileMergeUtil.mergeInputs(INPUT, true);
        if (SAMFileHeader.SortOrder.coordinate != ret.header.getSortOrder()) {
            throw new PicardException("Input files are not coordinate sorted");
        }
        return ret;
    }
	
	private ObjectCounter<String> filterBarcodesByNumReads (ObjectCounter<String> barcodes, int minNumReads) {
	
		ObjectCounter<String> result = new ObjectCounter<String>();
		for (String k: barcodes.getKeys()) {
			int count = barcodes.getCountForKey(k);
			if (count>=minNumReads) {
				result.setCount(k, count);
			}
		}
		log.info("Filtering barcodes by min non-core reads.  Started with [" + barcodes.getSize()+ "] ended with ["+ result.getSize()+"]");
		return (result);
	}
	
	private void retagReads (Iterator<SAMRecord> inputSam, SAMFileWriter writer, Map<String, String> collapsed, String tag, String outTag) {
		for (final SAMRecord r : new IterableAdapter<>(inputSam)) {
			pl.record(r);
			String s1 = r.getStringAttribute(tag);
			String newTag=collapsed.get(s1);
			if (newTag==null) newTag=s1;
			r.setAttribute(outTag, newTag);
			writer.addAlignment(r);
		}
	}
	
	private void retagReads (List<SAMRecord> batch, SAMFileWriter writer, Map<String, String> collapsed, String tag, String outTag) {		
		for (final SAMRecord r : batch) {
			pl.record(r);
			String s1 = r.getStringAttribute(tag);
			String newTag=collapsed.get(s1);
			
			if (newTag==null) newTag=s1;
			
			r.setAttribute(outTag, newTag);
			writer.addAlignment(r);
			pl.record(r);
		}
	}
	
	/**
	 * Collapse barcodes.  Does this for all barcodes.
	 * @param barcodes
	 * @return A map of each child barcode to it's parent.  Many keys will point to the same value.
	 */
	
	private Map<String, String> collapseBarcodes(ObjectCounter<String> barcodes, boolean findIndels, int editDistance) {
		List<String> barcodeList = barcodes.getKeysOrderedByCount(true);		
		Map<String, String> result = collapseBarcodes(barcodeList, barcodes, findIndels, editDistance);
		return (result);
	}
	
	
	/**
	 * Convenience method
	 * @param numReadsCore
	 * @param barcodes
	 * @param findIndels
	 * @param editDistance
	 * @return
	 */
	private Map<String, String> collapseBarcodes(Integer numReadsCore, Integer numCells, ObjectCounter<String> barcodes, boolean findIndels, int editDistance) {
		if (numReadsCore==null && numCells==null) return (collapseBarcodes(barcodes, findIndels, editDistance));
		// otherwise, select core barcodes and run.
		List<String> core=null;
		BarcodeListRetrieval u = new BarcodeListRetrieval();
		
		if (numReadsCore!=null) {
			core = u.getCoreBarcodesByReadCount(barcodes, numReadsCore);
		} 
		else if (numCells!=null) {
			core = u.getTopCoreBarcodesByReadCount (barcodes, numCells);
		}
		
		return (collapseBarcodes(core, barcodes, findIndels, editDistance));
	}
	
	
	public Map<String, String> collapseBarcodes(List<String> coreBarcodes, ObjectCounter<String> barcodes, boolean findIndels, int editDistance) {
		Map<String, String> result = new HashMap<String, String>();
		
		MapBarcodesByEditDistance med = new MapBarcodesByEditDistance(true, this.NUM_THREADS, 10000);
		Map<String, List<String>> r = med.collapseBarcodes(coreBarcodes, barcodes, findIndels, editDistance);
		for (String key: r.keySet()) {
			for (String value: r.get(key)) {
				result.put(value, key);
			}
		}
		return (result);
		
	}
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CollapseBarcodesInPlace().instanceMain(args));
	}
	
	
	/**
	 * Collapses a core set of barcodes.
	 * @param coreBarcodes
	 * @param barcodes
	 * @param findIndels
	 * @param editDistance
	 * @return
	 */
	/*
	public Map<String, String> collapseBarcodesOld(List<String> coreBarcodes, ObjectCounter<String> barcodes, boolean findIndels, int editDistance) {
		// don't allow side effects to modify input lists.
		coreBarcodes = new ArrayList<String>(coreBarcodes);
		barcodes = new ObjectCounter<String>(barcodes);
		
		Map<String, String> result = new HashMap<String, String>();
		int count = 0;
		int numBCCollapsed=0;
		List<BarcodeWithCount> barcodesWithCount=getBarcodesWithCounts(barcodes);
		List<String> barcodeList = EDUtils.getInstance().getBarcodes(barcodesWithCount);
		//int totalCount = barcodes.getTotalCount();
		
		int coreBarcodeCount=coreBarcodes.size();
		long startTime = System.currentTimeMillis();
		while (coreBarcodes.isEmpty()==false) {
			String b = coreBarcodes.get(0);
			count++;
			coreBarcodes.remove(b);
			barcodeList.remove(b);
			
			Set<String> closeBC=processSingle(b, barcodeList, findIndels, editDistance);
			numBCCollapsed+=closeBC.size();
			barcodeList.removeAll(closeBC);
			coreBarcodes.removeAll(closeBC);	
			//log.info("End remove all");
			for (String c: closeBC) {
				result.put(c, b);
			}
			
			if (count % this.REPORT_PROGRESS_INTERVAL == 0) {
				if (barcodes.getSize()>10000) log.info("Processed [" + count + "] records, totals BC Space left [" + barcodeList.size() +"]", " # collapsed this set [" + numBCCollapsed+"]");
				numBCCollapsed=0;
			}
		}
		long endTime = System.currentTimeMillis();
		long duration = (endTime - startTime)/1000;
		log.info("Collapse with [" + this.NUM_THREADS +"] threads took [" + duration + "] seconds to process");
		log.info("Started with core barcodes [" +coreBarcodeCount+  "] ended with [" + count + "] num collapsed [" +  (coreBarcodeCount-count) +"]");
		
		return (result);
	}
	
	
	private List<BarcodeWithCount> getBarcodesWithCounts (ObjectCounter<String> barcodes) {
		List<BarcodeWithCount> result = new ArrayList<BarcodeWithCount>();
		List<String> keys = barcodes.getKeysOrderedByCount(true);
		for (String k: keys) {
			BarcodeWithCount b = new BarcodeWithCount(k, barcodes.getCountForKey(k));
			result.add(b);
		}
		return (result);
	}
	
	
	public Set<String> processSingle(String barcode, List<String> comparisonBarcodes, boolean findIndels, int editDistance) {
		Set<String> closeBarcodes =null;
		if (this.NUM_THREADS>1) {
			closeBarcodes=cbt.getStringsWithinEditDistanceWithIndel(barcode, comparisonBarcodes, editDistance, findIndels);
		} else {
			// single threaded mode for now.  Maybe remove this later?  Not sure if single threaded is slower, probably is.
			if (findIndels) {
				closeBarcodes = EDUtils.getInstance().getStringsWithinEditDistanceWithIndel(barcode,comparisonBarcodes, editDistance);
			} else {
				closeBarcodes = EDUtils.getInstance().getStringsWithinEditDistance(barcode,comparisonBarcodes, editDistance);
			}	
		}
		return (closeBarcodes);
	}
	*/	
	
	
	/**
	 * Get all the reads for 1 tag.  This assumes that whatever tag is retrieved by the call to next() on the iterator gets the correct tag.
	 * When this method finishes, it should have peeked the next tag, but not called next to retrieve it.
	 * @param iter A peekable iterator of reads in tag order.
	 * @return A list of records, or null if the iterator is out.
	 */
	/*
	private List<SAMRecord> getReadsForNextTag (PeekableIterator<SAMRecord> iter, String tag, int readQuality, boolean filterDuplicates) {
		if (iter.hasNext()==false) return (null);
		List<SAMRecord> result = new ArrayList<SAMRecord>();
		
		SAMRecord r = iter.peek();
		String currentTag = r.getStringAttribute(tag);
		
		while (iter.hasNext()) {
			r=iter.peek();
			String nextTag = r.getStringAttribute(tag);
			// no gene info, you're done with all the records, and return null to signal the iterator is finished.
			if (nextTag==null) {
				return (null);
			}
			
			if (!nextTag.equals(currentTag)) {
				break;
			}
			// this is the same gene as before, keep grabbing records.
			// grab this record for "real" so peek gets the next record that might be in the same gene.
			iter.next();
			this.pl.record(r);
			
			// if there's filtering on read quality or PCR dupes, then filter on supplemental reads as well.
			if (readQuality>0 || filterDuplicates==true) {
				
				if (!r.isSecondaryOrSupplementary()&& r.getMappingQuality()>=readQuality) {
					if (filterDuplicates && r.getDuplicateReadFlag()==false) {
						result.add(r);
					}
					
				}
			} else {
				result.add(r);
			}
		}
		// if you didn't get any reads from this gene tag, then try again.
		if (result.size()==0) {
			return (getReadsForNextTag(iter, tag, readQuality, filterDuplicates));
		} else {
			return (result);
		}
		
	}
	*/
	/*
	public class CoreBCContainer {
		private String barcode;
		private List <String> childCoreBarcodes;
		private Set <String> childNonCoreBarcodes;
		
		int nextChildBarcodeIndex;
		
		public CoreBCContainer(String barcode) {
			this.barcode = barcode;
			this.childCoreBarcodes = new ArrayList<String>();
			this.childNonCoreBarcodes = new HashSet<String>();
			this.nextChildBarcodeIndex=0;
		}
		
		
		public void addChildCoreBarcode (String barcode) {
			if (!childCoreBarcodes.contains(barcode)) {
				childCoreBarcodes.add(barcode);
			}
		}
		
		public void addChildCoreBarcodes (Set<String> barcodes) {
			for (String b: barcodes) {
				addChildCoreBarcode(b);
			}
		}
		
		public void adChildNonCoreBarcode (String barcode) {
			childNonCoreBarcodes.add(barcode);
		}
		
		public void adChildNonCoreBarcodes (Set<String> barcodes) {
			childNonCoreBarcodes.addAll(barcodes);
		}
		
		/**
		 * If this has not yet been called on the object, hands off the first childCoreBarcode and then subsequent barcodes.
		 * @return the next child barcode, or null if no barcodes are left.
		 */
	/*
		public String getNextChildCoreBarcode () {
			if (!this.hasNextChildCoreBarcode()) return (null);
			String r = this.childCoreBarcodes.get(this.nextChildBarcodeIndex);
			this.nextChildBarcodeIndex++;
			return (r);
			 
		}
		
		public boolean hasNextChildCoreBarcode() {
			int listSize=this.childCoreBarcodes.size();
			if (this.nextChildBarcodeIndex>=listSize) return (false);
			return (true);
		}
		
		/**
		 * For this core barcode, does it have a child barcode that overlaps any other core barcode's child coreBC?
		 * If so, add that barcode's children to this barcode and search them as well.
		 * 
		 * @param other The list of core barcodes that were merged by this operation
		 * @return
		 */
	/*
		public Set<String> findChildOverlaps (Collection<CoreBCContainer> other) {
			Set<String> result = new HashSet<String>();
			// make a copy that I can safely mess with.
			Set<CoreBCContainer> o = new HashSet<CoreBCContainer>(other);
			
			while (this.hasNextChildCoreBarcode()) {
				String child = this.getNextChildCoreBarcode();
				for (CoreBCContainer c : o) {
					// ignore if the "other" and this are the same.
					if (c.barcode.equals(this.barcode)) {
						continue;
					}
					Set<String> otherChildren = new HashSet<String>(c.childCoreBarcodes);
					boolean hasOverlap = otherChildren.contains(child);
					if (hasOverlap) {
						this.addChildCoreBarcodes(otherChildren);
						result.add(c.barcode);
					}					
				}
			}
			
			return (result);
		}
		
	}
	*/
	
}
