package org.broadinstitute.dropseqrna.metrics;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.CustomBAMIterators;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.PrintStream;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@CommandLineProgramProperties(
        usage = "For a given BAM tag, how many unique values of a second BAM tag are present?",
        usageShort = "For a given BAM tag, how many unique values of a second BAM tag are present?",
        programGroup = DropSeq.class
)

public class BAMTagofTagCounts extends CommandLineProgram {

private static final Log log = Log.getInstance(BAMTagofTagCounts.class);

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.  Must be coordinate sorted. (???)")
	public File INPUT;
	
	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file of tag frequencies. This supports zipped formats like gz and bz2.")
	public File OUTPUT;
	
	@Option(doc="Primary Tag to extract")
	public String PRIMARY_TAG;
	
	@Option(doc="Secondary Tag to extract")
	public String SECONDARY_TAG;
	
	@Option(doc="Remove Singleton Results")
	public Boolean REMOVE_SINGLETONS=false;
	
	@Option(doc="Filter PCR Duplicates.  Defaults to true")
	public boolean FILTER_PCR_DUPLICATES=true;
	
	@Option(doc="Read quality filter.  Filters all reads lower than this mapping quality.  Defaults to 10.  Set to 0 to not filter reads by map quality.")
	public Integer READ_QUALITY=10;
	
	@Option(doc="If the secondary tag can occur multiple times, break it up with this delimiter.", optional=true)
	public String SECONDARY_DELIMITER;
	
	public static final int MAX_RECORDS_IN_RAM = 500000;
	
	@Override
	protected int doWork() {
		
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));
		
		writeHeader(out);
		
		TagOfTagResults<String,String> results= getResults(this.INPUT, this.PRIMARY_TAG, this.SECONDARY_TAG, this.FILTER_PCR_DUPLICATES, this.READ_QUALITY);
		
		for (String k: results.getKeys()) {
			Set<String> values = results.getValues(k);
			writeStats(k, values, out);
		}
		out.close();
		
		return(0);
	}
	
	public TagOfTagResults<String,String> getResults (File inputBAM, String primaryTag, String secondaryTag, boolean filterPCRDuplicates, Integer readQuality) {
		TagOfTagResults<String,String> result = new TagOfTagResults<String,String>();
		SamReader reader = SamReaderFactory.makeDefault().open(inputBAM);	
		CloseableIterator<SAMRecord> iter = CustomBAMIterators.getReadsInTagOrder(reader, primaryTag);
		CloserUtil.close(reader);
		String currentTag="";
		Set<String> otherTagCollection=new HashSet<String>();
		
		ProgressLogger progress = new ProgressLogger(log);
		
		while (iter.hasNext()) {
			SAMRecord r = iter.next();
			progress.record(r);
			// skip reads that don't pass filters.
			
			boolean discardResult=(filterPCRDuplicates && r.getDuplicateReadFlag()) || r.getMappingQuality()<readQuality || r.isSecondaryOrSupplementary();
			String data=r.getStringAttribute(secondaryTag);
			
			// short circuit if there's no tag for this read.
			if (data==null) discardResult=true; 
			
			String newTag = r.getStringAttribute(primaryTag);
			if (newTag==null) {
				newTag="";
			}
					
			// if you see a new tag.
			if (!currentTag.equals(newTag)) {
				// write out tag results, if any.
				if (!currentTag.equals("") && otherTagCollection.size()>0) {
					result.addEntries(currentTag, otherTagCollection);	
				}
				currentTag=newTag;
				otherTagCollection.clear();
				if (!discardResult) otherTagCollection=addTagToCollection (data,otherTagCollection);
					
			} else {
				// gather stats
				if (!discardResult) otherTagCollection=addTagToCollection (data,otherTagCollection);
			}
		}
		if (otherTagCollection.size()>0) {
			result.addEntries(currentTag, otherTagCollection);
		}
		
		return (result);
	}
	// TODO: Need to make SECONDARY_DELIMITER be a passed in variable 
	private Set<String> addTagToCollection (String data, Set<String> collection) {
		if (SECONDARY_DELIMITER!=null) {
			String [] d2= data.split(this.SECONDARY_DELIMITER);
            Collections.addAll(collection, d2);
		} else {
			collection.add(data);
		}
		return (collection);
	}
	
	private void writeStats (String tag, Set<String> otherTagCollection, PrintStream out) {
		if (!REMOVE_SINGLETONS || otherTagCollection.size()>1) {
			String otherTagList=StringUtils.join(otherTagCollection, ":");
			String [] line ={tag, otherTagCollection.size()+"", otherTagList};
			String h = StringUtils.join(line, "\t");
			out.println(h);
		}
		
	}
	
	private void writeHeader(PrintStream out) {
		String [] header = {"TAG", "COUNT", "TAG_LIST"};
		String h = StringUtils.join(header, "\t");
		out.println(h);
	}
	
	public CloseableIterator<SAMRecord> getReadsInTagOrder (File bamFile) {
		SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
		SAMSequenceDictionary dict= reader.getFileHeader().getSequenceDictionary();
		List<SAMProgramRecord> programs =reader.getFileHeader().getProgramRecords();
		
		final SAMFileHeader writerHeader = new SAMFileHeader();
        writerHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        writerHeader.setSequenceDictionary(dict);
        for (SAMProgramRecord spr : programs) {
        	writerHeader.addProgramRecord(spr);
        }

        log.info("Reading in records for TAG name sorting");
        final ProgressLogger progressLogger = new ProgressLogger(log, 1000000);
        final CloseableIterator<SAMRecord> result =
                SamRecordSortingIteratorFactory.create(writerHeader, reader.iterator(), new StringTagComparator(this.PRIMARY_TAG), progressLogger);

		log.info("Sorting finished.");
		return (result);
	}
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new BAMTagofTagCounts().instanceMain(args));
	}
	
	
}
