package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

/**
 * Filter a BAM file by various filters and output a new filtered BAM.
 * Mostly for use with BAM tags and the like to avoid making tons of temp files and samtools messiness.
 * @author nemesh
 *
 */
@CommandLineProgramProperties(usage = "Filters a BAM file by various qualities to produce a new subset of the BAM containing the reads of interest.",
        usageShort = "Filters a BAM file by various qualities to produce a new subset of the BAM containing the reads of interest.",
        programGroup = DropSeq.class)
public class FilterBAM extends CommandLineProgram{
	private final Log log = Log.getInstance(FilterBAM.class);

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.")
	public File INPUT;

	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output report")
	public File OUTPUT;
	
	@Option(doc = "Minimum mapping quality to consider the read", optional=true)
	public Integer MINIMUM_MAPPING_QUALITY = null;
	
	@Option(doc = "Should PCR duplicates be filtered?", optional=true)
	public boolean FILTER_PCR_DUPES=false;
	
	@Option (doc = "Retain primary reads only", optional=true)
	public boolean RETAIN_ONLY_PRIMARY_READS=false;
	
	@Option (doc="Retain reads that have at least this many M bases total in the cigar string.  This sums all the M's in the cigar string.", optional=true)
	public Integer SUM_MATCHING_BASES=null;
	
	@Option (doc="Soft match reference names that have this string. If multiple matches are specificed, they are OR'd together." +
			"This is the equivilent of a hard match with wrapped with .* on either side.", optional=true)
	public List<String> REF_SOFT_MATCHED_RETAINED=null;
	
	@Option (doc="Soft match and reject reference names that have this string.  If multiple matches are specified, they are OR'd together. " +
			"This is the equivilent of a hard match with wrapped with .* on either side. ", optional=true)
	public List<String> REF_SOFT_MATCHED_REJECTED=null;
		
	@Option (doc="Exact match reference names that have this string.  If multiple matches are specified, they are OR'd together." +
			"For example, '1' would retain only references that were exactly '1'. This method accepts regular expressions.", optional=true)
	public List<String> REF_HARD_MATCHED_RETAINED=null;
	
	@Option (doc="Exact match and reject reference names that have this string.  If multiple matches are specified, they are OR'd together." +
			"This method accepts regular expressions.", optional=true)
	public List<String> REF_HARD_MATCHED_REJECTED=null;

	@Option (doc="Retain reads that have these tags set with any value.  Can be set multiple times", optional=true)
	public List<String> TAG_RETAIN=null;
	
	@Option(doc="If multiple TAG_RETAIN flags are set, should the result be the union of the filters, or the intersect?  [UNION/INTERSECT].", optional=true)
	public String TAG_RETAIN_COMBINE_FLAG=null;
	
	@Option (doc="Reject reads that have these tags set with any value.  Can be set multiple times.", optional=true)
	public List<String> TAG_REJECT=null;
	
	@Option(doc="If multiple TAG_REJECT flags are set, should the result be the union of the filters, or the intersect?  [UNION/INTERSECT].", optional=true)
	public String TAG_REJECT_COMBINE_FLAG=null;
	
	//@Option (doc="File with one or more TAG:Value combinations, for example ZC:Z:AAACCCTTGGG.  Any read with any of the tags in the file will be retained.")
	
	private String UNION="UNION";
	private String INTERSECT="INTERSECT";
	
	private Map<MatchTypes, List<Pattern>> patterns;
	
	public enum MatchTypes {
		REF_SOFT_MATCHED_RETAINED, REF_SOFT_MATCHED_REJECTED, REF_HARD_MATCHED_RETAINED, REF_HARD_MATCHED_REJECTED
	}
	
	public FilterBAM() {
	}
	
	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		buildPatterns();
		
		SamReader in = SamReaderFactory.makeDefault().open(INPUT);
		
        SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(), true, OUTPUT);
		ProgressLogger progLog=new ProgressLogger(log);
		
		
		for (final SAMRecord r : in) {
			progLog.record(r);
			if (filterRead (r)==false) {
				out.addAlignment(r);
			}
		}
		CloserUtil.close(in);
		out.close();
		return (0);
	}
	
	void buildPatterns () {
		this.patterns = new HashMap<MatchTypes, List<Pattern>>();
		if (REF_SOFT_MATCHED_RETAINED!=null) {
			List<Pattern> p = buildPatterns(this.REF_SOFT_MATCHED_RETAINED, true);
			this.patterns.put(MatchTypes.REF_SOFT_MATCHED_RETAINED, p);
		}
		if (REF_SOFT_MATCHED_REJECTED!=null) {
			List<Pattern> p = buildPatterns(this.REF_SOFT_MATCHED_REJECTED, true);
			this.patterns.put(MatchTypes.REF_SOFT_MATCHED_REJECTED, p);
		}
		if (REF_HARD_MATCHED_RETAINED!=null) {
			List<Pattern> p = buildPatterns(this.REF_HARD_MATCHED_RETAINED, false);
			this.patterns.put(MatchTypes.REF_HARD_MATCHED_RETAINED, p);
		}
		if (REF_HARD_MATCHED_REJECTED!=null) {
			List<Pattern> p = buildPatterns(this.REF_HARD_MATCHED_REJECTED, false);
			this.patterns.put(MatchTypes.REF_HARD_MATCHED_REJECTED, p);
		}
	}
	
	private List<Pattern> buildPatterns (List<String> matches, boolean soft) {
		List<Pattern> result = new ArrayList<Pattern>();
		for (String s: matches) {
			Pattern pattern=null;
			if (soft) {
				pattern=Pattern.compile(".*"+s+".*");
			} else {
				pattern=Pattern.compile(s);
			}
			
			result.add(pattern);
		}
		return (result);
	}
	
	
	public boolean filterRead (SAMRecord r) {
		// process all the reject filters first.
		if (rejectOnMapQuality(r)) return (true);
		if (rejectPCRDuplicate(r)) return (true);
		if (rejectNonPrimaryReads(r)) return (true);
		if (rejectOnCigar(r)) return (true);
		if (rejectSoftMatch(r)) return(true);
		if (rejectHardMatch(r)) return(true);
		if (rejectOnTags(this.TAG_REJECT, r)) return (true);

		// process all the accept filters.  All must pass.
		boolean sm=acceptSoftMatch(r);
		boolean hm=acceptHardMatch(r);
		boolean at=acceptOnTags(this.TAG_RETAIN, r);
		boolean all = sm && hm && at;
		return (!all);
		
	}
	
	/**
	 * Rejects reads if the sum of the cigar string bases is less than M_BASES_IN_CIGAR, reject the read.
	 * Don't process if M_BASES_IN_CIGAR is -1.
	 * @param r
	 * @return return false if the sum of the matching bases in the cigar is greater than the threshold.
	 */
	boolean rejectOnCigar(SAMRecord r) {		
		if (this.SUM_MATCHING_BASES==null) return (false);
		Cigar c = r.getCigar();
		int count=0;
		for (CigarElement ce: c.getCigarElements()) {
			if (ce.getOperator()==CigarOperator.M) {
				count+=ce.getLength();
			}
		}
		if (count>=this.SUM_MATCHING_BASES) {
			return false;
		}
		return true;
	}
	
	boolean rejectOnMapQuality (SAMRecord r) {
		if (this.MINIMUM_MAPPING_QUALITY==null) return (false);
		if (r.getMappingQuality() < this.MINIMUM_MAPPING_QUALITY) return (true);
		return (false);
	}
	
	boolean rejectPCRDuplicate (SAMRecord r) {
		if (this.FILTER_PCR_DUPES && r.getDuplicateReadFlag()) return (true);
		return (false);
	}
	
	boolean rejectNonPrimaryReads (SAMRecord r) {
		if (this.RETAIN_ONLY_PRIMARY_READS && r.isSecondaryOrSupplementary()) return (true);
		return (false);
	}
	
	boolean rejectSoftMatch(SAMRecord r) {
		if (this.REF_SOFT_MATCHED_REJECTED==null || this.REF_SOFT_MATCHED_REJECTED.isEmpty()) return false;
		boolean hasMatch = matchReference(this.patterns.get(MatchTypes.REF_SOFT_MATCHED_REJECTED), r);
		//boolean hasMatch = softMatchReference(this.REF_SOFT_MATCHED_REJECTED, r);
		return (hasMatch);
	}
	
	boolean rejectHardMatch(SAMRecord r) {
		if (this.REF_HARD_MATCHED_REJECTED==null || this.REF_HARD_MATCHED_REJECTED.isEmpty()) return false;
		//boolean hasMatch = exactMatchReference(this.REF_HARD_MATCHED_REJECTED, r);
		boolean hasMatch = matchReference(this.patterns.get(MatchTypes.REF_HARD_MATCHED_REJECTED), r);
		return (hasMatch);
	}
	
	boolean acceptSoftMatch(SAMRecord r) {
		if (this.REF_SOFT_MATCHED_RETAINED==null || this.REF_SOFT_MATCHED_RETAINED.isEmpty()) return true;
		//boolean hasMatch = softMatchReference(this.REF_SOFT_MATCHED_RETAINED, r);
		boolean hasMatch = matchReference(this.patterns.get(MatchTypes.REF_SOFT_MATCHED_RETAINED), r);
		return (hasMatch);
	}
	
	boolean acceptHardMatch(SAMRecord r) {
		if (this.REF_HARD_MATCHED_RETAINED==null || this.REF_HARD_MATCHED_RETAINED.isEmpty()) return true;
		//boolean hasMatch = exactMatchReference(this.REF_HARD_MATCHED_RETAINED, r);
		boolean hasMatch = matchReference(this.patterns.get(MatchTypes.REF_HARD_MATCHED_RETAINED), r);
		return (hasMatch);
	}
	
	
	/**
	 * Does the reference have a soft match for any of the proposed matches?  
	 * @param matches the list of possible matches.  All matches are wrapped with .* on either side.
	 * @param r the read
	 * @return return true if the record matches any of the filters.
	 */
	boolean softMatchReference (List<String> matches, SAMRecord r) {
		String refName = r.getReferenceName();
		for (String match: matches) {
			boolean m = refName.matches(".*"+match+".*");
			if (m) return (m);
		}
		return (false);
		
	}
	
	/**
	 * For matching, using patterns instead of string regex.
	 * @param patterns
	 * @param r
	 * @return
	 */
	boolean matchReference(List<Pattern> patterns, SAMRecord r) {
		String refName = r.getReferenceName();
		for (Pattern p: patterns) {
			Matcher m = p.matcher(refName);
			if (m.find()) {
				return true;
			}
		}
		return (false);
	}
	
	
	
	/**
	 * Does the reference have a soft match for any of the proposed matches?  
	 * @param matches the list of possible matches
	 * @param r the read
	 * @return return true if the record matches any of the filters.
	 */
	boolean exactMatchReference (List<String> matches, SAMRecord r) {
		String refName = r.getReferenceName();
		for (String match: matches) {
			boolean m = refName.matches(match);
			if (m) return (m);
		}
		return (false);
		
	}
	
	/**
	 * Reject on any of the tags if they are set.
	 * @param tags
	 * @param r
	 * @return
	 */
	boolean rejectOnTags (List<String> tags, SAMRecord r) {
		if (tags==null || tags.isEmpty()) return(false);
	
		for (String tag: tags) {
			Object v = r.getAttribute(tag);
			// if not null and union, return true if any are true.
			if (v!=null && this.TAG_REJECT_COMBINE_FLAG!=null && this.TAG_REJECT_COMBINE_FLAG.equals(this.UNION)) return true;
			// if null and intersect, return false
			if (v==null && this.TAG_REJECT_COMBINE_FLAG!=null && this.TAG_REJECT_COMBINE_FLAG.equals(this.INTERSECT)) return false;
			// if not null and single iteration, return true
			if (v!=null && this.TAG_REJECT_COMBINE_FLAG==null) return true;
		}
		return (false);
	}
	
	/**
	 * 
	 * @param tags
	 * @param r
	 * @return
	 */
	boolean acceptOnTags (List<String> tags, SAMRecord r) {
		if (tags==null || tags.isEmpty()) return(true);
		for (String tag: tags) {
			Object v = r.getAttribute(tag);
			// if not null and union, return true if any are true.
			if (v!=null && this.TAG_RETAIN_COMBINE_FLAG!=null && this.TAG_RETAIN_COMBINE_FLAG.equals(this.UNION)) return true;
			// if null and intersect, return false
			if (v==null && this.TAG_RETAIN_COMBINE_FLAG!=null && this.TAG_RETAIN_COMBINE_FLAG.equals(this.INTERSECT)) return false;
			// if not null and single iteration, return true
			if (v!=null && this.TAG_RETAIN_COMBINE_FLAG==null) return true;
		}
		return (false);
	}
	
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new FilterBAM().instanceMain(args));
	}
}
