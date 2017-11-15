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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.util.*;
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
	
	@Option (doc="Soft match reference names that have this string. If multiple matches are specified, they are OR'd together." +
			"This is the equivalent of a hard match with wrapped with .* on either side.", optional=true)
	public List<String> REF_SOFT_MATCHED_RETAINED=null;
	
	@Option (doc="Soft match and reject reference names that have this string.  If multiple matches are specified, they are OR'd together. " +
			"This is the equivalent of a hard match with wrapped with .* on either side. ", optional=true)
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

	@Option(doc="Edit contig names so that a contig that starts with one of these prefixes has the prefix stripped.", optional = true)
    public List<String> STRIP_REF_PREFIX;

	@Option(doc="Edit sequence dictionary and remove any contig that has been filtered by reference name filtering. " +
            "A read with mate alignment info in which mate is aligned to a contig that has been removed will be changed " +
            "to have an unmapped mate.")
    public boolean DROP_REJECTED_REF = false;
	
	//@Option (doc="File with one or more TAG:Value combinations, for example ZC:Z:AAACCCTTGGG.  Any read with any of the tags in the file will be retained.")
	
	private static final String UNION="UNION";
	private static final String INTERSECT="INTERSECT";
	
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

		SAMFileHeader fileHeader = editSequenceDictionary(in.getFileHeader().clone());
		SamHeaderUtil.addPgRecord(fileHeader, this);
		SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(fileHeader, true, OUTPUT);
		ProgressLogger progLog=new ProgressLogger(log);

		final boolean sequencesRemoved = fileHeader.getSequenceDictionary().getSequences().size() != in.getFileHeader().getSequenceDictionary().getSequences().size();
		
		for (final SAMRecord r : in) {
			progLog.record(r);
			if (!filterRead(r)) {
                String sequenceName = stripReferencePrefix(r.getReferenceName());
                String mateSequenceName = null;
                if (r.getMateReferenceIndex() != -1) {
                    mateSequenceName = stripReferencePrefix(r.getMateReferenceName());
                }
			    if (sequencesRemoved || sequenceName != null) {
			        if (sequenceName == null) {
			            sequenceName = r.getReferenceName();
                    }
                    // Even if sequence name has not been edited, if sequences have been removed, need to set
                    // reference name again to invalidate reference index cache.
                    r.setReferenceName(sequenceName);
                }
                if (r.getMateReferenceIndex() != -1 && (sequencesRemoved || mateSequenceName != null)) {
                    if (mateSequenceName == null) {
                        mateSequenceName = r.getMateReferenceName();
                    }
                    // It's possible that the mate was mapped to a reference sequence that has been dropped from
                    // the sequence dictionary.  If so, set the mate to be unmapped.
                    if (fileHeader.getSequenceDictionary().getSequence(mateSequenceName) != null) {
                        r.setMateReferenceName(mateSequenceName);
                    } else {
                        r.setMateUnmappedFlag(true);
                        r.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                        r.setMateAlignmentStart(0);
                    }
                }
				out.addAlignment(r);
			}
		}
		CloserUtil.close(in);
		out.close();
		return (0);
	}

    private SAMFileHeader editSequenceDictionary(SAMFileHeader fileHeader) {
	    if (DROP_REJECTED_REF || !STRIP_REF_PREFIX.isEmpty()) {
	        // Make mutable copy of sequence list
            final ArrayList<SAMSequenceRecord> sequences = new ArrayList<>(fileHeader.getSequenceDictionary().getSequences());
            final ListIterator<SAMSequenceRecord> it = sequences.listIterator();
            while (it.hasNext()) {
                final SAMSequenceRecord sequence = it.next();
                if (DROP_REJECTED_REF && filterReference(sequence.getSequenceName())) {
                    it.remove();
                } else {
                    final String editedSequenceName = stripReferencePrefix(sequence.getSequenceName());
                    if (editedSequenceName != null) {
                        it.set(cloneWithNewName(sequence, editedSequenceName));
                    }
                }
            }
            fileHeader.getSequenceDictionary().setSequences(sequences);
        }
	    return fileHeader;
    }

    private SAMSequenceRecord cloneWithNewName(SAMSequenceRecord sequence, String editedSequenceName) {
        final SAMSequenceRecord ret = new SAMSequenceRecord(editedSequenceName, sequence.getSequenceLength());
        for (Map.Entry<String, String> entry : sequence.getAttributes()) {
            if (entry.getKey().equals(SAMSequenceRecord.SEQUENCE_NAME_TAG)) {
                ret.setAttribute(SAMSequenceRecord.SEQUENCE_NAME_TAG, editedSequenceName);
            } else {
                ret.setAttribute(entry.getKey(), entry.getValue());
            }
        }
        return ret;
    }

    /**
     * @param refName sequence name which may have a prefix stripped from it
     * @return edited refName, or null if no edit was done.
     */
    private String stripReferencePrefix(final String refName) {
        for (final String prefix : STRIP_REF_PREFIX) {
            if (refName.startsWith(prefix)) {
                return refName.substring(prefix.length());
            }
        }
        return null;
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
			final Pattern pattern;
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
		if (filterReference(r.getReferenceName())) return(true);
		if (rejectOnTags(this.TAG_REJECT, r)) return (true);

		if (!acceptOnTags(this.TAG_RETAIN, r)) return true;
		return false;
	}

	private boolean filterReference(final String refName) {
        return (rejectSoftMatch(refName) || rejectHardMatch(refName) || (!acceptSoftMatch(refName)) || (!acceptHardMatch(refName)));
    }
	
	/**
	 * Rejects reads if the sum of the cigar string bases is less than M_BASES_IN_CIGAR, reject the read.
	 * Don't process if M_BASES_IN_CIGAR is -1.
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

    private boolean rejectOnMapQuality (SAMRecord r) {
		if (this.MINIMUM_MAPPING_QUALITY==null) return (false);
		if (r.getMappingQuality() < this.MINIMUM_MAPPING_QUALITY) return (true);
		return (false);
	}

    private boolean rejectPCRDuplicate (SAMRecord r) {
		if (this.FILTER_PCR_DUPES && r.getDuplicateReadFlag()) return (true);
		return (false);
	}
	
	boolean rejectNonPrimaryReads (SAMRecord r) {
		if (this.RETAIN_ONLY_PRIMARY_READS && r.isSecondaryOrSupplementary()) return (true);
		return (false);
	}

    private boolean rejectSoftMatch(final String refName) {
		if (this.REF_SOFT_MATCHED_REJECTED==null || this.REF_SOFT_MATCHED_REJECTED.isEmpty()) return false;
		boolean hasMatch = matchReference(this.patterns.get(MatchTypes.REF_SOFT_MATCHED_REJECTED), refName);
		//boolean hasMatch = softMatchReference(this.REF_SOFT_MATCHED_REJECTED, r);
		return (hasMatch);
	}

    private boolean rejectHardMatch(final String refName) {
		if (this.REF_HARD_MATCHED_REJECTED==null || this.REF_HARD_MATCHED_REJECTED.isEmpty()) return false;
		//boolean hasMatch = exactMatchReference(this.REF_HARD_MATCHED_REJECTED, r);
		boolean hasMatch = matchReference(this.patterns.get(MatchTypes.REF_HARD_MATCHED_REJECTED), refName);
		return (hasMatch);
	}

    private boolean acceptSoftMatch(final String refName) {
		if (this.REF_SOFT_MATCHED_RETAINED==null || this.REF_SOFT_MATCHED_RETAINED.isEmpty()) return true;
		//boolean hasMatch = softMatchReference(this.REF_SOFT_MATCHED_RETAINED, r);
		boolean hasMatch = matchReference(this.patterns.get(MatchTypes.REF_SOFT_MATCHED_RETAINED), refName);
		return (hasMatch);
	}

    private boolean acceptHardMatch(final String refName) {
		if (this.REF_HARD_MATCHED_RETAINED==null || this.REF_HARD_MATCHED_RETAINED.isEmpty()) return true;
		//boolean hasMatch = exactMatchReference(this.REF_HARD_MATCHED_RETAINED, r);
		boolean hasMatch = matchReference(this.patterns.get(MatchTypes.REF_HARD_MATCHED_RETAINED), refName);
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
			if (refName.matches(".*"+match+".*")) {
                return true;
            }
		}
		return (false);
		
	}
	
	/**
	 * For matching, using patterns instead of string regex.
	 */
    private boolean matchReference(List<Pattern> patterns, final String refName) {
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
			if (refName.matches(match)) {
                return true;
            }
		}
		return (false);
		
	}
	
	/**
	 * Reject on any of the tags if they are set.
	 */
	boolean rejectOnTags (List<String> tags, SAMRecord r) {
		if (tags==null || tags.isEmpty()) return(false);
	
		for (String tag: tags) {
			Object v = r.getAttribute(tag);
			// if not null and union, return true if any are true.
			if (v!=null && this.TAG_REJECT_COMBINE_FLAG!=null && this.TAG_REJECT_COMBINE_FLAG.equals(UNION)) return true;
			// if null and intersect, return false
			if (v==null && this.TAG_REJECT_COMBINE_FLAG!=null && this.TAG_REJECT_COMBINE_FLAG.equals(INTERSECT)) return false;
			// if not null and single iteration, return true
			if (v!=null && this.TAG_REJECT_COMBINE_FLAG==null) return true;
		}
		return (false);
	}

    private boolean acceptOnTags (List<String> tags, SAMRecord r) {
		if (tags==null || tags.isEmpty()) return(true);
		for (String tag: tags) {
			Object v = r.getAttribute(tag);
			// if not null and union, return true if any are true.
			if (v!=null && this.TAG_RETAIN_COMBINE_FLAG!=null && this.TAG_RETAIN_COMBINE_FLAG.equals(UNION)) return true;
			// if null and intersect, return false
			if (v==null && this.TAG_RETAIN_COMBINE_FLAG!=null && this.TAG_RETAIN_COMBINE_FLAG.equals(INTERSECT)) return false;
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
