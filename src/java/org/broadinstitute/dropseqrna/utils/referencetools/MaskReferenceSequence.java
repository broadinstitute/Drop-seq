/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils.referencetools;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.utils.FastaSequenceFileWriter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReferenceProgramGroup;

@CommandLineProgramProperties(
        summary = "Change reference sequence fasta to ignore contigs or portions thereof.  The bases of these contigs are set to N.",
        oneLineSummary = "Modify reference sequence fasta contig sequence.",
        programGroup = ReferenceProgramGroup.class
)
public class MaskReferenceSequence extends CommandLineProgram {

	private final Log log = Log.getInstance(MaskReferenceSequence.class);

	@Override
	protected boolean requiresReference() {
		return true;
	}

	@Argument(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output FASTA file to write.",  optional=false)
    public File OUTPUT;

	@Argument (doc="The number of bases per line in the output file")
	public Integer OUTPUT_LINE_LENGTH=50;

	@Argument (doc="A contig name to ignore, or a partial name that will be used as a pattern.  If the contig name contains any of the strings in this list it will be set to N.  This option can be used multiple times", mutex={"INTERVALS"})
	public List<String> CONTIG_PATTERN_TO_IGNORE;

	@Argument (doc="A file containing one or more intervals that will have their bases set to N. This file is in Interval format - tab seperated with fields: chr start end strand name\"", mutex={"CONTIG_PATTERN_TO_IGNORE"})
	public File INTERVALS;

	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(this.REFERENCE_SEQUENCE);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		// validate that an index is present for the reference sequence, since it's required.
		final ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE, true, true);
		if (!ref.isIndexed())
			throw new IllegalStateException ("Input fasta must be indexed.  You can do this by using samtools faidx to create an index");

		FastaSequenceFileWriter writer = new FastaSequenceFileWriter(OUTPUT, OUTPUT_LINE_LENGTH);
		if (this.CONTIG_PATTERN_TO_IGNORE!=null && !this.CONTIG_PATTERN_TO_IGNORE.isEmpty()) processByWholeContig(ref, writer, this.CONTIG_PATTERN_TO_IGNORE);
		if (this.INTERVALS!=null) processByPartialContig(ref, writer, this.INTERVALS);

		CloserUtil.close(ref);
		CloserUtil.close(writer);
		return 0;

	}

	private void processByPartialContig (final ReferenceSequenceFile ref, final FastaSequenceFileWriter writer, final File intervalListFile) {
		SAMSequenceDictionary sd = ref.getSequenceDictionary();
		// validate that the intervals and the reference have the same sequence dictionary.
		IntervalList iList = IntervalList.fromFile(intervalListFile);
		iList.getHeader().getSequenceDictionary().assertSameDictionary(sd);
		// map the intervals to a map to each contig.
		Map<String, List<Interval>> intervalsPerContig = getIntervalsForContig(iList);

		for (SAMSequenceRecord r: sd.getSequences()) {
			String contig = r.getSequenceName();
			log.info("Processing partial contig " + contig);
			// this list can be null.
			List<Interval> intervalsToMask = intervalsPerContig.get(contig);
			ReferenceSequence rs = ref.getSequence(contig);
			writeSequence(rs, intervalsToMask, writer);
		}
	}



	private Map<String, List<Interval>> getIntervalsForContig (final IntervalList iList) {
		Map<String,List<Interval>> result = new HashMap<>();
		for (Interval i: iList.getIntervals()) {
			String contig = i.getContig();
			List<Interval> r= result.get(contig);
			if (r==null) {
				r = new ArrayList<>();
				result.put(contig, r);
			}
			r.add(i);
		}
		return result;
	}


	private void processByWholeContig (final ReferenceSequenceFile ref, final FastaSequenceFileWriter writer, final List<String> contigPatternToIgnore) {
		SAMSequenceDictionary sd = ref.getSequenceDictionary();
		Set<String> contigsToIgnore = selectContigsToIgnore(sd, contigPatternToIgnore);
		// write out each contig.

		for (SAMSequenceRecord r: sd.getSequences()) {
			String contig = r.getSequenceName();
			log.info("Processing complete contig " + contig);
			ReferenceSequence rs = ref.getSequence(contig);
			boolean setSequenceToN = contigsToIgnore.contains(contig);
			writeSequence(rs, setSequenceToN, writer);
		}

	}

	private void writeSequence (final ReferenceSequence rs, final List<Interval> intervalsToMask, final FastaSequenceFileWriter writer) {
		String sequence = rs.getBaseString();
		if (intervalsToMask!=null && intervalsToMask.size()>0) {
			char [] seqArray = sequence.toCharArray();
			for (Interval i: intervalsToMask)
				for (int pos=i.getStart()-1; pos<i.getEnd(); pos++)
					seqArray[pos]='N';
			sequence=new String (seqArray);
		}
		writer.writeSequence(rs.getName(), sequence);
	}

	private void writeSequence (final ReferenceSequence rs, final boolean setSequenceToN, final FastaSequenceFileWriter writer) {
		String sequence = rs.getBaseString();
		if (setSequenceToN) sequence = StringUtils.repeat("N", sequence.length());
		writer.writeSequence(rs.getName(), sequence);
	}

	private Set<String> selectContigsToIgnore (final SAMSequenceDictionary sd, final List<String> patterns) {
		Set<String> result = new HashSet<>();
		for (SAMSequenceRecord r:  sd.getSequences()) {
			String contigName = r.getSequenceName();
			if (containsSubString(contigName, patterns))
				result.add(contigName);
		}
		return result;
	}

	private boolean containsSubString (final String contigName, final List<String> patterns) {
		for (String p: patterns)
			if (contigName.contains(p)) return true;
		return false;
	}


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new MaskReferenceSequence().instanceMain(args));
	}

}
