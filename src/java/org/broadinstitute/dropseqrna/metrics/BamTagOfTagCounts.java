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
package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.io.PrintStream;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.CustomBAMIterators;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        summary = "For a given BAM tag, how many unique values of a second BAM tag are present?",
        oneLineSummary = "For a given BAM tag, how many unique values of a second BAM tag are present?",
        programGroup = DropSeq.class
)

public class BamTagOfTagCounts extends CommandLineProgram {

private static final Log log = Log.getInstance(BamTagOfTagCounts.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.  Must be coordinate sorted. (???)")
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file of tag frequencies. This supports zipped formats like gz and bz2.")
	public File OUTPUT;

	@Argument(doc="Primary Tag to extract")
	public String PRIMARY_TAG;

	@Argument(doc="Secondary Tag to extract")
	public String SECONDARY_TAG;

	@Argument(doc="Remove Singleton Results")
	public Boolean REMOVE_SINGLETONS=false;

	@Argument(doc="Filter PCR Duplicates.  Defaults to true")
	public boolean FILTER_PCR_DUPLICATES=true;

	@Argument(doc="Read quality filter.  Filters all reads lower than this mapping quality.  Defaults to 10.  Set to 0 to not filter reads by map quality.")
	public Integer READ_QUALITY=10;

	@Argument(doc="If the secondary tag can occur multiple times, break it up with this delimiter.", optional=true)
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

	public TagOfTagResults<String,String> getResults (final File inputBAM, final String primaryTag, final String secondaryTag, final boolean filterPCRDuplicates, final Integer readQuality) {
		TagOfTagResults<String,String> result = new TagOfTagResults<>();
		SamReader reader = SamReaderFactory.makeDefault().open(inputBAM);
		CloseableIterator<SAMRecord> iter = CustomBAMIterators.getReadsInTagOrder(reader, primaryTag);
		CloserUtil.close(reader);
		String currentTag="";
		Set<String> otherTagCollection=new HashSet<>();

		ProgressLogger progress = new ProgressLogger(log);

		while (iter.hasNext()) {
			SAMRecord r = iter.next();
			progress.record(r);
			// skip reads that don't pass filters.
			boolean discardResult=(filterPCRDuplicates && r.getDuplicateReadFlag()) || r.getMappingQuality()<readQuality || r.isSecondaryOrSupplementary();
			// short circuit if read is below the quality to be considered.
			if (discardResult) continue;

			Object d = r.getAttribute(secondaryTag);
			// short circuit if there's no tag for this read.
			if (d==null) continue;

			String data=null;

			if (d instanceof String)
				data=r.getStringAttribute(secondaryTag);
			else if (d instanceof Integer)
				data=Integer.toString(r.getIntegerAttribute(secondaryTag));


			String newTag = r.getStringAttribute(primaryTag);
			if (newTag==null)
				newTag="";

			// if you see a new tag.
			if (!currentTag.equals(newTag)) {
				// write out tag results, if any.
				if (!currentTag.equals("") && otherTagCollection.size()>0)
					result.addEntries(currentTag, otherTagCollection);
				currentTag=newTag;
				otherTagCollection.clear();
				if (!discardResult) otherTagCollection=addTagToCollection (data,otherTagCollection);

			} else // gather stats
			if (!discardResult) otherTagCollection=addTagToCollection (data,otherTagCollection);
		}
		if (otherTagCollection.size()>0)
			result.addEntries(currentTag, otherTagCollection);

		return (result);
	}

	private Set<String> addTagToCollection (final String data, final Set<String> collection) {
		if (SECONDARY_DELIMITER!=null) {
			String [] d2= data.split(this.SECONDARY_DELIMITER);
            Collections.addAll(collection, d2);
		} else
			collection.add(data);
		return (collection);
	}

	private void writeStats (final String tag, final Set<String> otherTagCollection, final PrintStream out) {
		if (!REMOVE_SINGLETONS || otherTagCollection.size()>1) {
			String otherTagList=StringUtils.join(otherTagCollection, ":");
			String [] line ={tag, otherTagCollection.size()+"", otherTagList};
			String h = StringUtils.join(line, "\t");
			out.println(h);
		}

	}

	private void writeHeader(final PrintStream out) {
		String [] header = {"TAG", "COUNT", "TAG_LIST"};
		String h = StringUtils.join(header, "\t");
		out.println(h);
	}


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new BamTagOfTagCounts().instanceMain(args));
	}


}
