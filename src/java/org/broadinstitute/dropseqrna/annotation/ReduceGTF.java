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
package org.broadinstitute.dropseqrna.annotation;

import java.io.File;
import java.io.PrintStream;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.cmdline.MetaData;
import org.broadinstitute.dropseqrna.utils.DropSeqSamUtil;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

/**
 * GTF files are annoyingly complex with a poor definition of what data is in them.
 * So hey, why not write a file parser.
 * This program reduces the GTF file in to a simplier, easier to parse format, while simultaneously allowing for data to be filtered.
 * @author nemesh
 *
 */

@CommandLineProgramProperties(
        usage = "GTF files are annoyingly complex with a poor definition of what data is in them. So hey, why not write a file parser. This program reduces the GTF file in to a simplier, easier to parse format, while simultaneously allowing for data to be filtered.",
        usageShort = "Parse and simplify a GTF file into an easier to use format.",
        programGroup = MetaData.class
)
public class ReduceGTF extends CommandLineProgram {

    private static final Log log = Log.getInstance(ReduceGTF.class);
    private static final List<String> DEFAULT_FEATURE_TYPES = CollectionUtil.makeList("gene", "transcript", "exon");

    private static final List<String> DEFAULT_IGNORED_FUNC_TYPES = CollectionUtil.makeList(
            "pseudogene", "polymorphic_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "IG_C_pseudogene",
            "IG_J_pseudogene", "IG_V_pseudogene");
    private static final String NA = "NA";

    @Option(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, doc="The reference sequence dictionary." +
            "  Only chromosomes found in this file AND the GTF file will be retained.")
	public File SEQUENCE_DICTIONARY;

	@Option(doc="The GTF file to reduce")
	public File GTF;

	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,doc="The output reduced GTF file.")
	public File OUTPUT;

	@Option(doc="Feature type(s) to extract. Only lines of the GTF that have these feature types will be extracted.  " +
            "This is the 3rd field of the GTF file, some examples of standard feature types are CDS, start_codon, stop_codon, and exon. ")
	public List<String> FEATURE_TYPE = DEFAULT_FEATURE_TYPES;

	@Option(doc="Functional type(s) to ignore.  These are values in the FUNCTIONAL_FIELD column in the GTF file.")
	public List<String> IGNORE_FUNC_TYPE = DEFAULT_IGNORED_FUNC_TYPES;

	@Option(doc="Enhance this reduced GTF file with genes,transcripts,introns, and consensus introns.  This is real " +
            "handy when your GTF file only defines exons, but has the transcript and gene IDs they belong to.")
	public boolean ENHANCE_GTF=true;

    private SAMSequenceDictionary dict;
    private boolean initialized = false;

	@Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
        IOUtil.assertFileIsReadable(GTF);
        IOUtil.assertFileIsWritable(this.OUTPUT);
        initialize();

        FilteringGTFParser parser = parseGTF();

        try {
            PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));
            writeHeader(out);

            // if no enhancement needed, just write out the results.
            if (!ENHANCE_GTF)
				writeRecords (out, parser);
			else {

                EnhanceGTFRecords e = new EnhanceGTFRecords();
                List<GTFRecord> records = e.enhanceGTFRecords(parser);
                Collections.sort(records, new GenomicOrderComparator(dict));
                writeRecords(out, records);
            }
            out.close();
            return 0;
        } finally {
            CloserUtil.close(parser);
        }
    }

    private void initialize() {
        if (!initialized) {
            dict = DropSeqSamUtil.loadSequenceDictionary(SEQUENCE_DICTIONARY);
            initialized = true;
        }
    }

    FilteringGTFParser parseGTF() {
        // Need to initialize for testing.
        initialize();
        return new FilteringGTFParser(this.GTF);
    }


	private class FilteringGTFParser extends FilteredIterator<GTFRecord> {
        final Set<String> featureTypes = new HashSet<>(FEATURE_TYPE);
        final Set<String> ignoredFunctionalTypes = new HashSet<>(IGNORE_FUNC_TYPE);

        private FilteringGTFParser(final File gtf) {
            super(new GTFParser(gtf, ValidationStringency.STRICT));
        }

        @Override
        public boolean filterOut(final GTFRecord rec) {
            return ignoredFunctionalTypes.contains(rec.getTranscriptType()) ||
                    !featureTypes.contains(rec.getFeatureType()) ||
                    dict.getSequence(rec.getChromosome()) == null;
        }
    }


	private void writeHeader (final PrintStream out) {
		String [] line = {"chr", "start", "end", "strand", "gene_name", "gene_id", "transcript_name", "transcript_id",
                "transcriptType", "annotationType"};
		String h = StringUtils.join(line, "\t");
		out.println(h);
	}

	private void writeRecords (final PrintStream out, final Iterable<GTFRecord> records) {
		for (GTFRecord r: records)
			writeLine(r, out);
	}

	private void writeLine (final GTFRecord r, final PrintStream out) {
		if (r==null) return;
		String [] line={r.getChromosome(),new Integer(r.getStart()).toString(), new Integer(r.getEnd()).toString(), r.getStrandAsString(), r.getGeneName(), r.getGeneID(),
				r.getTranscriptName(), r.getTranscriptID(), r.getTranscriptType(), r.getFeatureType()};
        for (int i = 0; i < line.length; ++i)
			if (line[i] == null || line[i].isEmpty())
				line[i] = NA;
		String h = StringUtils.join(line, "\t");
		out.println(h);
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new ReduceGTF().instanceMain(args));
	}


}
