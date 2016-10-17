package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;

import java.io.File;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.broadinstitute.dropseqrna.cmdline.MetaData;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.referencetools.ReferenceUtils;

import picard.annotation.Gene;
import picard.annotation.Gene.Transcript;
import picard.annotation.Gene.Transcript.Exon;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        usage = "Given a GTF file and a reference sequence, produce a report containing the %GC and length of each gene.  "
        		+ "GC is calculated for each gene by finding the unique set of base positions overlapping an exon and counting [G/C] bases compared to the total number of bases."
        		+ "Length is calculated by computing the length of each transcript for the gene, and taking the median value",
        usageShort = "Calculate GC content and length for genes",
        programGroup = MetaData.class
)

public class GatherGeneGCLength extends CommandLineProgram {

    private static final Log log = Log.getInstance(GatherGeneGCLength.class);

    @Option(doc="The GTF file containing gene models to generate length and GC metrics from")
	public File GTF;

	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,doc="The output report containg the genes and GC/Length metrics.")
	public File OUTPUT;

	@Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="The reference fasta")
    public File REFERENCE;

	// store a copy of this and re-use instead of constructing new ones.
	private DescriptiveStatistics stats = new DescriptiveStatistics();

	// format output percentages.
	private DecimalFormat percentageFormat = new DecimalFormat("###.#");

	@Override
    protected int doWork() {

        IOUtil.assertFileIsReadable(GTF);
        IOUtil.assertFileIsWritable(this.OUTPUT);
        PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));
        writeHeader(out);

        ReferenceSequenceFileWalker refFileWalker = new ReferenceSequenceFileWalker(REFERENCE);

        SAMSequenceDictionary dict= refFileWalker.getSequenceDictionary();

        OverlapDetector<Gene> geneOverlapDetector= GeneAnnotationReader.loadAnnotationsFile(GTF, dict);

        List<SAMSequenceRecord> records = dict.getSequences();

		for (SAMSequenceRecord record: records) {
			String seqName = record.getSequenceName();
			int seqIndex=dict.getSequenceIndex(seqName);
			ReferenceSequence fastaRef=refFileWalker.get(seqIndex);

			// get the genes for this contig.
			Interval i = new Interval(seqName, 1, record.getSequenceLength());
			Collection< Gene> genes = geneOverlapDetector.getOverlaps(i);
			for (Gene g: genes) {
				GCResult gc = calculateGCContent(g, fastaRef, dict);
				int medianLength = getMedianTranscriptLength(g);
				writeResult(g, gc, medianLength, out);
			}
		}
		CloserUtil.close(refFileWalker);
		CloserUtil.close(out);
        return 0;
	}



	/**
	 * For a GC record and a fasta sequence, calculate the GC content.
	 * Builds intervals of the unique sequences overlapped by exons, calculates the GC content for each, and aggregates results.
	 * @param gtfRecord
	 * @param fastaRef
	 * @return
	 */
	private GCResult calculateGCContent(final Gene gene, final ReferenceSequence fastaRef, final SAMSequenceDictionary dict) {
		// make an interval list.
		SAMFileHeader h = new SAMFileHeader();
		h.setSequenceDictionary(dict);
		h.setSortOrder(SAMFileHeader.SortOrder.unsorted);
		IntervalList intervalList  = new IntervalList(h);

		for (Transcript t : gene)
			for (Exon e: t.exons)
				intervalList.add(new Interval (gene.getContig(), e.start, e.end, gene.isNegativeStrand(), gene.getName()));

		List<Interval> uniqueIntervals = IntervalList.getUniqueIntervals(intervalList, false);

		// track aggregated GC.
		GCResult result = new GCResult(0, 0);

		for (Interval i: uniqueIntervals) {
			GCResult gcResultInterval = calculateGCContent(i, fastaRef);
			result.incrementGCCount(gcResultInterval.getGcCount());
			result.incrementRegionLength(gcResultInterval.getRegionLength());
		}
		return result;
	}

	/**
	 * Get the median transcript length fromt the gene model.
	 * @param gene
	 * @return
	 */
	private int getMedianTranscriptLength (final Gene gene) {
		List<Double> transcriptLengths= new ArrayList<Double>();

		for (Transcript t : gene) {
			double length = (t.transcriptionEnd - t.transcriptionStart) +1;
			transcriptLengths.add(length);
		}

		// calculate the median.
		Collections.sort(transcriptLengths);
		int numTranscripts=transcriptLengths.size();
		int middle = numTranscripts/2;
		if (numTranscripts%2==1)
			return (int) Math.round(transcriptLengths.get(middle));
		else {
			double result = (transcriptLengths.get(middle-1) + transcriptLengths.get(middle))/2;
			return (int) Math.round(result);
		}

	}


	private GCResult calculateGCContent (final Interval interval, final ReferenceSequence fastaRef) {
		//ensure that the start and end positions occur within the genome.
		String seq=ReferenceUtils.getSequence (fastaRef.getBases(), interval);
		seq=seq.toUpperCase();
		char [] seqArray=seq.toCharArray();
		int gcCount=0;
		int length=0;
		for (char c: seqArray) {
			length++;
			if (c=='G' || c=='C')
				gcCount++;
		}
		GCResult result = new GCResult(length, gcCount);

		return (result);
	}

	private class GCResult {
		private int regionLength=0;
		private int gcCount=0;

		public GCResult (final int regionLength, final int gcCount) {
			this.regionLength=regionLength;
			this.gcCount=gcCount;
		}

		public int getRegionLength() {
			return regionLength;
		}

		public int getGcCount() {
			return gcCount;
		}

		public double getGCPercent() {
			double result = ((double) this.gcCount / (double) this.regionLength)*100;
			return (result);
		}

		public void incrementRegionLength(final int length) {
			this.regionLength+=length;
		}

		public void incrementGCCount (final int count) {
			this.gcCount+=count;
		}

		@Override
		public String toString() {
			return "";
		}
	}

	 private void writeHeader(final PrintStream out) {
		 String [] header = {"GENE", "CHR", "START", "END", "PCT_GC_UNIQUE_EXON_BASES", "MEDIAN_TRANSCRIPT_LENGTH"};
		 String h = StringUtils.join(header, "\t");
		 out.println(h);
	 }

	private void writeResult(final Gene gene, final GCResult gc,
			final int medianTranscriptLength, final PrintStream out) {

		String[] line = {gene.getName(), gene.getContig(), Integer.toString(gene.getStart()), Integer.toString(gene.getEnd()), percentageFormat.format(gc.getGCPercent()), Integer.toString(medianTranscriptLength)};
		String h = StringUtils.join(line, "\t");
		out.println(h);
	}


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new GatherGeneGCLength().instanceMain(args));
	}

}
