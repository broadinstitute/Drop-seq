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
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.broadinstitute.dropseqrna.cmdline.MetaData;
import org.broadinstitute.dropseqrna.utils.FastaSequenceFileWriter;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.referencetools.ReferenceUtils;

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
import htsjdk.samtools.util.SequenceUtil;
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

	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,doc="The output report containg the genes and GC/Length metrics.  Output at the gene level, using the median values across transcripts.")
	public File OUTPUT;

	@Option(doc="The output report containg the genes and GC/Length metrics at the transcript level.", optional=true)
	public File OUTPUT_TRANSCRIPT_LEVEL;

	@Option(doc="The sequences of each transcript", optional=true)
	public File OUTPUT_TRANSCRIPT_SEQUENCES;

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

        PrintStream outTranscript = null;
        if (this.OUTPUT_TRANSCRIPT_LEVEL!=null) {
			outTranscript = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT_TRANSCRIPT_LEVEL));
			writeHeaderTranscript(outTranscript);
        }

        FastaSequenceFileWriter  outSequence = null;

        if (this.OUTPUT_TRANSCRIPT_SEQUENCES!=null) {
        	IOUtil.assertFileIsWritable(this.OUTPUT_TRANSCRIPT_SEQUENCES);
			outSequence = new FastaSequenceFileWriter (this.OUTPUT_TRANSCRIPT_SEQUENCES);
        }
        ReferenceSequenceFileWalker refFileWalker = new ReferenceSequenceFileWalker(REFERENCE);

        SAMSequenceDictionary dict= refFileWalker.getSequenceDictionary();
        if (dict==null) {
        	CloserUtil.close(refFileWalker);
        	throw new IllegalArgumentException("Reference file" + this.REFERENCE.getAbsolutePath()+" is missing a dictionary file [.dict].  Please make one!");
        }

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
				List<GCResult> gcList = calculateGCContentGene(g, fastaRef, dict);
				if (this.OUTPUT_TRANSCRIPT_LEVEL!=null)
					writeResultTranscript(gcList, outTranscript);
				GCIsoformSummary summary = new GCIsoformSummary(g, gcList);
				if (this.OUTPUT_TRANSCRIPT_SEQUENCES!=null)
					writeTranscriptSequence(g, fastaRef, dict, outSequence);

				GCResult gc = calculateGCContentUnionExons(g, fastaRef, dict);

				writeResult(gc, summary, out);
			}
		}
		CloserUtil.close(refFileWalker);
		CloserUtil.close(out);
		if (this.OUTPUT_TRANSCRIPT_LEVEL!=null) CloserUtil.close(outTranscript);
		if (this.OUTPUT_TRANSCRIPT_SEQUENCES!=null) CloserUtil.close(outSequence);
        return 0;
	}

	/**
	 * For a GC record and a fasta sequence, calculate the GC content.
	 * Builds intervals of the unique sequences overlapped by exons, calculates the GC content for each, and aggregates results.
	 * @param gtfRecord
	 * @param fastaRef
	 * @return
	 */
	private GCResult calculateGCContentUnionExons(final Gene gene, final ReferenceSequence fastaRef, final SAMSequenceDictionary dict) {
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
		GCResult result = new GCResult(0, 0, 0);

		for (Interval i: uniqueIntervals) {
			GCResult gcResultInterval = calculateGCContentExon(i, fastaRef);
			result.increment(gcResultInterval);
		}
		return result;
	}


	private List<GCResult> calculateGCContentGene (final Gene gene, final ReferenceSequence fastaRef, final SAMSequenceDictionary dict) {
		List<GCResult> result = new ArrayList<>();

		for (Transcript t : gene) {
			String seq=getTranscriptSequence(t, fastaRef, dict);
			GCResult gc = new GCResult(seq);
			gc.setTranscript(t);

			// check for GQuadruplexes.
			List<GQuadruplex> gq = GQuadruplex.find(t.name, seq);
			gc.incrementGQuadruplexCount(gq.size());
			result.add(gc);
		}
		return (result);
	}


	public static double getMedian (final double [] data) {
		Arrays.sort(data);
		int numTranscripts=data.length;
		int middle = numTranscripts/2;
		if (numTranscripts%2==1) {
			double result = data[middle];
			return (result);

		}
		else {
			double result = (data[middle-1] + data[middle])/2;
			return result;

		}
	}

	public void writeTranscriptSequence (final Gene gene, final ReferenceSequence fastaRef, final SAMSequenceDictionary dict, final FastaSequenceFileWriter outSequence ) {

		for (Transcript t : gene) {
			String seqName=gene.getName()+" " + t.name;
			String sequence=getTranscriptSequence(t, fastaRef, dict);
			outSequence.writeSequence(seqName, sequence);
		}
	}

	/**
	 * For a given transcript model, gather up the exons in genomic order, get their sequences, stitch them together,
	 * upper case them, then reverse compliment if the gene this transcript belongs to is on the negative strand.
	 * @param transcript
	 * @param fastaRef
	 * @param dict
	 * @return
	 */
	public String getTranscriptSequence (final Transcript transcript, final ReferenceSequence fastaRef, final SAMSequenceDictionary dict) {

		StringBuilder b = new StringBuilder();

		for (Exon e: transcript.exons) {
			Interval i= new Interval ( transcript.getGene().getContig(), e.start, e.end, transcript.getGene().isNegativeStrand(), transcript.getGene().getName());
			String seq=ReferenceUtils.getSequence (fastaRef.getBases(), i);
			b.append(seq);
		}
		// build the sequence in genomic order, upper case, reverse compliment if needed.
		String finalSeq = b.toString();
		finalSeq=finalSeq.toUpperCase();
		// if the sequence is on the negative strand, reverse compliemnt.
		if (transcript.getGene().isNegativeStrand()) finalSeq=SequenceUtil.reverseComplement(finalSeq);
		// log.info(">" + transcript.name + "\n" + finalSeq);
		return (finalSeq);
	}


	private GCResult calculateGCContentExon (final Interval interval, final ReferenceSequence fastaRef) {
		//ensure that the start and end positions occur within the genome.
		String seq=ReferenceUtils.getSequence (fastaRef.getBases(), interval);
		// if the sequence is on the negative strand, reverse compliemnt.
		if (interval.isNegativeStrand()) seq=SequenceUtil.reverseComplement(seq);
		GCResult result = new GCResult(seq);
		return (result);
	}


	public class GCIsoformSummary {

		private List<GCResult> transcriptGCList;
		private Gene gene;

		public GCIsoformSummary (final Gene gene, final List<GCResult> transcriptGCList) {
			this.transcriptGCList=transcriptGCList;
			this.gene=gene;
		}

		public double getMedianGC () {
			double result = getMedian(transcriptGCList.stream().mapToDouble(GCResult::getGCPercent).toArray());
			return (result);
		}

		public double getMedianG () {
			double result = getMedian(transcriptGCList.stream().mapToDouble(GCResult::getGPercent).toArray());
			return (result);
		}

		public double getMedianC () {
			double result = getMedian(transcriptGCList.stream().mapToDouble(GCResult::getCPercent).toArray());
			return (result);
		}

		public int getMedianTranscriptLength() {
			double result = getMedian(transcriptGCList.stream().mapToDouble(GCResult::getRegionLength).toArray());
			return (int) Math.round(result);
		}

		public int getMedianGQuadruplexes() {
			return (int) Math.round(getMedian(transcriptGCList.stream().mapToDouble(GCResult::getNumGQuadruplexesObserved).toArray()));
		}

		public int getNumTranscripts () {
			return this.transcriptGCList.size();
		}

		public Gene getGene () {
			return this.gene;
		}

		@Override
		public String toString () {
			StringBuilder b = new StringBuilder();
			b.append(this.gene.toString());
			b.append(" %GC [" + percentageFormat.format(this.getMedianGC())+"]");
			b.append(" %G [" + percentageFormat.format(this.getMedianG())+"]");
			b.append(" %C [" + percentageFormat.format(this.getMedianC())+"]");
			b.append(" median GQuadruplex [" + this.getMedianGQuadruplexes() +"]");
			b.append(" Length [" + this.getMedianTranscriptLength()+"]");

			return b.toString();
		}

	}

	public class GCResult {
		private int regionLength=0;
		private Transcript transcript=null;

		private int gCount=0;
		private int cCount=0;
		private int numGQuadruplexesObserved=0;

		public GCResult (final String sequence) {
			String seq=sequence.toUpperCase();
			char [] seqArray=seq.toCharArray();
			regionLength=seqArray.length;
			for (char c: seqArray) {
				if (c=='C') this.incrementC(1);
				if (c=='G') this.incrementG(1);
			}
		}

		public GCResult (final int regionLength, final int cCount, final int gCount) {
			this.regionLength=regionLength;
			this.cCount+=cCount;
			this.gCount+=gCount;
		}

		public int getRegionLength() {
			return regionLength;
		}

		public int getGcCount() {
			return gCount+cCount;
		}

		public double getGCPercent() {
			double result = (this.getGcCount() / (double) this.regionLength)*100;
			return (result);
		}

		public double getCPercent () {
			double result = ((double) this.cCount / (double) this.regionLength)*100;
			return (result);
		}

		public double getGPercent () {
			double result = ((double) this.gCount / (double) this.regionLength)*100;
			return (result);
		}

		public void incrementRegionLength(final int length) {
			this.regionLength+=length;
		}

		public void incrementG (final int count) {
			this.gCount+=count;
		}

		public void incrementC (final int count) {
			this.cCount+=count;
		}

		public void increment(final GCResult other) {
			this.cCount+=other.cCount;
			this.gCount+=other.gCount;
			this.regionLength+=other.regionLength;
		}

		public void incrementGQuadruplexCount(final int count) {
			this.numGQuadruplexesObserved+=count;
		}

		public int getNumGQuadruplexesObserved() {
			return this.numGQuadruplexesObserved;
		}

		public Transcript getTranscript() {
			return transcript;
		}

		public void setTranscript(final Transcript transcript) {
			this.transcript = transcript;
		}


		@Override
		public boolean equals(final Object obj) {
			   if (obj == null)
				return false;
			   if (obj == this)
				return true;
			   if (obj.getClass() != getClass())
				return false;
			   GCResult rhs = (GCResult) obj;
			   return new EqualsBuilder()
			                 .appendSuper(super.equals(obj))
			                 .append(gCount, rhs.gCount)
			                 .append(cCount, rhs.cCount)
			                 .append(regionLength, rhs.regionLength)
			                 .append(numGQuadruplexesObserved, rhs.numGQuadruplexesObserved)
			                 .isEquals();
		}

		@Override
		public String toString() {
			return "Length [" +this.regionLength +"] %G [" +percentageFormat.format(this.getGPercent()) +"] %C [" + percentageFormat.format(this.getCPercent()) +"] %GC ["+percentageFormat.format(this.getGCPercent()) +"]" + "G-Quadruplexes [" + this.numGQuadruplexesObserved +"]";
		}
	}



	 private void writeHeader(final PrintStream out) {
		 String [] header = {"GENE", "CHR", "START", "END", "PCT_GC_UNIQUE_EXON_BASES", "PCT_GC_ISOFORM_AVERAGE", "PCT_C_ISOFORM_AVERAGE", "PCT_G_ISOFORM_AVERAGE", "MEDIAN_TRANSCRIPT_LENGTH", "NUM_TRANSCRIPTS", "MEDIAN_GQUADRUPLEXES"};
		 String h = StringUtils.join(header, "\t");
		 out.println(h);
	 }

	private void writeResult(final GCResult gc, final GCIsoformSummary summary, final PrintStream out) {

		Gene gene = summary.getGene();
		String[] line = {gene.getName(), gene.getContig(), Integer.toString(gene.getStart()), Integer.toString(gene.getEnd()), percentageFormat.format(gc.getGCPercent()),
				percentageFormat.format(summary.getMedianGC()), percentageFormat.format(summary.getMedianC()), percentageFormat.format(summary.getMedianG()),
				Integer.toString(summary.getMedianTranscriptLength()), Integer.toString(summary.getNumTranscripts()), Integer.toString(summary.getMedianGQuadruplexes())};
		String h = StringUtils.join(line, "\t");
		out.println(h);
	}

	private void writeHeaderTranscript(final PrintStream out) {
		 String [] header = {"TRANSCRIPT", "CHR", "START", "END", "PCT_GC", "PCT_C", "PCT_G", "TRANSCRIPT_LENGTH", "NUM_GQUADRUPLEXES"};
		 String h = StringUtils.join(header, "\t");
		 out.println(h);
	 }

	public void writeResultTranscript (final List<GCResult> gcList, final PrintStream out) {
		for (GCResult gc: gcList)
			writeResultTranscript(gc, out);
	}

	public void writeResultTranscript(final GCResult gc, final PrintStream out) {
		String [] line = {gc.getTranscript().name, gc.getTranscript().getGene().getContig(), Integer.toString(gc.getTranscript().start()), Integer.toString(gc.getTranscript().end()),
				percentageFormat.format(gc.getGCPercent()), percentageFormat.format(gc.getCPercent()), percentageFormat.format(gc.getGPercent()),
				Integer.toString(gc.getRegionLength()), Integer.toString(gc.getNumGQuadruplexesObserved())};

		String h = StringUtils.join(line, "\t");
		out.println(h);

	}


	/**
	 * Get the median transcript length fromt the gene model.
	 * @param gene
	 * @return
	 */
	/*
	private int getMedianTranscriptLength (final Gene gene) {
		List<Double> transcriptLengths= new ArrayList<Double>();

		for (Transcript t : gene) {
			double length = (t.transcriptionEnd - t.transcriptionStart) +1;
			transcriptLengths.add(length);
		}

		// calculate the median.
		double median = getMedian(transcriptLengths);
		return (int) Math.round(median);

	}
	*/
	/*
	private double getMedian (final List<Double> data) {
		Collections.sort(data);
		int numTranscripts=data.size();
		int middle = numTranscripts/2;
		if (numTranscripts%2==1) {
			double result = data.get(middle);
			return (result);
			//return (int) Math.round();
		}
		else {
			double result = (data.get(middle-1) + data.get(middle))/2;
			return result;
			//return (int) Math.round(result);
		}
	}
	*/


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new GatherGeneGCLength().instanceMain(args));
	}

}
