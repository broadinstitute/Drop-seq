package org.broadinstitute.dropseqrna.readtrimming;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.util.Arrays;

import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(usage = "", usageShort = "", programGroup = DropSeq.class)
public class PolyATrimmer extends CommandLineProgram {

private final Log log = Log.getInstance(PolyATrimmer.class);
	
	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.")
	public File INPUT;
	
	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM file")
	public File OUTPUT;
	
	@Option(doc = "The output summary statistics", optional=true)
	public File OUTPUT_SUMMARY;
		
	@Option(doc="How many mismatches are acceptable in the sequence.")
	public Integer MISMATCHES=1; 
	
	@Option(doc="How many bases of polyA qualifies as a run of A's.")
	public Integer NUM_BASES=6; 
	
	@Option (doc="The tag to set for trimmed reads.  This tags the first base to exclude in the read.  37 would mean to retain the first 36 bases.")
	public String TRIM_TAG="ZP";
	
	private Integer readsTrimmed=0;
	private Histogram<Integer> numBasesTrimmed= new Histogram<Integer>();
	
	@Override
	protected int doWork() {
		
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		final ProgressLogger progress = new ProgressLogger(log);
		
		SamReader bamReader = SamReaderFactory.makeDefault().open(INPUT);
		SAMFileHeader header = bamReader.getFileHeader();
        SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);
        PolyAFinder polyAFinder = new PolyAFinder(this.NUM_BASES, this.MISMATCHES);
        
        for (SAMRecord r: bamReader) {        
        	SAMRecord rr= hardClipPolyAFromRecord(r, polyAFinder);
        	writer.addAlignment(rr);
        	progress.record(r);
        }
        CloserUtil.close(bamReader);
        writer.close();
        log.info("Number of reads trimmed: ", this.readsTrimmed);
        if (this.OUTPUT_SUMMARY!=null) writeSummary(this.numBasesTrimmed); 
        
		return 0;
	}
	
	/**
	 * Hard clip out reads that have polyA runs.  Finds the longest sequence of polyAs, and clips all bases that occur in or after that.
	 * @param r The read to trim
	 * @param polyAFinder The polyAFinder configured to find polyA runs at a min number of bases
	 * and with at most some number of errors.
	 * @return
	 */
	SAMRecord hardClipPolyAFromRecord (SAMRecord r, PolyAFinder polyAFinder) {
		
		int readLength=r.getReadLength();
		
		String readString = r.getReadString();
		int polyAStart =polyAFinder.getPolyAStart(readString);
		
		// short circuit if the template wasn't found.
		if (polyAStart==-1) {
			return (r);
		} 
		this.readsTrimmed++;
		
		this.numBasesTrimmed.increment(polyAStart);
		// terrible luck.  your read is just a wad of A's.
		if (polyAStart==0) {
			// attempt a work around for reads that would be 0 length after trimming.
			// instead of trimming the barcode to a 0 length read, set the base qualities to be low.
			byte [] value= new byte [readLength];
			Arrays.fill(value, (byte) 3);
			r.setBaseQualities(value);
			return (r);
			// log.info("STOP");
		}
		
		byte [] read = r.getReadBases();
		read=Arrays.copyOfRange(read, 0, polyAStart);
		r.setReadBases(read);
		// String after = r.getReadString();
		
		byte [] quality = r.getBaseQualities();
		quality=Arrays.copyOfRange(quality, 0, polyAStart);
		r.setBaseQualities(quality);
		r.setAttribute("ZP", polyAStart+1);
		return (r);
	}
	
	private void writeSummary (Histogram<Integer> h) {
		
		MetricsFile<TrimMetric, Integer> mf = new MetricsFile<TrimMetric, Integer>();
		mf.addHistogram(h);
		TrimMetric tm=new TrimMetric(h);
		mf.addMetric(tm);
		mf.write(this.OUTPUT_SUMMARY);
	}
	
	public class TrimMetric extends MetricBase {
		public Double mean;
		public Double stdev;
		
		public TrimMetric (Histogram<Integer> h) {
			mean=h.getMean();
			stdev=h.getStandardDeviation();
		}
		
		public Double getMean() {
			return mean;
		}

		public Double getStdev() {
			return stdev;
		}		
	}
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new PolyATrimmer().instanceMain(args));
	}

	
}


