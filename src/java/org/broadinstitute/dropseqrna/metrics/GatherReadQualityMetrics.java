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
package org.broadinstitute.dropseqrna.metrics;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        usage = "Calculate the number of reads that are in the BAM, that are mapped, mapped + HQ, mapped + HQ + not PCR duplicated",
        usageShort = "Calculate reads that are in the BAM at different mapping qualities.",
        programGroup = DropSeq.class
)
public class GatherReadQualityMetrics extends CommandLineProgram {
	private final Log log = Log.getInstance(GatherReadQualityMetrics.class);

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.  Must be coordinate sorted.")
	public File INPUT;

	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The file to write stats to.")
	public File OUTPUT;
	
	@Option(doc="The minimum map quality for a read to be considered high quality")
	public Integer MAP_QUALITY=10;
	
	@Option(doc="Optionally aggregate reads by a tag and output per-tag metrics.  The map quality scores histogram will still be computed globally.", optional=true)
	public String TAG=null;

    @Option(doc="Include non-PF reads when gathering metrics")
    public boolean INCLUDE_NON_PF_READS=false;
	
	private String GLOBAL="all";
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new GatherReadQualityMetrics().instanceMain(args));
	}

	/**
	 * Do the work after command line has been parsed. RuntimeException may be
	 * thrown by this method, and are reported appropriately.
	 * 
	 * @return program exit status.
	 */
	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		Map<String, ReadQualityMetrics> metricsMap = gatherMetrics(INPUT);
		MetricsFile<ReadQualityMetrics, Integer> outFile = new MetricsFile<ReadQualityMetrics, Integer>();
		outFile.addHistogram(metricsMap.get(this.GLOBAL).getHistogram());
		for (ReadQualityMetrics metrics: metricsMap.values()) {
			outFile.addMetric(metrics);
		}
		BufferedWriter w = IOUtil.openFileForBufferedWriting(OUTPUT);
		outFile.write(w);
		// close properly.
		try {
			w.close();
		} catch (IOException io) {
			throw new TranscriptomeException("Problem writing file", io);
		}
		return 0;
	}

	public Map<String, ReadQualityMetrics> gatherMetrics(final File inputSamOrBamFile) {
		ProgressLogger p = new ProgressLogger(this.log);
		Map<String, ReadQualityMetrics> result = new HashMap<String, ReadQualityMetrics>();
		
		ReadQualityMetrics globalMetrics = new ReadQualityMetrics(this.MAP_QUALITY, this.GLOBAL, true);
				
		SamReader in = SamReaderFactory.makeDefault().open(INPUT);
		
		for (final SAMRecord r : in) {
            if (!r.getReadFailsVendorQualityCheckFlag() || INCLUDE_NON_PF_READS) {
                p.record(r);

                globalMetrics.addRead(r);
                // gather per tag metrics if required.
                result = addMetricsPerTag(r, result);
            }
		}
		
		CloserUtil.close(in);
		result.put(this.GLOBAL, globalMetrics); 
		return (result);
	}
	
	private Map<String, ReadQualityMetrics> addMetricsPerTag (SAMRecord r, Map<String, ReadQualityMetrics> result) {
		// short circuit 1 - no aggregate requested, don't bother.
		if (this.TAG==null) return (result);
		
		// short circuit 2 - no tag on the read, don't bother.
		String tag = r.getStringAttribute(this.TAG);
		if (tag==null) return (result);
		
		ReadQualityMetrics m = result.get(tag);
		if (m==null) {
			m = new ReadQualityMetrics(this.MAP_QUALITY, tag, false);
			result.put(tag, m);
		}
		m.addRead(r);
		return (result);
	}
	
	
	
	

}
