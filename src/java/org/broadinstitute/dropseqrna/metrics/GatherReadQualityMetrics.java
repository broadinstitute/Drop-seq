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

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        summary = "Calculate the number of reads that are in the BAM, that are mapped, mapped + HQ, mapped + HQ + not PCR duplicated. " +
                "The output lines are ordered by cell barcodes, with the first line containing the summary counts for all the cell barcodes",
        oneLineSummary = "Calculate reads that are in the BAM at different mapping qualities.",
        programGroup = DropSeq.class
)
public class GatherReadQualityMetrics extends CommandLineProgram {
	private final Log log = Log.getInstance(GatherReadQualityMetrics.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.  Must be coordinate sorted.")
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The file to write stats to.")
	public File OUTPUT;

	@Argument(shortName="READ_MQ", doc = "Minimum mapping quality to include the read in the analysis. Set to 0 to not filter reads by map quality.")
	public int MINIMUM_MAPPING_QUALITY = 10;

	@Argument(doc="Optionally aggregate reads by a tag and output per-tag metrics.  The map quality scores histogram will still be computed globally.", optional=true)
	public String TAG=null;

    @Argument(doc="Include non-PF reads when gathering metrics")
    public boolean INCLUDE_NON_PF_READS=false;

    // The key used to output the global metrics collected for all the cell barcodes
    public final static String GLOBAL = "all";

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
		MetricsFile<ReadQualityMetrics, Integer> outFile = new MetricsFile<>();
		outFile.addHistogram(metricsMap.get(this.GLOBAL).getHistogram());
        // Make sure the GLOBAL metrics is added first
        outFile.addMetric(metricsMap.remove(this.GLOBAL));
		for (ReadQualityMetrics metrics: metricsMap.values())
			outFile.addMetric(metrics);
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
        Map<String, ReadQualityMetrics> result = new TreeMap<>();

		ReadQualityMetrics globalMetrics = new ReadQualityMetrics(this.MINIMUM_MAPPING_QUALITY, this.GLOBAL, true);

		SamReader in = SamReaderFactory.makeDefault().open(INPUT);

		for (final SAMRecord r : in)
			if (!r.getReadFailsVendorQualityCheckFlag() || INCLUDE_NON_PF_READS) {
                p.record(r);

                globalMetrics.addRead(r);
                // gather per tag metrics if required.
                result = addMetricsPerTag(r, result);
            }

		CloserUtil.close(in);
		result.put(this.GLOBAL, globalMetrics);
		return (result);
	}

	private Map<String, ReadQualityMetrics> addMetricsPerTag (final SAMRecord r, final Map<String, ReadQualityMetrics> result) {
		// short circuit 1 - no aggregate requested, don't bother.
		if (this.TAG==null) return (result);

		// short circuit 2 - no tag on the read, don't bother.
		String tag = r.getStringAttribute(this.TAG);
		if (tag==null) return (result);

		ReadQualityMetrics m = result.get(tag);
		if (m==null) {
			m = new ReadQualityMetrics(this.MINIMUM_MAPPING_QUALITY, tag, false);
			result.put(tag, m);
		}
		m.addRead(r);
		return (result);
	}





}
