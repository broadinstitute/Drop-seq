/*
 * MIT License
 *
 * Copyright 2020 Broad Institute
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
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

/**
 * @author skashin
 *
 */
@CommandLineProgramProperties(summary = "Merges histogram reports for a given tag",
        oneLineSummary = "Merges histogram reports for a given tag",
        programGroup = DropSeq.class)
public class MergeBamTagHistograms extends CommandLineProgram {

	private static final Log log = Log.getInstance(MergeBamTagHistograms.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input histogram reports to be merged.", minElements = 1)
	public List<File> INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The merged histogram report.")
	public File OUTPUT;

	@Override
	protected int doWork() {
		IOUtil.assertFileIsWritable(OUTPUT);

        List<String> headers = null;
        ObjectCounter <String> mergedCounter = null;
        for (File reportFile : INPUT) {
            Report report = readReportFile(reportFile);
            if (mergedCounter == null) {
                headers = report.getHeaders();
                mergedCounter = report.getCounter();
            } else {
                headers.addAll(report.getHeaders());
                mergedCounter.increment(report.getCounter());
            }
        }

		PrintStream writer = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));
        for (String header : headers) {
            writer.println(header);
        }
		writeHeader(writer);
        BamTagHistogram.writeHistogram(mergedCounter, writer);

		return 0;
	}

	public void writeHeader(final PrintStream writer) {
		List<String> headerList = new ArrayList<>();
        for (File reportFile : INPUT)
		    headerList.add("INPUT=" + reportFile.getAbsolutePath());
		String header = StringUtils.join(headerList, "\t");
		writer.print("#");
		writer.println(header);
	}

    /**
     * @param reportFile The report file to read.
     */
    public static Report readReportFile(final File reportFile) {
        List<String> headers = new ArrayList<>();
        ObjectCounter <String> counter = new ObjectCounter<>();

        for (String line : IOUtil.readLines(reportFile)) {
            line = line.trim();
            if (line.startsWith("#")) {
                headers.add(line);
                continue;
            }
            String[] strLine = line.split("\t");
            int count = Integer.parseInt(strLine[0]);
            String key = strLine[1];
            counter.incrementByCount(key, count);
        }

        return new Report(headers, counter);
    }

    static class Report {
        private List<String> headers;
        private ObjectCounter<String> counter;

        Report(List<String> headers, ObjectCounter<String> counter) {
            this.headers = headers;
            this.counter = counter;
        }

        public List<String> getHeaders() {
            return headers;
        }

        public ObjectCounter<String> getCounter() {
            return counter;
        }
    }
}
