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

import org.apache.commons.lang.StringUtils;
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

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input histogram reports")
	public List<File> INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output merged histogram report file.")
	public File OUTPUT;

	@Override
	protected int doWork() {
		IOUtil.assertFileIsWritable(OUTPUT);

        ObjectCounter <String> mergedCounter = null;
        for (File reportFile : INPUT) {
            ObjectCounter <String> counter = ObjectCounter.readReportFile(reportFile);
            if (mergedCounter == null) {
                mergedCounter = counter;
            } else {
                mergedCounter.increment(counter);
            }
        }
        if (mergedCounter == null) {
            throw new RuntimeException("mergedCounter is null");
        }

		PrintStream writer = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));
		writeHeader(writer);

        List<String> keysByCount = mergedCounter.getKeysOrderedByCount(true);
        for (String key : keysByCount) {
			int count = mergedCounter.getCountForKey(key);
			String[] fields = {count+"", key};
			String result = StringUtils.join(fields, "\t");
			writer.println(result);
			writer.flush();
		}
		writer.close();

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
}
