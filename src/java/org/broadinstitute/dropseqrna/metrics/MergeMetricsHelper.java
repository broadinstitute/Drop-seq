/*
 * MIT License
 *
 * Copyright 2021 Broad Institute
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

import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * Helper class to make implementing metrics mergers easier.
 * Histogram hasn't been included, but could be.
 * Currently this only supports files containing a single metric, but can be enhanced as necessary.
 * @param <METRIC_CLASS>
 */
public class  MergeMetricsHelper<METRIC_CLASS extends MetricBase> {

 /**
  * Merge input metric files where each input contains a single metric object
  * @param OUTPUT output file to be written
  * @param INPUT input files to be read and accumulated.
  * @param ctor creates the object to be merged into, with non-accumulating fields populated
  * @param merge function for merging a source metric into target merge(target, source)
  * @param outputMetricsFile If null, one is created to write metrics to OUTPUT.  Caller can provide a MetricsFile
  *                          that has been populated, e.g. by CommandLineProgram.getMetricsFile()
  */
 public void mergeSingleMetric(final File OUTPUT, final List<File> INPUT, final boolean DELETE_INPUTS,
                               Supplier<? extends METRIC_CLASS> ctor, BiConsumer<METRIC_CLASS, METRIC_CLASS> merge,
                               MetricsFile outputMetricsFile) {
  IOUtil.assertFileIsWritable(OUTPUT);
  try {
   final METRIC_CLASS outMetric = ctor.get();
   for (File file : INPUT) {
    MetricsFile<METRIC_CLASS, Integer> metricsFile = new MetricsFile<>();
    metricsFile.read(IOUtil.openFileForBufferedReading(file));
    for (Header header : metricsFile.getHeaders()) {
     outputMetricsFile.addHeader(header);
    }
    final METRIC_CLASS inMetric = metricsFile.getMetrics().getFirst();
    merge.accept(outMetric, inMetric);
   }

   if (outputMetricsFile == null) {
    outputMetricsFile = new MetricsFile<>();
   }
   outputMetricsFile.addMetric(outMetric);
   outputMetricsFile.write(OUTPUT);
  } catch (Exception e) {
   throw new RuntimeException(e);
  }
  if (DELETE_INPUTS) {
   final List<String> couldNotBeDeleted = INPUT.stream().filter(f -> !f.delete()).map(File::getAbsolutePath).collect(Collectors.toList());
   if (!couldNotBeDeleted.isEmpty()) {
    throw new RuntimeIOException("Unabled to delete files: " + StringUtil.join(", ", couldNotBeDeleted));
   }
  }
 }
}
