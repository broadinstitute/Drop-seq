/*
 * MIT License
 *
 * Copyright 2023 Broad Institute
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
package org.broadinstitute.dropseqrna.beadsynthesis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.BaseRange;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.FilterBamByTag;
import org.broadinstitute.dropseqrna.utils.readiterators.UnsortedMergingSamRecordIterator;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

@CommandLineProgramProperties(
        summary = "Produce a histogram of the unique sequences for a given {read number, base range}.",
        oneLineSummary = "Produce a histogram of the unique sequences for a given {read number, base range}.",
        programGroup = DropSeq.class
)
public class CountBarcodeSequences
extends CommandLineProgram {
 @Argument(doc="One of more paired-end BAMs or bam_lists.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
 public List<File> INPUT;

 @Argument(doc="From which read in the pair should the barcode sequence be extracted. Must be 1 or 2.")
 public int BARCODED_READ;

 @Argument(doc="Base range to extract, seperated by a dash.  E.g 1-4.  Can extract multiple ranges by separating them by a colon.  For example 1-4:17-22 extracts the first 4 bases, then the 17-22 bases, and glues the sequence together into a single sequence for a tag.")
 public String BASE_RANGE;

 @Argument(doc="Histogram of barcode sequences.", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
 public File OUTPUT;

 @Argument(doc="List of expected barcodes, one per line.  if this is specified, only barcodes in this list will be emitted.",
 optional = true)
 public File ALLOWED_BARCODES;

 @Argument(doc="Add this value to the observed count of every barcode in the allow list.  The effect is to make " +
         "every barcode in the allow list observed at least once.")
 public int ALLOWLIST_PSEUDOCOUNT = 0;

 private final Log log = Log.getInstance(CountBarcodeSequences.class);

 @Override
 protected String[] customCommandLineValidation() {
  IOUtil.assertFileIsWritable(OUTPUT);
  if (ALLOWED_BARCODES != null) {
   IOUtil.assertFileIsReadable(ALLOWED_BARCODES);
  }
  final ArrayList<String> list = new ArrayList<>(1);
  if (BARCODED_READ != 1 && BARCODED_READ != 2) {
   list.add("Invalid BARCODED_READ value: " + BARCODED_READ);
  }
  return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
 }

 @Override
 protected int doWork() {
  INPUT = FileListParsingUtils.expandFileList(INPUT);
  if (INPUT.isEmpty()) {
   throw new RuntimeException("At least one INPUT must be specified.");
  }
  INPUT.forEach(IOUtil::assertFileIsReadable);
  List<BaseRange> baseRanges = BaseRange.parseBaseRange(this.BASE_RANGE);
  final int totalRangeSize = BaseRange.getTotalRangeSize(baseRanges);
  final Set<String> allowedBarcodes;
  if (ALLOWED_BARCODES != null) {
   allowedBarcodes = FilterBamByTag.readValues(ALLOWED_BARCODES);
   allowedBarcodes.forEach(allowedBarcode -> {
            if (allowedBarcode.length() != totalRangeSize) {
             throw new RuntimeException(String.format("BASE_RANGES implies length of %d, but value '%s' in ALLOWED_BARCODES file %s has length %d",
                     totalRangeSize, allowedBarcode, ALLOWED_BARCODES, allowedBarcode.length()));
            }
           }
   );
  } else {
   allowedBarcodes = null;
  }
  final SamReaderFactory readerFactory = SamReaderFactory.makeDefault();
  final List<SamReader> readers = INPUT.stream().map(readerFactory::open).toList();
  // Header is not important, but we need one to iterate
  final SAMFileHeader header = readers.getFirst().getFileHeader();
  final CloseableIterator<SAMRecord> it = new UnsortedMergingSamRecordIterator(header, readers);
  final Histogram<String> histogram = new Histogram<>("SEQUENCE", "COUNT");

  if (allowedBarcodes != null && ALLOWLIST_PSEUDOCOUNT > 0) {
   allowedBarcodes.stream().forEach(bc -> histogram.increment(bc, 1));
  }

  CountBarcodeSequenceMetrics metrics = new CountBarcodeSequenceMetrics();
  ProgressLogger progress = new ProgressLogger(this.log, 10000000);
  while (it.hasNext()) {
   final SAMRecord rec = it.next();
   progress.record(rec);
   if (!rec.getReadPairedFlag()) {
    throw new RuntimeException("Unexpected unpaired read: " + rec);
   }
   if ((rec.getFirstOfPairFlag() && (BARCODED_READ == 2)) || (rec.getSecondOfPairFlag() && (BARCODED_READ == 1))) {
    continue;
   }
   ++metrics.NUM_READ_PAIRS;
   final String barcode = BaseRange.getSequenceForBaseRange(baseRanges, rec.getReadString());
   if (allowedBarcodes != null && !allowedBarcodes.contains(barcode)) {
    ++metrics.NUM_NOT_ALLOWED;
    continue;
   }
   histogram.increment(barcode);
  }
  CloserUtil.close(it);
  final MetricsFile<CountBarcodeSequenceMetrics, String> metricsFile = getMetricsFile();
  metricsFile.addMetric(metrics);
  metricsFile.addHistogram(histogram);
  final Writer metricsWriter = IOUtil.openFileForBufferedWriting(OUTPUT); // Enables gzip writing
  metricsFile.write(metricsWriter);
  try {
   metricsWriter.close();
  } catch (IOException e) {
   throw new RuntimeIOException(e);
  }
  return 0;
 }

 public static class CountBarcodeSequenceMetrics extends MetricBase {
  public long NUM_READ_PAIRS;
  public long NUM_NOT_ALLOWED;
 }
}
