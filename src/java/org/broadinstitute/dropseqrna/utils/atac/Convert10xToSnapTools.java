/*
 * MIT License
 *
 * Copyright 2022 Broad Institute
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

package org.broadinstitute.dropseqrna.utils.atac;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableOnceIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(
    summary = "Converts 10X ATAC-seq reads to SnapTools format. "
        + "The barcodes are prepended to the read name. For more information see: "
        + "https://github.com/r3fang/SnapATAC/wiki/FAQs#10X_snap "
        + "The functionality is similar to snaptools dex-fastq but also includes filtering of "
        + "reads based on the barcode quality score, removing the 10X spacer, and matching the "
        + "barcodes to an allowlist.",
    oneLineSummary = "Converts 10X ATAC-seq data to SnapTools format",
    programGroup = org.broadinstitute.dropseqrna.cmdline.DropSeq.class
)
public class Convert10xToSnapTools extends CommandLineProgram {

  @Argument(
      doc = "The input 10X ATAC-seq fastq file containing the read data.",
      shortName = "I"
  )
  public File INPUT;

  @Argument(
      doc = "The input 10X ATAC-seq fastq file containing the barcode, usually the 'R2' file.",
      shortName = "BF"
  )
  public File BARCODE_FASTQ;

  @Argument(
      doc = "The minimum PHRED base quality for barcodes.  Any read for barcodes, excluding their "
          + "spacer, with bases below this quality will be dropped.",
      shortName = "BQ",
      optional = true
  )
  public int BARCODE_QUALITY_THRESHOLD = 10;

  @Argument(
      doc = "The list of 10X ATAC barcodes.",
      shortName = "BA"
  )
  public File BARCODE_ATAC;

  @Argument(
      doc = "The list of 10X Gene Expression (GEX) barcodes. If this is not provided, then "
          + "the ATAC barcodes will be used for both ATAC and GEX.",
      shortName = "BG",
      optional = true
  )
  public File BARCODE_GEX;

  @Argument(
      doc = "The output file containing the SnapTools-formatted fastq.",
      shortName = "O"
  )
  public File OUTPUT;

  private static final Log LOG = Log.getInstance(Convert10xToSnapTools.class);

  /**
   * The orientation of the 10X ATAC barcodes relative to the spacers.
   * For more information see:
   * <a href="https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/algorithms/overview#atac_bc"
   * >https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/algorithms/overview#atac_bc</a>
   */
  public enum Orientation {
    FIRST_BASES_FORWARD,
    FIRST_BASES_REVERSE,
    LAST_BASES_FORWARD,
    LAST_BASES_REVERSE
  }

  @Override
  protected int doWork() {
    LOG.info("Loading barcode lists");
    final Map<String, String> originalBarcodeAllowlist =
        loadOriginalBarcodeAllowlist(BARCODE_ATAC, BARCODE_GEX);
    final Map<String, String> mismatchBarcodeAllowlist =
        getMismatchBarcodeAllowlist(originalBarcodeAllowlist);
    final int barcodeLength = getBarcodeLength(originalBarcodeAllowlist);

    final Orientation orientation = detectBarcodeOrientation(
        BARCODE_FASTQ,
        originalBarcodeAllowlist,
        barcodeLength
    );

    prefixInputReads(
        INPUT,
        BARCODE_FASTQ,
        OUTPUT,
        mismatchBarcodeAllowlist,
        orientation,
        barcodeLength,
        BARCODE_QUALITY_THRESHOLD
    );

    return 0;
  }

  /**
   * Returns the barcode allowlist with each ATAC to GEX mapping. The ATAC barcodes will contain Ns
   * for the mismatched bases.
   */
  private static Map<String, String> loadOriginalBarcodeAllowlist(
      final File barcodeAtac,
      final File barcodeGex
  ) {
    try {
      // Resolve the allowlist to enable detection of gzipped file extensions.
      IOUtil.assertFileIsReadable(barcodeAtac);
      IOUtil.assertFileIsReadable(barcodeGex);
      final Map<String, String> barcodeAllowlist = new HashMap<>();
      try (
          final IterableOnceIterator<String> atacReader = IOUtil.readLines(barcodeAtac);
          final IterableOnceIterator<String> gexReader = IOUtil.readLines(barcodeGex)
      ) {
        while (atacReader.hasNext()) {
          if (!gexReader.hasNext()) {
            throw new SAMException("Ran out of GEX barcode records before ATAC records.");
          }
          final String atacBarcode = atacReader.next();
          final String gexBarcode = gexReader.next();
          if (barcodeAllowlist.put(atacBarcode, gexBarcode) != null) {
            throw new SAMException("ATAC barcode duplicated: " + atacBarcode);
          }
        }
        if (gexReader.hasNext()) {
          throw new SAMException("Ran out of ATAC barcode records before GEX records.");
        }
      }
      return barcodeAllowlist;
    } catch (final IOException e) {
      throw new SAMException(e);
    }
  }

  /**
   * Creates a pre-populated list of barcodes with Ns for the mismatched bases.
   */
  private static Map<String, String> getMismatchBarcodeAllowlist(
      final Map<String, String> originalBarcodeAllowlist
  ) {
    final Map<String, String> barcodeAllowlist = new HashMap<>();
    for (final Map.Entry<String, String> entry : originalBarcodeAllowlist.entrySet()) {
      final String atacBarcode = entry.getKey();
      final String gexBarcode = entry.getValue();
      for (int i = 0; i < atacBarcode.length(); i++) {
        final String atacBarcodeMismatch =
            atacBarcode.substring(0, i) + "N" + atacBarcode.substring(i + 1);
        if (barcodeAllowlist.put(atacBarcodeMismatch, gexBarcode) != null) {
          throw new SAMException(
              String.format(
                  "Barcode mismatch %s generated for %s already existed in the allowlist. "
                      + "This occurs when two barcodes differ by a single base. "
                      + "Please validate your barcode list and try again.",
                  atacBarcodeMismatch,
                  atacBarcode
              )
          );
        }
      }
    }
    return barcodeAllowlist;
  }

  /**
   * Returns the ATAC barcode length.
   */
  private static int getBarcodeLength(final Map<String, String> barcodeAllowlist) {
    return barcodeAllowlist
        .keySet()
        .stream()
        .findFirst()
        .orElseThrow(() -> new SAMException("No barcodes found in allowlist."))
        .length();
  }

  /**
   * Detects the orientation of the barcodes and spacers across all the barcodes.
   */
  private static Orientation detectBarcodeOrientation(
      final File barcodeFastq,
      final Map<String, String> originalBarcodeAllowlist,
      final int barcodeLength
  ) {
    LOG.info("Retrieving orientation of barcodes and spacers");
    final ProgressLogger progress = new ProgressLogger(LOG, 10000000, "Processed", "barcodes");
    final long[] counts = new long[Orientation.values().length];
    try (final FastqReader barcodeRecords = new FastqReader(barcodeFastq)) {
      for (final FastqRecord barcodeRecord : barcodeRecords) {
        progress.record("", 0);
        final int barcodeReadLength = barcodeRecord.getReadLength();
        if (barcodeReadLength < barcodeLength) {
          continue;
        }
        for (final Orientation orientation : Orientation.values()) {
          final String barcode = getAtacBarcode(barcodeRecord, orientation, barcodeLength);
          if (originalBarcodeAllowlist.containsKey(barcode)) {
            counts[orientation.ordinal()]++;
          }
        }
      }
    }
    LOG.info("Barcode orientations:");
    for (final Orientation orientation : Orientation.values()) {
      LOG.info(String.format("  %s: %d", orientation, counts[orientation.ordinal()]));
    }
    Orientation orientation = null;
    long maxCount = -1;
    for (final Orientation o : Orientation.values()) {
      if (counts[o.ordinal()] > maxCount) {
        orientation = o;
        maxCount = counts[o.ordinal()];
      }
    }
    LOG.info("Detected barcode orientation: " + orientation);
    return orientation;
  }

  /**
   * Outputs the input reads with the barcodes prefixed to the read name.
   */
  private static void prefixInputReads(
      final File input,
      final File barcodeFastq,
      final File output,
      final Map<String, String> mismatchBarcodeAllowlist,
      final Orientation orientation,
      final int barcodeLength,
      final int barcodeQualityThreshold
  ) {
    LOG.info("Prefixing input reads with barcodes");
    final char barcodeQualityThresholdFastq = SAMUtils.phredToFastq(barcodeQualityThreshold);
    final ProgressLogger progress = new ProgressLogger(LOG, 10000000, "Processed", "reads");
    int writtenBarcodes = 0;
    try (
        final FastqReader fastqReader = new FastqReader(input);
        final FastqReader barcodeRecords = new FastqReader(barcodeFastq);
        final FastqWriter fastqWriter = new BasicFastqWriter(output)
    ) {
      for (final FastqRecord inputRecord : fastqReader) {
        progress.record("", 0);
        if (!barcodeRecords.hasNext()) {
          throw new SAMException("Ran out of barcode records before read records.");
        }
        final FastqRecord barcodeRecord = barcodeRecords.next();
        final int barcodeReadLength = barcodeRecord.getReadLength();

        // Remove reads that are too short to contain a barcode
        if (barcodeReadLength < barcodeLength) {
          continue;
        }
        // Remove reads that have a barcode with a low quality score
        if (hasLowQuality(
            barcodeRecord,
            orientation,
            barcodeLength,
            barcodeQualityThresholdFastq
        )) {
          continue;
        }

        final String atacBarcode = getAtacBarcode(barcodeRecord, orientation, barcodeLength);
        final String gexBarcode = getGexBarcode(atacBarcode, mismatchBarcodeAllowlist);

        if (gexBarcode == null) {
          continue;
        }

        final FastqRecord outputRecord = new FastqRecord(
            gexBarcode + ":" + inputRecord.getReadName(),
            inputRecord.getReadString(),
            inputRecord.getBaseQualityHeader(),
            inputRecord.getBaseQualityString()
        );

        fastqWriter.write(outputRecord);
        writtenBarcodes++;
      }

      LOG.info(
          String.format(
              "Processed %d total reads. Wrote %d and skipped %d.",
              progress.getCount(),
              writtenBarcodes,
              progress.getCount() - writtenBarcodes
          )
      );

      if (barcodeRecords.hasNext()) {
        throw new SAMException("Ran out of read records before barcode records.");
      }
    }
  }

  /**
   * Returns true if any of the bases in the barcode are below the quality threshold.
   */
  private static boolean hasLowQuality(
      final FastqRecord barcodeRecord,
      final Orientation orientation,
      final int barcodeLength,
      final char barcodeQualityThresholdFastq
  ) {
    final String barcodeQuality = barcodeRecord.getBaseQualityString();
    final int start =
        orientation == Orientation.FIRST_BASES_FORWARD || orientation == Orientation.FIRST_BASES_REVERSE
            ? 0
            : barcodeQuality.length() - barcodeLength;
    final int end = start + barcodeLength;
    for (int i = start; i < end; i++) {
      if (barcodeQuality.charAt(i) < barcodeQualityThresholdFastq) {
        return true;
      }
    }
    return false;
  }

  /**
   * Returns the ATAC barcode from the barcode record based on the given orientation.
   */
  @SuppressWarnings("DuplicateExpressions")
  private static String getAtacBarcode(
      final FastqRecord barcodeRecord,
      final Orientation orientation,
      final int barcodeLength
  ) {
    final String barcodeSequence = barcodeRecord.getReadString();
    switch (orientation) {
      case FIRST_BASES_FORWARD:
        return barcodeSequence.substring(0, barcodeLength);
      case FIRST_BASES_REVERSE:
        return SequenceUtil.reverseComplement(barcodeSequence.substring(0, barcodeLength));
      case LAST_BASES_FORWARD:
        return barcodeSequence.substring(barcodeSequence.length() - barcodeLength);
      case LAST_BASES_REVERSE:
        return SequenceUtil.reverseComplement(
            barcodeSequence.substring(barcodeSequence.length() - barcodeLength)
        );
      default:
        throw new SAMException("Unexpected barcode orientation: " + orientation);
    }
  }

  /**
   * Return the GEX barcode for the given ATAC barcode with up to a single mismatch.
   */
  private static String getGexBarcode(
      final String atacBarcode,
      final Map<String, String> barcodeMismatchAllowlist
  ) {
    for (int i = 0; i < atacBarcode.length(); i++) {
      final String atacBarcodeMismatch =
          atacBarcode.substring(0, i) + "N" + atacBarcode.substring(i + 1);
      final String gexBarcode = barcodeMismatchAllowlist.get(atacBarcodeMismatch);
      if (gexBarcode != null) {
        return gexBarcode;
      }
    }

    return null;
  }

  public static void main(String[] args) {
    System.exit(new Convert10xToSnapTools().instanceMain(args));
  }
}
