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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class TagBamWithReadSequenceExtendedTest {

	@Test
	public void testHardClip () {
		SAMRecord r = new SAMRecord(null);
		r.setReadString("AAAAAAAAAATTTTTT"); //10A6T
		byte [] quals = {10,10,10,0,0,0,10,10,10,10,0,0,0,0,10,10}; // 0 qual at 4-6, 11-14
		r.setBaseQualities(quals);
		List<BaseRange> ranges = new ArrayList<BaseRange>();
		ranges.add(new BaseRange (4,6));
		ranges.add(new BaseRange (11,14));

		r=TagBamWithReadSequenceExtended.hardClipBasesFromRead(r, ranges);

		String newSeq = r.getReadString();
		byte [] newQuals = r.getBaseQualities();
		String expectedSeq ="AAAAAAATT";
		byte [] expectedQuals = {10,10,10,10,10,10,10,10,10};

		Assert.assertEquals(expectedSeq.length(), newSeq.length());
		Assert.assertEquals(expectedSeq, newSeq);
		for (int i=0; i<expectedQuals.length;i++)
			Assert.assertEquals(expectedQuals[i], newQuals[i]);
	}

	private static class SamRecords {
		final File samFile = File.createTempFile("TagBamWithReadSequenceExtended.", ".sam");
		final SAMFileHeader header;
		final ArrayList<SAMRecord> records = new ArrayList<>();
		SamRecords(boolean paired, boolean sorted) throws IOException {
		    samFile.deleteOnExit();
			final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(sorted, sorted? SAMFileHeader.SortOrder.queryname: SAMFileHeader.SortOrder.unsorted);
			if (paired) {
				builder.addUnmappedPair("pair1");
				builder.addUnmappedPair("pair2");
			} else {
				builder.addUnmappedFragment("frag1");
				builder.addUnmappedFragment("frag2");
			}
			final SamReader reader = builder.getSamReader();
			header = reader.getFileHeader();
			final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(header, true, samFile, null);
			for (final SAMRecord rec : reader) {
				records.add(rec);
				writer.addAlignment(rec);
			}
			writer.close();
		}
	}

	private Map.Entry<String, String> getBarcodeAndNonBarcode(byte[] bases, boolean singleRange,
                                                              int rangeStart, int rangeEnd,
                                                              int secondRangeStart, int secondRangeEnd) {
	    final StringBuilder barcode = new StringBuilder();
        final StringBuilder nonBarcode = new StringBuilder();
        for (int i = 0; i < bases.length; ++i) {
            final char theChar = (char)bases[i];
            if ((i+1 >= rangeStart && i+1 <= rangeEnd) ||
                    (!singleRange && i+1 >= secondRangeStart && i+1 <= secondRangeEnd)) {
                barcode.append(theChar);
            } else {
                nonBarcode.append(theChar);
            }
        }
        return new Map.Entry<String, String>() {
            @Override
            public String getKey() {
                return barcode.toString();
            }

            @Override
            public String getValue() {
                return nonBarcode.toString();
            }

            @Override
            public String setValue(String value) {
                return null;
            }
        };
    }

    private enum ReadTreatment {
        Discard, Clip, Keep
    }

    private enum WhichRead {
        First, Second, Unpaired
    }

    @Test(dataProvider = "testTagBamWithReadSequenceExtendedDataProvider")
	public void testTagBamWithReadSequenceExtended(
	        ReadTreatment readTreatment,
            WhichRead whichRead,
            boolean singleRange,
            boolean sorted,
            boolean goodQuality,
            boolean tagBarcodedRead
    ) throws IOException {
        System.err.println(String.format("read treatment: %s; which read: %s; single range: %s; sorted: %s; good quality: %s;" +
                "tag barcoded read: %s", readTreatment, whichRead, singleRange, sorted, goodQuality, tagBarcodedRead));
        final SamRecords samRecords = new SamRecords(whichRead != WhichRead.Unpaired, sorted);
        final TagBamWithReadSequenceExtended clp = new TagBamWithReadSequenceExtended();
        clp.INPUT = samRecords.samFile;
        clp.TAG_FULL_QUALITY = "XY";
        clp.OUTPUT = File.createTempFile("tagged.", ".sam");
        clp.OUTPUT.deleteOnExit();
        clp.SUMMARY = File.createTempFile("TagBamWithReadSequenceExtended.", ".tag_summary.txt");
        clp.SUMMARY.deleteOnExit();
        int rangeStart = 1;
        int rangeEnd = 4;
        clp.BASE_RANGE = rangeStart + "-" + rangeEnd;
        int secondRangeStart = 10;
        int secondRangeEnd = 20;
        if (!singleRange) {
            clp.BASE_RANGE += ":" + secondRangeStart + "-" + secondRangeEnd;
        }
        if (whichRead == WhichRead.Second) {
            clp.BARCODED_READ = 2;
        } else {
            clp.BARCODED_READ = 1;
        }
        clp.TAG_BARCODED_READ = tagBarcodedRead;
        clp.DISCARD_READ = (readTreatment == ReadTreatment.Discard);
        clp.HARD_CLIP_BASES = (readTreatment == ReadTreatment.Clip);
        if (goodQuality) {
            clp.BASE_QUALITY = 0;
            clp.NUM_BASES_BELOW_QUALITY = 35; // SamRecordSetBuilder creates 36-base reads
        } else {
            clp.BASE_QUALITY = 52; // SamRecordSetBuilder assigns random qualitys up to 50
            clp.NUM_BASES_BELOW_QUALITY = 2;
        }
        Assert.assertEquals(0, clp.doWork());
        final ArrayList<SAMRecord> taggedRecords = new ArrayList<>();
        final SamReader samReader = SamReaderFactory.makeDefault().open(clp.OUTPUT);
        samReader.iterator().stream().forEach(taggedRecords::add);
        CloserUtil.close(samReader);

        if (readTreatment == ReadTreatment.Discard) {
            Assert.assertEquals(taggedRecords.size(), 2);
        } else {
            Assert.assertEquals(taggedRecords.size(), whichRead == WhichRead.Unpaired? 2: 4);
        }
        // Just check the first pair or fragment
        final SAMRecord barcodedRead;

        final int firstIndex;
        final int secondIndex;
        if (readTreatment == ReadTreatment.Discard) {
            barcodedRead = null;
            firstIndex = -1;
            secondIndex = -1;
        } else if (whichRead == WhichRead.Unpaired) {
            firstIndex = -1;
            secondIndex = -1;
            barcodedRead = taggedRecords.get(0);
        } else {
            if (taggedRecords.get(0).getFirstOfPairFlag()) {
                firstIndex = 0;
                secondIndex = 1;
            } else {
                firstIndex = 1;
                secondIndex = 0;
            }
            barcodedRead = (whichRead == WhichRead.Second ? taggedRecords.get(secondIndex) : taggedRecords.get(firstIndex));
        }
        final SAMRecord originalBarcodedRead = (whichRead == WhichRead.Second? samRecords.records.get(1): samRecords.records.get(0));
        final Map.Entry<String, String> barcodeAndNonBarcode =
                getBarcodeAndNonBarcode(originalBarcodedRead.getReadBases(), singleRange, rangeStart, rangeEnd, secondRangeStart, secondRangeEnd);
        final String barcode = barcodeAndNonBarcode.getKey();
        final String nonBarcode = barcodeAndNonBarcode.getValue();

        if (readTreatment == ReadTreatment.Clip) {
            // Confirm that clipped read has the right stuff
            Assert.assertEquals(barcodedRead.getReadString(), nonBarcode);
        } else if (readTreatment == ReadTreatment.Keep) {
            // Confirm that read is not clipped.
            Assert.assertEquals(barcodedRead.getReadString(), originalBarcodedRead.getReadString());
        }
        final SAMRecord taggedRead;
        if (whichRead == WhichRead.Unpaired) {
            taggedRead = taggedRecords.get(0);
        } else if (readTreatment == ReadTreatment.Discard) {
            taggedRead = taggedRecords.get(0);
        } else if (tagBarcodedRead) {
            if (whichRead == WhichRead.First) {
                taggedRead = taggedRecords.get(firstIndex);
            } else {
                taggedRead = taggedRecords.get(secondIndex);
            }
        } else {
            if (whichRead == WhichRead.First) {
                taggedRead = taggedRecords.get(secondIndex);
            } else {
                taggedRead = taggedRecords.get(firstIndex);
            }
        }
        // Confirm that barcode tag contains the barcode
        Assert.assertEquals(taggedRead.getStringAttribute(clp.TAG_NAME), barcode);

        // Check quality flag
        if (goodQuality) {
            Assert.assertNull(taggedRead.getAttribute(clp.TAG_QUALITY));
        } else {
            int barcodeBasesBelowQualityThreshold = barcode.length();
            Assert.assertEquals(taggedRead.getIntegerAttribute(clp.TAG_QUALITY).intValue(), barcodeBasesBelowQualityThreshold);
        }

        // check that the reconstructed TAG_QUALITY from TAG_FULL_QUALITY equals expected TAG_QUALITY
        byte[] qual = (byte[]) taggedRead.getAttribute(clp.TAG_FULL_QUALITY);
        int numBadBases = 0;
        for (int i = 0; i < qual.length; i++) {
            byte q = qual[i];
            if (q < clp.BASE_QUALITY)
                numBadBases++;
        }
        if (numBadBases>=clp.NUM_BASES_BELOW_QUALITY) {
            Assert.assertEquals(taggedRead.getIntegerAttribute(clp.TAG_QUALITY).intValue(), numBadBases);
        }else {
            Assert.assertNull(taggedRead.getAttribute(clp.TAG_QUALITY));
        }
	}

	@DataProvider(name="testTagBamWithReadSequenceExtendedDataProvider")
    public Object[][] testTagBamWithReadSequenceExtendedDataProvider() {
	    final ArrayList<Object[]> ret = new ArrayList<>();
	    final boolean[] booleanValues = new boolean[]{true, false};
	    for (final WhichRead which: WhichRead.values()) {
	        final ReadTreatment[] possibleTreatements;
	        if (which == WhichRead.Unpaired) {
	            possibleTreatements = new ReadTreatment[]{ReadTreatment.Clip, ReadTreatment.Keep};
            } else {
	            possibleTreatements = ReadTreatment.values();
            }
            for (final ReadTreatment treatment: possibleTreatements) {
	            for (final boolean singleRange: booleanValues) {
	                for (final boolean sorted: booleanValues) {
	                    for (final boolean goodQuality: booleanValues) {
	                        final boolean[] tagBarcodeReadChoices;
	                        if (treatment == ReadTreatment.Discard) {
	                            tagBarcodeReadChoices = new boolean[]{false};
                            } else if (which == WhichRead.Unpaired) {
	                            tagBarcodeReadChoices = new boolean[]{true};
                            } else {
	                            tagBarcodeReadChoices = booleanValues;
                            }
                            for (final boolean tagBarcodeRead: tagBarcodeReadChoices) {
	                            ret.add(new Object[]{
	                                    treatment,
                                        which,
                                        singleRange,
                                        sorted,
                                        goodQuality,
                                        tagBarcodeRead
                                });
                            }
                        }
                    }
                }
            }
        }
        return ret.toArray(new Object[ret.size()][]);
    }
}
