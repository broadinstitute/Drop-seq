/*
 * MIT License
 *
 * Copyright 2019 Broad Institute
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

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.*;
import org.broadinstitute.dropseqrna.utils.editdistance.MapBarcodesByEditDistance;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionIteratorWrapper;
import picard.cmdline.StandardOptionDefinitions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "Reads a BAM and analyzes groups of reads that share the same value for " +
        "COLLAPSE_TAG.  Within such a group, reads for which UNCOLLAPSED_TAG value is the same as COLLAPSE_TAG value " +
        "are put into the parent subgroup.  The remaining reads are grouped by UNCOLLAPSED_TAG value into child " +
        "subgroups.  For each read, the values for the COUNT_TAG form a tuple.  For each subgroup, the set of unique " +
        "tuples are accumulated.  For each child subgroup, the size of the intersection of the tuple set with the " +
        "parent tuple set is computed.  A metrics file is produced for each child subgroup.",
        oneLineSummary = "Computes overlap of COUNT_TAG values between parent and child groups of reads.",
        programGroup = DropSeq.class)
public class ComputeUMISharing
        extends GeneFunctionCommandLineBase {

    private static final Log log = Log.getInstance(ComputeUMISharing.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input SAM or BAM file to analyze.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The metrics file to be written.")
    public File OUTPUT;

    @Argument(doc="The tag used to groups reads.")
    public String COLLAPSE_TAG;

    @Argument(doc="The tag used to create subgroups of reads within a group.")
    public String UNCOLLAPSED_TAG;

    @Argument(minElements = 1, doc="The tag(s) whose values are used to create the set of unique tuples in each subgroup.")
    public List<String> COUNT_TAG;

    @Argument(doc="The edit distance for comparing COUNT_TAG values when comparing parent and child tuple sets.  " +
            "the nth EDIT_DISTANCE corresponds to the nth COUNT_TAG.  If there are fewer EDIT_DISTANCEs than " +
            "COUNT_TAGs, the remainder are assumed to have EDIT_DISTANCE=0.")
    public List<Integer> EDIT_DISTANCE;

    @Argument(doc = "Should indels be considered in edit distance calculations?  Doing this correctly is far slower " +
            "than a simple edit distance test, but is a more aggressive method that may be useful in some situations.")
    public boolean FIND_INDELS=false;

    @Argument(doc="Number of threads to use.  Only relevant if there is non-zero EDIT_DISTANCE.")
    public int NUM_THREADS=1;


    private ParentEditDistanceMatcher parentEditDistanceMatcher;

    public static void main(final String[] argv) {
        new ComputeUMISharing().instanceMainWithExit(argv);
    }

    public ComputeUMISharing() {
        LOCUS_FUNCTION_LIST.clear();
    }

    @Override
    protected String[] customCommandLineValidation() {
        String[] superMessages = super.customCommandLineValidation();
        if (EDIT_DISTANCE.size() > COUNT_TAG.size()) {
            superMessages = CustomCommandLineValidationHelper.makeValue(superMessages,
                    Collections.singletonList(
                            String.format("More values specified for EDIT_DISTANCE(%d) than COUNT_TAG(%d)",
                                    EDIT_DISTANCE.size(), COUNT_TAG.size())));
        }
        if (EDIT_DISTANCE.stream().anyMatch((editDistance) -> editDistance < 0)) {
            superMessages = CustomCommandLineValidationHelper.makeValue(superMessages,
                    Collections.singletonList("EDIT_DISTANCEs must be >= 0"));
        }
        return superMessages;
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        // Make sure this is modifiable
        EDIT_DISTANCE = new ArrayList<>(EDIT_DISTANCE);
        while (EDIT_DISTANCE.size() < COUNT_TAG.size()) {
            EDIT_DISTANCE.add(0);
        }
        parentEditDistanceMatcher = new ParentEditDistanceMatcher();

        SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
        final ProgressLogger progressLogger = new ProgressLogger(log, 1000000, "Sorting");
        Iterator<SAMRecord> iter = reader.iterator();
        if (LOCUS_FUNCTION_LIST.size() > 0) {
            iter = new GeneFunctionIteratorWrapper(iter, this.GENE_NAME_TAG,
                    this.GENE_STRAND_TAG, this.GENE_FUNCTION_TAG, false, this.STRAND_STRATEGY,
                    this.LOCUS_FUNCTION_LIST);
        }

        CloseableIterator<SAMRecord> sortedIter = SortingIteratorFactory.create(SAMRecord.class, iter,
                PARENT_CHILD_COMPARATOR, new BAMRecordCodec(reader.getFileHeader()), MAX_RECORDS_IN_RAM,
                (SortingIteratorFactory.ProgressCallback<SAMRecord>) progressLogger::record);

        PeekableIterator<List<SAMRecord>> subgroupIterator =
                new PeekableIterator<List<SAMRecord>>(new GroupingIterator<SAMRecord>(
                        new ProgressLoggingIterator(sortedIter, new ProgressLogger(log, 1000000, "Grouping")),
                        GROUPING_COMPARATOR));

        MetricsFile<UmiSharingMetrics, Integer> outFile = getMetricsFile();
        List<SAMRecord> parentSubgroup = null;
        Set<TagValues> parentTuples = new HashSet<>();

        while (subgroupIterator.hasNext()) {
            if (parentSubgroup == null ||
            !parentSubgroup.get(0).getAttribute(COLLAPSE_TAG).equals(subgroupIterator.peek().get(0).getAttribute(COLLAPSE_TAG))) {
                parentSubgroup = subgroupIterator.next();
                parentTuples.clear();
                parentSubgroup.stream().map(parentEditDistanceMatcher::getValues).forEach(parentTuples::add);

            } else {
                final Set<TagValues> childTuples = new HashSet<>();
                final List<SAMRecord> childSubgroup = subgroupIterator.next();
                childSubgroup.stream().map(parentEditDistanceMatcher::getValues).forEach(childTuples::add);
                final UmiSharingMetrics metrics = new UmiSharingMetrics();
                metrics.PARENT = parentSubgroup.get(0).getAttribute(COLLAPSE_TAG).toString();
                metrics.CHILD = childSubgroup.get(0).getAttribute(UNCOLLAPSED_TAG).toString();
                metrics.NUM_PARENT = parentTuples.size();
                metrics.NUM_CHILD = childTuples.size();
                metrics.NUM_SHARED = computeNumShared(parentTuples, childTuples);
                metrics.FRAC_SHARED = metrics.NUM_SHARED/(double)metrics.NUM_CHILD;
                outFile.addMetric(metrics);
            }
        }
        BufferedWriter w = IOUtil.openFileForBufferedWriting(OUTPUT);
        outFile.write(w);
        try {
            w.close();
        } catch (IOException e) {
            throw new RuntimeIOException("Problem writing " + OUTPUT.getAbsolutePath(), e);
        }
        CloserUtil.close(reader);
        return 0;
    }

    private int computeNumShared(final Set<TagValues> parentTuples, final Set<TagValues> childTuples) {
        return (int)childTuples.stream().filter((childTuple) -> parentEditDistanceMatcher.isShared(childTuple, parentTuples)).count();
    }

    private static class TagValues {
        final Object[] zeroEditDistanceValues;
        final String[] nonZeroEditDistanceValues;

        TagValues(final int numZeroEditDistanceTags, final int numNonZeroEditDistanceTags) {
            zeroEditDistanceValues = new Object[numZeroEditDistanceTags];
            nonZeroEditDistanceValues = new String[numNonZeroEditDistanceTags];
        }

        boolean zeroEditDistanceEquals(final TagValues other) {
            return Arrays.equals(this.zeroEditDistanceValues, other.zeroEditDistanceValues);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            TagValues tagValues = (TagValues) o;
            return Arrays.equals(zeroEditDistanceValues, tagValues.zeroEditDistanceValues) &&
                    Arrays.equals(nonZeroEditDistanceValues, tagValues.nonZeroEditDistanceValues);
        }

        @Override
        public int hashCode() {
            int result = Arrays.hashCode(zeroEditDistanceValues);
            result = 31 * result + Arrays.hashCode(nonZeroEditDistanceValues);
            return result;
        }
    }

    private class ParentEditDistanceMatcher {
        private final int[] zeroEditDistanceIndices;
        private final int[] nonZeroEditDistanceIndices;
        private final MapBarcodesByEditDistance med = new MapBarcodesByEditDistance(false, NUM_THREADS, 0);

        ParentEditDistanceMatcher() {
            final ArrayList<Integer> zeroEDs = new ArrayList<>();
            final ArrayList<Integer> nonZeroEDs = new ArrayList<>();
            for (int i = 0; i < EDIT_DISTANCE.size(); ++i) {
                if (EDIT_DISTANCE.get(i) == 0) {
                    zeroEDs.add(i);
                } else {
                    nonZeroEDs.add(i);
                }
            }
            zeroEditDistanceIndices = zeroEDs.stream().mapToInt(i->i).toArray();
            nonZeroEditDistanceIndices = nonZeroEDs.stream().mapToInt(i->i).toArray();

        }

        TagValues getValues(final SAMRecord rec) {
            final TagValues ret = new TagValues(zeroEditDistanceIndices.length, nonZeroEditDistanceIndices.length);
            for (int i = 0; i < zeroEditDistanceIndices.length; ++i) {
                ret.zeroEditDistanceValues[i] = rec.getAttribute(COUNT_TAG.get(zeroEditDistanceIndices[i]));
            }
            for (int i = 0; i < nonZeroEditDistanceIndices.length; ++i) {
                ret.nonZeroEditDistanceValues[i] = rec.getStringAttribute(COUNT_TAG.get(nonZeroEditDistanceIndices[i]));
            }
            return ret;
        }

        boolean isShared(final TagValues childTuple, final Set<TagValues> parentTuples) {
            if (nonZeroEditDistanceIndices.length == 0) {
                return parentTuples.contains(childTuple);
            }

            // Accumulate parents for which ED-zero values match the child.
            final List<TagValues> exactMatchParentTuples = parentTuples.stream().
                    filter(childTuple::zeroEditDistanceEquals).collect(Collectors.toList());

            if (exactMatchParentTuples.isEmpty()) {
                return false;
            }

            if (nonZeroEditDistanceIndices.length == 1) {
                // Fast way for single tag
                final List<String> parentValues = exactMatchParentTuples.stream().map((tuple) -> tuple.nonZeroEditDistanceValues[0]).collect(Collectors.toList());
                return matchesWithinEditDistance(childTuple.nonZeroEditDistanceValues[0], parentValues, EDIT_DISTANCE.get(nonZeroEditDistanceIndices[0]));
            }

            // Unfortunately, if more than one ED-nonzero tag, needs to test each individually.
            for (final TagValues parentTuple : exactMatchParentTuples) {
                boolean matched = true;
                for (int i = 0; i < nonZeroEditDistanceIndices.length; ++i) {
                    if (!matchesWithinEditDistance(
                            childTuple.nonZeroEditDistanceValues[i],
                            parentTuple.nonZeroEditDistanceValues[i],
                            EDIT_DISTANCE.get(nonZeroEditDistanceIndices[i]))) {
                        matched = false;
                        break;
                    }
                }
                // Matched on every ED-nonzero tag
                if (matched) {
                    return true;
                }
            }
            return false;
        }

        private boolean matchesWithinEditDistance(final String barcode, final List<String> comparisonBarcodes, final int editDistance) {
            return !med.processSingleBarcodeMultithreaded(barcode, comparisonBarcodes, FIND_INDELS, editDistance).isEmpty();
        }

        private boolean matchesWithinEditDistance(final String barcode, final String comparisonBarcode, final int editDistance) {
            return matchesWithinEditDistance(barcode, Collections.singletonList(comparisonBarcode), editDistance);
        }

    }

    private final Comparator<SAMRecord> GROUPING_COMPARATOR = new Comparator<SAMRecord>() {
        @Override
        public int compare(final SAMRecord r1, final SAMRecord r2) {
            final String collapse = r1.getAttribute(COLLAPSE_TAG).toString();
            int cmp = collapse.compareTo(r2.getAttribute(COLLAPSE_TAG).toString());
            if (cmp != 0) {
                return cmp;
            }
            final String r1Uncollapsed = r1.getAttribute(UNCOLLAPSED_TAG).toString();
            final String r2Uncollapsed = r2.getAttribute(UNCOLLAPSED_TAG).toString();
            cmp = r1Uncollapsed.compareTo(r2Uncollapsed);
            if (cmp == 0) {
                return 0;
            }
            if (collapse.equals(r1Uncollapsed)) {
                return -1;
            }
            if (collapse.equals(r2Uncollapsed)) {
                return 1;
            }
            return r1Uncollapsed.compareTo(r2Uncollapsed);
        }
    };

    private final Comparator<SAMRecord> PARENT_CHILD_COMPARATOR = (r1, r2) -> {
        int cmp = GROUPING_COMPARATOR.compare(r1, r2);
        if (cmp != 0) {
            return cmp;
        }
        return r1.getReadName().compareTo(r2.getReadName());
    };

}
