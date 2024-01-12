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

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.metrics.umisharing.ParentEditDistanceMatcher;
import org.broadinstitute.dropseqrna.metrics.umisharing.ParentEditDistanceMatcher.TagValues;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.ProgressLoggingIterator;
import org.broadinstitute.dropseqrna.utils.SortingIteratorFactory;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionIteratorWrapper;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.RuntimeIOException;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(summary = "Calculate UMI Sharing in a flexible way between parent / child entities.  "
		+ "After running a collapse algorithm that merges entities (for example cell barcodes), "
		+ "calculate the UMIs that are shared between the parent entity and each of the children that has been collapsed. "
		+ "For example, if cell barcodes have been collapsed into a new tag, compare the parent cell barcode (the new "
		+ "cell barcode all results were collapsed into) to each of the child cell barcodes (all of the cells that were collapsed into the main one). "
		+ "UMI Sharing is calcuated in a flexible way by gathering the unique COUNT_TAG values for each entity, then determining what fraction of "
		+ "those counts in the child can be explain as coming from the parent counts at an edit distance.  For example, for a child cell barcode, gather counts of "
		+ "the gene, gene strand and molecular barcode, and allow an edit distance of 0 for the gene and gene strand, and 1 for the molecular barcode, "
		+ "then determine the fraction of UMIs in the child barcode that could have arisen from the parent barcode. "
		+ "\n\nIMPLEMENTATION DETAILS: Reads a BAM and analyzes groups of reads that share the same value for " +
        "COLLAPSE_TAG.  Within such a group, reads for which UNCOLLAPSED_TAG value is the same as COLLAPSE_TAG value " +
        "are put into the parent subgroup.  The remaining reads are grouped by UNCOLLAPSED_TAG value into child " +
        "subgroups.  For each read, the values for the COUNT_TAG form a tuple.  For each subgroup, the set of unique " +
        "tuples are accumulated.  For each child subgroup, the size of the intersection of the tuple set with the " +
        "parent tuple set is computed.  A metrics file is produced for each child subgroup.",
        oneLineSummary = "Computes UMI sharing between uncollapsed and collapsed sets of reads.",
        programGroup = DropSeq.class)
public class ComputeUMISharing
        extends GeneFunctionCommandLineBase {

    private static final Log log = Log.getInstance(ComputeUMISharing.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input SAM or BAM file to analyze.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The metrics file to be written.")
    public File OUTPUT;

    @Argument(doc="This is the tag that has been collapsed by some other process, for example a cell barcode tag after a cell barcode collapse process has been performed on the data.  The tag used to groups reads.")
    public String COLLAPSE_TAG;

    @Argument(doc="The tag before collapse took place.  For example, the cell barcode tag before a cell barcode collapse process was performed on the data.  The tag used to create subgroups of reads within a group, and compares the UMI Sharing of each uncollapsed tag to it's parent collapsed tag.")
    public String UNCOLLAPSED_TAG;

    @Argument(minElements = 1, doc="One or more tags used to define how a UMI should be counted.  For a cell barcode, this would be the gene, gene strand, and UMI tags. The tag(s) whose values are used to create the set of unique tuples in each subgroup.")
    public List<String> COUNT_TAG;

    @Argument(doc="The edit distance of the COUNT_TAGS that definies sharing across COLLAPSE and UNCOLLAPSED TAGs. Since there can be multiple COUNT_TAG elements, this allows you to define the edit distance threshold for each COUNT_TAG individually in the same order as the COUNT_TAGs.  "
    		+ "For example, if you defined a UMI by the gene, gene strand, and UMI TAGS, you might want edit distances of 0,0,1, so that on the molecular barcode could have an inexact match. "+ 
            "If there are fewer EDIT_DISTANCEs than COUNT_TAGs, the remainder are assumed to have EDIT_DISTANCE=0.")
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
        parentEditDistanceMatcher = new ParentEditDistanceMatcher(this.COUNT_TAG, this.EDIT_DISTANCE, this.FIND_INDELS, this.NUM_THREADS);

        SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
        final ProgressLogger progressLogger = new ProgressLogger(log, 1000000, "Sorting");
        Iterator<SAMRecord> iter = reader.iterator();
        if (LOCUS_FUNCTION_LIST.size() > 0) {
            iter = new GeneFunctionIteratorWrapper(iter, this.GENE_NAME_TAG,
                    this.GENE_STRAND_TAG, this.GENE_FUNCTION_TAG, false, this.STRAND_STRATEGY,
                    this.LOCUS_FUNCTION_LIST, this.FUNCTIONAL_STRATEGY);
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
                parentTuples = parentEditDistanceMatcher.getValues(parentSubgroup);                
            } else {                
                final List<SAMRecord> childSubgroup = subgroupIterator.next();
                final Set<TagValues> childTuples = parentEditDistanceMatcher.getValues(childSubgroup);                
                final UmiSharingMetrics metrics = new UmiSharingMetrics();
                metrics.PARENT = parentSubgroup.get(0).getAttribute(COLLAPSE_TAG).toString();
                metrics.CHILD = childSubgroup.get(0).getAttribute(UNCOLLAPSED_TAG).toString();
                metrics.NUM_PARENT = parentTuples.size();
                metrics.NUM_CHILD = childTuples.size();
                metrics.NUM_SHARED = parentEditDistanceMatcher.computeNumShared(parentTuples, childTuples);
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
