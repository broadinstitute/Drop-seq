/*
 * MIT License
 *
 * Copyright 2026 Broad Institute
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
package org.broadinstitute.dropseqrna.utils.readiterators.umisort;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import static org.broadinstitute.dropseqrna.utils.readiterators.umisort.UMISortCodecIO.*;

/**
 * Where the reads of one UMI aligned: contig, strand, the smallest alignment start and largest alignment end, the
 * number of reads with at least one gap (N in the CIGAR), and for each gap the number of reads that have it.
 * <p>
 * Built from one read, then {@link #merge}d with the summaries of the other reads of the same UMI.  Every value is
 * either a sum, a min or a max, so the result does not depend on the order in which the reads are merged.  The number
 * of reads is not held here: it is the UMI's read count.
 */
public final class ReadIntervalSummary {

    /** A gap between two alignment blocks: the 1-based reference position before the gap and the one after it. */
    public record Gap(int start, int end) implements Comparable<Gap> {
        @Override
        public int compareTo(final Gap o) {
            final int c = Integer.compare(start, o.start);
            return c != 0 ? c : Integer.compare(end, o.end);
        }
    }

    private final String contig;
    private boolean positiveStrand;
    private int positionMin;
    private int positionMax;
    private int numGappedReads;
    // Null until a gapped read is seen; most reads have no gaps.  Sorted by position.
    private TreeMap<Gap, Integer> gaps;

    private ReadIntervalSummary(final String contig, final boolean positiveStrand, final int positionMin,
                                final int positionMax, final int numGappedReads, final TreeMap<Gap, Integer> gaps) {
        this.contig = contig;
        this.positiveStrand = positiveStrand;
        this.positionMin = positionMin;
        this.positionMax = positionMax;
        this.numGappedReads = numGappedReads;
        this.gaps = gaps;
    }

    /** @return the summary of one aligned read. */
    public static ReadIntervalSummary fromRead(final SAMRecord r) {
        final List<Gap> readGaps = getGaps(r);
        TreeMap<Gap, Integer> gaps = null;
        if (!readGaps.isEmpty()) {
            gaps = new TreeMap<>();
            for (final Gap g : readGaps)
                gaps.merge(g, 1, Integer::sum);
        }
        return new ReadIntervalSummary(r.getContig(), !r.getReadNegativeStrandFlag(), r.getAlignmentStart(),
                r.getAlignmentEnd(), readGaps.isEmpty() ? 0 : 1, gaps);
    }

    /**
     * Adds the reads summarized by other to this summary, and returns this.  Reads of one UMI share a gene, so they
     * share a contig; the contig of this summary is kept.  The strand is positive if any read is on the positive strand.
     */
    public ReadIntervalSummary merge(final ReadIntervalSummary other) {
        positiveStrand |= other.positiveStrand;
        positionMin = Math.min(positionMin, other.positionMin);
        positionMax = Math.max(positionMax, other.positionMax);
        numGappedReads = Math.addExact(numGappedReads, other.numGappedReads);
        if (other.gaps != null) {
            if (gaps == null)
                gaps = new TreeMap<>();
            for (final Map.Entry<Gap, Integer> e : other.gaps.entrySet())
                gaps.merge(e.getKey(), e.getValue(), Math::addExact);
        }
        return this;
    }

    public String getContig() {
        return contig;
    }

    public boolean isPositiveStrand() {
        return positiveStrand;
    }

    /** @return the smallest 1-based alignment start of the reads. */
    public int getPositionMin() {
        return positionMin;
    }

    /** @return the largest 1-based inclusive alignment end of the reads. */
    public int getPositionMax() {
        return positionMax;
    }

    public int getNumGappedReads() {
        return numGappedReads;
    }

    /** @return each gap and the number of reads with it, ordered by position. */
    public Map<Gap, Integer> getGaps() {
        return gaps == null ? Collections.emptyMap() : Collections.unmodifiableMap(gaps);
    }

    void encode(final DataOutputStream out) throws IOException {
        writeString(out, contig);
        out.writeBoolean(positiveStrand);
        writeVarInt(out, positionMin);
        writeVarInt(out, positionMax);
        writeVarInt(out, numGappedReads);
        writeVarInt(out, gaps == null ? 0 : gaps.size());
        if (gaps != null)
            for (final Map.Entry<Gap, Integer> e : gaps.entrySet()) {
                writeVarInt(out, e.getKey().start());
                writeVarInt(out, e.getKey().end());
                writeVarInt(out, e.getValue());
            }
    }

    static ReadIntervalSummary decode(final DataInputStream in) throws IOException {
        final String contig = readString(in);
        final boolean positiveStrand = in.readBoolean();
        final int positionMin = readVarInt(in);
        final int positionMax = readVarInt(in);
        final int numGappedReads = readVarInt(in);
        final int numGaps = readVarInt(in);
        TreeMap<Gap, Integer> gaps = null;
        if (numGaps > 0) {
            gaps = new TreeMap<>();
            for (int i = 0; i < numGaps; i++)
                gaps.put(new Gap(readVarInt(in), readVarInt(in)), readVarInt(in));
        }
        return new ReadIntervalSummary(contig, positiveStrand, positionMin, positionMax, numGappedReads, gaps);
    }

    @Override
    public String toString() {
        return "ReadIntervalSummary{" + contig + ":" + positionMin + "-" + positionMax + ", positiveStrand=" +
                positiveStrand + ", numGappedReads=" + numGappedReads + ", gaps=" + getGaps() + "}";
    }

    /** The gaps of one read, one per N (skipped region) in the CIGAR. */
    private static List<Gap> getGaps(final SAMRecord r) {
        final List<Integer> alignmentBlocksPrecedingGaps = getIndicesOfAlignmentBlocksPrecedingGaps(r);
        if (alignmentBlocksPrecedingGaps.isEmpty())
            return Collections.emptyList();
        final List<Gap> ret = new ArrayList<>(alignmentBlocksPrecedingGaps.size());
        final List<AlignmentBlock> alignmentBlocks = r.getAlignmentBlocks();
        for (int alignmentBlockIndex : alignmentBlocksPrecedingGaps) {
            final AlignmentBlock precedingAlignmentBlock = alignmentBlocks.get(alignmentBlockIndex);
            if (alignmentBlockIndex == alignmentBlocks.size() - 1) {
                throw new IllegalStateException("Strange CIGAR: " + r.getCigarString() +
                        " for read " + r.getReadName());
            }
            ret.add(new Gap(precedingAlignmentBlock.getReferenceStart() + precedingAlignmentBlock.getLength() - 1,
                    alignmentBlocks.get(alignmentBlockIndex + 1).getReferenceStart()));
        }
        return ret;
    }

    private static List<Integer> getIndicesOfAlignmentBlocksPrecedingGaps(final SAMRecord r) {
        ArrayList<Integer> ret = null;
        int alignmentBlockIndex = -1;
        for (final CigarElement ce : r.getCigar().getCigarElements()) {
            if (ce.getOperator() == CigarOperator.MATCH_OR_MISMATCH) {
                ++alignmentBlockIndex;
            } else if (ce.getOperator() == CigarOperator.SKIPPED_REGION) {
                if (ret == null) {
                    ret = new ArrayList<>();
                }
                ret.add(alignmentBlockIndex);
            }
        }
        return ret == null ? Collections.emptyList() : ret;
    }
}
