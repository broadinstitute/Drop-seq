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

/**
 * The minimal information {@link org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator} needs from a read once it has passed the SAMRecord filters:
 * cell, gene and UMI, plus (only when low map quality reads are kept) the read name, the read strand and
 * the gene strand / function tags that are needed to resolve reads that map to more than one place.
 * Compare with a SAMRecord: this is spilled to disk by {@link UMIReadCodec} in a few dozen bytes
 * instead of a full BAM record.
 * <p>
 * Instances are immutable. Nullable fields: readName, geneStrand, geneFunction.
 */
public final class UMIRead {
    public final String readName;
    public final boolean negativeStrand;
    public final String cell;
    public final String gene;
    public final String geneStrand;
    public final String geneFunction;
    public final String umi;

    public UMIRead(final String readName, final boolean negativeStrand, final String cell, final String gene,
                   final String geneStrand, final String geneFunction, final String umi) {
        this.readName = readName;
        this.negativeStrand = negativeStrand;
        this.cell = cell;
        this.gene = gene;
        this.geneStrand = geneStrand;
        this.geneFunction = geneFunction;
        this.umi = umi;
    }

    /** @return a copy of this read with the gene, strand and function replaced. A null function keeps the current one. */
    public UMIRead withGeneAssignment(final String gene, final String geneStrand, final String geneFunction) {
        return new UMIRead(readName, negativeStrand, cell, gene, geneStrand,
                geneFunction == null ? this.geneFunction : geneFunction, umi);
    }

    @Override
    public String toString() {
        return "UMIRead{" + readName + ", cell=" + cell + ", gene=" + gene + ", strand=" + geneStrand +
                ", function=" + geneFunction + ", umi=" + umi + ", negativeStrand=" + negativeStrand + "}";
    }
}
