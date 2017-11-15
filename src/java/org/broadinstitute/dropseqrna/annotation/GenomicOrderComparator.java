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
package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.lang.builder.CompareToBuilder;

import java.util.Comparator;

public class GenomicOrderComparator implements Comparator <GTFRecord> {
	SAMSequenceDictionary refDict;
	
	public GenomicOrderComparator (SAMSequenceDictionary refDict) {
		this.refDict = refDict;
	}
	
	public int compare(GTFRecord g1, GTFRecord g2) {
		int i1 = refDict.getSequenceIndex(g1.getChromosome());
		int i2 = refDict.getSequenceIndex(g2.getChromosome());
		return new CompareToBuilder().append(i1, i2)
				.append(g1.getStart(), g2.getStart())
				.append(g1.getEnd(), g2.getEnd())
                .append(g1.getTranscriptType(), g2.getTranscriptType())
				.toComparison();
	}
}
