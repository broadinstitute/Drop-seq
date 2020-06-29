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
package org.broadinstitute.dropseqrna.censusseq;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import picard.util.BasicInputParser;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class SNPSampleRecord {

	private final Interval interval;
	private final char refBase;
	private final char altBase;
	private final int refCount;
	private final int altCount;
	private final String donorName;
	private final String genotype;
	private final List<Byte> qualities;
	
	public SNPSampleRecord (final Interval interval, final char refBase, final char altBase, final int refCount, final int altCount, final String donorName, final String genotype) {
		this(interval, refBase, altBase, refCount, altCount, donorName, genotype, null);
	}
	
	public SNPSampleRecord (final Interval interval, final char refBase, final char altBase, final int refCount, final int altCount, final String donorName, final String genotype, final List<Byte> qualities) {
		this.interval=interval;
		this.refBase=refBase;
		this.altBase=altBase;
		this.refCount=refCount;
		this.altCount=altCount;
		this.donorName=donorName.intern();
		this.genotype=genotype.intern();
		this.qualities=qualities;
	}

	public Interval getInterval() {
		return interval;
	}

	public int getRefCount() {
		return refCount;
	}

	public int getAltCount() {
		return altCount;
	}

	public String getDonorName() {
		return donorName;
	}

	public String getGenotype() {
		return genotype;
	}

	public char getRefBase() {
		return refBase;
	}

	public char getAltBase() {
		return altBase;
	}

	public static List<SNPSampleRecord> parseFile (final File f) {
		IOUtil.assertFileIsReadable(f);
		List<SNPSampleRecord> result = new ArrayList<SNPSampleRecord>();
		BasicInputParser parser = new BasicInputParser(false, 8, f);
		// get rid of header.
		if (parser.hasNext())
			parser.next();

		while(parser.hasNext()) {
			String [] line =parser.next();
			// String [] header = {"CHR", "POS", "REF_BASE", "ALT", "REF", "ALT_COUNT", "DONOR", "GENOTYPE"};
			String chr = line[0];
			int pos = Integer.parseInt(line[1]);
			char refBase = line[2].charAt(0);
			char altBase = line[3].charAt(0);
			int refCount = Integer.parseInt(line[4]);
			int altCount = Integer.parseInt(line[5]);
			String donor = line[6];
			String genotype = line[7];
			SNPSampleRecord r = new SNPSampleRecord(new Interval(chr, pos, pos), refBase, altBase, refCount, altCount, donor, genotype, null);
			result.add(r);
		}
		parser.close();
		return (result);
	}

	public List<Byte> getQualities() {
		return qualities;
	}



}
