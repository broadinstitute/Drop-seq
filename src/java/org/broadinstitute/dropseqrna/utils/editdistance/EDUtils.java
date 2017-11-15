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
package org.broadinstitute.dropseqrna.utils.editdistance;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

public class EDUtils {

	// private final Log log = Log.getInstance(EDUtils.class);

	private static EDUtils singletonInstance;


	public static EDUtils getInstance() {
		if (null == singletonInstance)
			singletonInstance = new EDUtils();
		return singletonInstance;
	}

	private EDUtils() {
	}

	public Set<String> getStringsWithinEditDistanceWithIndel(final String baseString,
			final List<String> comparisonStrings, final int editDistance) {
		Set<String> result = new HashSet<String>();
		for (String b : comparisonStrings) {
			int ed = LevenshteinDistance.getIndelSlidingWindowEditDistance(baseString, b);
			if (ed <= editDistance)
				result.add(b);
		}
		return (result);
	}

	public Set<String> getStringsWithinEditDistance(final String baseString,
			final List<String> comparisonStrings, final int editDistance) {
		Set<String> result = new HashSet<String>();
		for (String b : comparisonStrings) {
			int ed = HammingDistance.getHammingDistance(baseString, b);
			if (ed <= editDistance)
				result.add(b);
		}
		return (result);
	}


	public Set<String> getStringsWithinHammingDistance(final String baseString,
			final List<String> comparisonStrings, final int editDistance) {
		Set<String> result = new HashSet<String>();
		for (String b : comparisonStrings) {
			boolean flag = HammingDistance.greaterThanHammingDistance(baseString, b, editDistance);
			if (flag==false)
				result.add(b);
		}
		return (result);
	}


	/**
	 *
	 * @param aFile The input file to read.  2 columns, the number of observations of that barcode followed by the barcode sequence. Tab seperated.
	 * @return a list of Barcodes with counts.
	 */
	public static ObjectCounter <String> readBarCodeFile(final File aFile) {
		ObjectCounter <String> result = new ObjectCounter<String>();

		try {
			BufferedReader input = new BufferedReader(new FileReader(aFile));
			try {
				String line = null; // not declared within while loop
				while ((line = input.readLine()) != null) {
					line=line.trim();
					String[] strLine = line.split("\t");
					int count = Integer.parseInt(strLine[0]);
					String barcode = strLine[1].toUpperCase();
					result.incrementByCount(barcode, count);
				}
			} finally {
				input.close();
			}
		} catch (IOException ex) {
			throw new TranscriptomeException("Could not read file: "
					+ aFile.toString());
		}

		return (result);
	}


}
