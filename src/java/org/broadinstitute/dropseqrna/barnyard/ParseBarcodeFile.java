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
package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.util.IOUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.dropseqrna.TranscriptomeException;

import picard.util.BasicInputParser;

public class ParseBarcodeFile {

	/**
	 * A single column file with no header, contains the cell barcodes.
	 * @param input
	 * @return
	 */
	public static List<String > readCellBarcodeFile (final File input) {
		IOUtil.assertFileIsReadable(input);
		List<String> result = new ArrayList<String>();
		BasicInputParser parser = new BasicInputParser(false, 1, input);
		while(parser.hasNext()) {
			String [] line =parser.next();
			result.add(line[0]);
		}
		parser.close();
		return (result);
	}

	/**
	 * Parses a tab delimited file with 2 columns and a header.
	 * The header is CLUSTER BARCODE.  The header is required.
	 * The two columns are a cluster identifier and a barcode per line of the file.
	 * @param input
	 * @return A map from each CLUSTER identifier to the set of barcodes related to it.
	 */
	public static Map<String, Set<String>> readCellClusterFile (final File input) {
		IOUtil.assertFileIsReadable(input);
		Map<String, Set<String>> result = new HashMap<String, Set<String>>();
		BasicInputParser parser = new BasicInputParser(false, 2, input);
		if (parser.hasNext()) {
			String [] header = parser.next();
			if (!header[0].equals("CLUSTER") || !header[1].equals("BARCODE"))
				throw new TranscriptomeException("The expected header for a cluster file is CLUSTER    BARCODE (tab separaterd)");
		}

		while(parser.hasNext()) {
			String [] line =parser.next();
			String cluster = line[0];
			String barcode = line[1];
			Set<String> c = result.get(cluster);
			if (c==null)
				c=new HashSet<String>();
			c.add(barcode);
			result.put(cluster, c);
		}
		return result;
	}
}
