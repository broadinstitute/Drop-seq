/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
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
