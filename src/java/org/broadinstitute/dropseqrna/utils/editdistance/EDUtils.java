package org.broadinstitute.dropseqrna.utils.editdistance;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.TranscriptomeException;

public class EDUtils {
	
	// private final Log log = Log.getInstance(EDUtils.class);
	
	private static EDUtils singletonInstance;
	
	
	public static EDUtils getInstance() {
		if (null == singletonInstance) {
			singletonInstance = new EDUtils();
		}
		return singletonInstance;
	}

	private EDUtils() {
	}
	
	public Set<String> getStringsWithinEditDistanceWithIndel(String baseString,
			List<String> comparisonStrings, int editDistance) {
		Set<String> result = new HashSet<String>();
		for (String b : comparisonStrings) {
			int ed = LevenshteinDistance.getIndelSlidingWindowEditDistance(baseString, b);
			if (ed <= editDistance)
				result.add(b);
		}
		return (result);
	}
	
	public Set<String> getStringsWithinEditDistance(String baseString,
			List<String> comparisonStrings, int editDistance) {
		Set<String> result = new HashSet<String>();
		for (String b : comparisonStrings) {
			int ed = HammingDistance.getHammingDistance(baseString, b);
			if (ed <= editDistance)
				result.add(b);
		}
		return (result);
	}
	
	
	public Set<String> getStringsWithinHammingDistance(String baseString,
			List<String> comparisonStrings, int editDistance) {
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
	public List<BarcodeWithCount> readBarCodeFile(File aFile) {
		List<BarcodeWithCount> result = new ArrayList<BarcodeWithCount>();

		try {
			BufferedReader input = new BufferedReader(new FileReader(aFile));
			try {
				String line = null; // not declared within while loop
				while ((line = input.readLine()) != null) {
					line=line.trim();
					String[] strLine = line.split("\t");
					int count = Integer.parseInt(strLine[0]);
					BarcodeWithCount bc = new BarcodeWithCount(strLine[1].toUpperCase(),
							count);
					result.add(bc);
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

	class BarcodeWithCountComparator implements Comparator<BarcodeWithCount> {
		@Override
		public int compare(BarcodeWithCount arg0, BarcodeWithCount arg1) {
		// TODO Auto-generated method stub
		int count =  arg1.getCount() - arg0.getCount();
		if (count != 0)
			return (count);
		return (0);
		}
	}
	
	public List<String> getBarcodes(List<BarcodeWithCount> barcodes) {
		List<String> result = new ArrayList<String>(barcodes.size());
		for (BarcodeWithCount b : barcodes) {
			result.add(b.getBarcode());
		}
		return (result);
	}
	
	public List<BarcodeWithCount> sortBC (List<BarcodeWithCount> barcodes) {
		Collections.sort(barcodes, new BarcodeWithCountComparator());
		return (barcodes);
	}
	
	
}