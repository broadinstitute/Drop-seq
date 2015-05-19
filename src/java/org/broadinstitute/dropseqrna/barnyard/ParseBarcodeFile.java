package org.broadinstitute.dropseqrna.barnyard;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import picard.util.BasicInputParser;

public class ParseBarcodeFile {

	public static List<String > readCellBarcodeFile (File input) {
		List<String> result = new ArrayList<String>();
		BasicInputParser parser = new BasicInputParser(false, 1, input);
		while(parser.hasNext()) {
			String [] line =parser.next();
			result.add(line[0]);
		}
		
		return (result);
	}
}
