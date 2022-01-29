package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.GdacAlleleFrequency;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.GdacAlleleFrequencyReader;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.GdacAlleleFrequencyWriter;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.Interval;


public class GdacAlleleFrequencyReaderTest {
	
	private final File ALLELE_FREQ_FILE= new File ("testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_5_cell_2_snp_testdata_retagged.dac.allele_freq.txt");
	
	@Test
	public void testParseMap() {
		
		GdacAlleleFrequencyReader reader = new GdacAlleleFrequencyReader(this.ALLELE_FREQ_FILE);
		
		Map<Interval, Double> result = reader.getUmiAlleleFrequencyMap();
		Assert.assertEquals(2, result.keySet().size());
		
		Assert.assertEquals(0.446, result.get(new Interval("HUMAN_1", 76227022,76227022)), 0.0001);
		Assert.assertEquals(0.329, result.get(new Interval("HUMAN_1", 150199123,150199123)), 0.0001);
				
	}
	
	
	// construct a few GdacAlleleFrequency objects, write them to disk, read them back and compare the two objects for equality.
	// tests both the reader and writer.
	@Test
	public void testRoundTrip () throws IOException {
		List<GdacAlleleFrequency> list = new ArrayList<>();
		list.add(new GdacAlleleFrequency(new Interval("chr1", 1, 1), 'A', 'T', 10, 5, 4, 2));
		list.add(new GdacAlleleFrequency(new Interval("chr1", 3, 3), 'C', 'T', 20, 7, 2, 1));
		list.add(new GdacAlleleFrequency(new Interval("chr2", 5, 5), 'G', 'T', 5, 14, 2, 6));
		
		
		File out = File.createTempFile("GdacAlleleFrequencyReaderTest.", ".roundtrip.txt");
		out.deleteOnExit();
		
		// write output.
		GdacAlleleFrequencyWriter w = new GdacAlleleFrequencyWriter(out);
		w.writeHeader();
		list.stream().forEach(x -> w.writeLine(x));
		w.close();
		
		// read from output.
		GdacAlleleFrequencyReader r = new GdacAlleleFrequencyReader(out);
		List<GdacAlleleFrequency> input = r.parseFile();
		
		// assert both lists are equal.
		Assert.assertEquals(list, input);
		
	}
	
	@Test
	public void testMerge() {
		GdacAlleleFrequency t1 = new GdacAlleleFrequency(new Interval("chr1", 1, 1), 'A', 'T', 10, 5, 4, 2);
		GdacAlleleFrequency t2 = new GdacAlleleFrequency(new Interval("chr1", 1, 1), 'A', 'T', 5, 2, 1, 1);
		
		GdacAlleleFrequency expected = new GdacAlleleFrequency(new Interval("chr1", 1, 1), 'A', 'T', 15, 7, 5, 3);
		GdacAlleleFrequency m = t1.merge(t2);
		
		Assert.assertEquals(m, expected);
		
		double rr = t1.getReadRatio();
		Assert.assertEquals(rr, 0.3333, 0.001);
		double ur = t1.getUMIRatio();
		Assert.assertEquals(ur, 0.3333, 0.001);
		
		int cmp = t1.compareTo(t2);
		Assert.assertEquals(cmp, 0);		
		
	}
	
	
}
