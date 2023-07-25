package org.broadinstitute.dropseqrna.utils;

import static org.testng.Assert.assertEquals;

import java.io.File;
import java.util.List;

import org.testng.annotations.Test;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;

public class FilteredReadsMetricTest {
	
	
	static final File METRICS_FILE_1=new File ("testdata/org/broadinstitute/dropseq/utils/13.organism_filter_metrics");
	static final File METRICS_FILE_2=new File ("testdata/org/broadinstitute/dropseq/utils/14.organism_filter_metrics");
	
	@Test
	public void testMerge() {
		FilteredReadsMetric m1 = new FilteredReadsMetric();
		m1.READS_ACCEPTED++;
		m1.READS_REJECTED++;
		
		FilteredReadsMetric m2 = new FilteredReadsMetric();
		m2.READS_ACCEPTED++;
		
		m1.merge(m2);
		
		assertEquals(m1.READS_ACCEPTED, 2);
		assertEquals(m1.READS_REJECTED, 1);
		
		
	}
	
	@Test
	public void testMergeTwoFiles () {
		MetricsFile<FilteredReadsMetric, String> m1 = new MetricsFile<FilteredReadsMetric, String>();
		m1.read(IOUtil.openFileForBufferedReading(METRICS_FILE_1));
		
		MetricsFile<FilteredReadsMetric, String> m2 = new MetricsFile<FilteredReadsMetric, String>();
		m2.read(IOUtil.openFileForBufferedReading(METRICS_FILE_2));
		
		FilteredReadsMetric frm1= m1.getMetrics().get(0);
		FilteredReadsMetric frm2= m2.getMetrics().get(0);
		frm1.merge(frm2);
		
		assertEquals(frm1.READS_ACCEPTED, 74570376);
		assertEquals(frm1.READS_REJECTED, 2728571);
	}
}
