package org.broadinstitute.dropseqrna.censusseq;

import java.io.File;
import java.util.List;

import org.broadinstitute.dropseqrna.censusseq.SNPSampleRecord;
import org.testng.Assert;
import org.testng.annotations.Test;

public class SNPSampleRecordTest {

	File ROLL_CALL_VERBOSE_FILE = new File ("testdata/org/broadinstitute/dropseq/censusseq/RollCall.verbose.txt.gz");

	@Test
	public void testParseFile() {
		List<SNPSampleRecord> result = SNPSampleRecord.parseFile(ROLL_CALL_VERBOSE_FILE);
		Assert.assertNotNull(result);
		// the 5th record:
		//   CHR      POS REF_BASE ALT REF ALT_COUNT          DONOR GENOTYPE
		//   22 16120030        T   C   1         0 WA7_P33_140529      het

		SNPSampleRecord r= result.get(4);
		testRecordsSame("22", 16120030, 'T', 'C', 1, 0, "WA7_P33_140529", "het", r);

		// the 50th record:
		// 22 16379041        C   T   1         0 HUES74_P7_150201      het

		r= result.get(49);
		testRecordsSame("22", 16379041, 'C', 'T', 1, 0, "HUES74_P7_150201", "het", r);

		// the 1000th record:
		// 22 18222350        G   A   1         0 Genea42_P19_150107      het

		r= result.get(499);
		testRecordsSame("22", 18222350, 'G', 'A', 1, 0, "Genea42_P19_150107", "het", r);

	}

	private void testRecordsSame (final String chr, final int pos, final char refBase, final char altBase, final int refCount, final int altCount, final String donorName, final String genotype, final SNPSampleRecord r) {
		Assert.assertEquals(chr, r.getInterval().getContig());
		Assert.assertEquals(pos, r.getInterval().getStart());
		Assert.assertEquals(r.getRefBase(), refBase);
		Assert.assertSame(altBase, r.getAltBase());
		Assert.assertEquals(refCount, r.getRefCount());
		Assert.assertEquals(altCount, r.getAltCount());
		Assert.assertEquals(donorName, r.getDonorName());
		Assert.assertEquals(genotype, r.getGenotype());
	}
}
