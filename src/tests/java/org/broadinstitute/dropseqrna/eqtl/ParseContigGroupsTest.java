package org.broadinstitute.dropseqrna.eqtl;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.dropseqrna.eqtl.ParseContigGroups;
import org.testng.Assert;
import org.testng.annotations.Test;

public class ParseContigGroupsTest {
	private final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/eqtl/");
	private final File contigGroupFile=new File(TEST_DATA_DIR, "GRCh38_GRCm38.contig_groups.yaml");
	
	@Test
	public void testGetContigGroupMapSimple() {
		Map<String,Set<String>> map = ParseContigGroups.getContigGroupMapSimple(contigGroupFile);
		Map<String,Set<String>> expected = getExpectedResults();
		
		Assert.assertNotNull(map);
		Assert.assertEquals(expected,  map);				
	}
	
	/*
	@Test
	public void testGetContigGroupMap() {
		Map<String,Set<String>> map = ParseContigGroups.getContigGroupMap(test);
		Map<String,Set<String>> expected = getExpectedResults();
		
		Assert.assertNotNull(map);
		Assert.assertEquals(expected,  map);				
	}
	*/
	
	private Map<String, Set<String>> getExpectedResults() {
		Map<String, Set<String>> result = new HashMap<>();
		result.put("X", new HashSet<String>(Arrays.asList("chrX")));
		result.put("Y", new HashSet<String>(Arrays.asList("chrY")));
		result.put("MT", new HashSet<String>(Arrays.asList("chrM")));
		result.put("non-autosome", new HashSet<String>(Arrays.asList("chrX","chrY", "chrM")));
		HashSet<String> autosomes = new HashSet<String>();
		
		for (int i=1; i<=22; i++) {
			autosomes.add("chr"+i);
		}		
		result.put("autosome", autosomes);
		return result;
	}
}
