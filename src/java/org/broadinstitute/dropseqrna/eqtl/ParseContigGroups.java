package org.broadinstitute.dropseqrna.eqtl;

import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.yaml.snakeyaml.Yaml;

import htsjdk.samtools.util.IOUtil;

public class ParseContigGroups {

	/**
	 * Build a map from a contig group name to all the contigs contained in that
	 * group A contig may be contained in multiple groups.
	 * 
	 * @param f
	 * @return
	 */
	public static Map<String, Set<String>> getContigGroupMapSimple(File f) {

		final Yaml yaml = new Yaml();
		Map<String, Object> yamlMap = (Map<String, Object>) yaml.load(IOUtil.openFileForReading(f));

		Map<String, Set<String>> result = new HashMap<String, Set<String>>();

		for (Object key : yamlMap.keySet()) {
			String contig = (String) key;
			Set<String> groups = new HashSet<String>();
			;

			Object value = yamlMap.get(key);
			if (value instanceof String) {
				String v = (String) value;
				groups.add(v);
			}
			if (value instanceof Collection) {
				Collection v = (Collection) value;
				for (Object vv : v) {
					groups.add((String) vv);
				}
			}

			for (String group : groups) {
				Set<String> contigSet = result.get(group);
				if (contigSet == null) {
					contigSet = new HashSet<String>();
					result.put(group, contigSet);
				}
				contigSet.add(contig);
			}

		}

		return result;
	}

	/*
	public static Map<String, Set<String>> getContigGroupMap(File f) {
		Constructor constructor = new Constructor();

		final Yaml yaml = new Yaml(constructor);

		// keys are contigs, values are either single elements or lists of groups.
		Iterable<Object> iter = yaml.loadAll(IOUtil.openFileForReading(f));
		for (Object o : iter) {
			ContigAnnotation ca = (ContigAnnotation) o;
			
		}
		
		Map<String, Set<String>> result = null;
		return (result);
	}

	public class ContigAnnotationList {
		List<ContigAnnotation> list;

		public List<ContigAnnotation> getList() {
			return list;
		}

		public void setList(List<ContigAnnotation> list) {
			this.list = list;
		}
		
	}
	public class ContigAnnotation {

		public ContigAnnotation (String contig, List<String> annotations) {
			this.contig=contig;
			this.annotations=annotations;
		}
		
		public ContigAnnotation() {
			this.annotations=new ArrayList<String>();
		}
		
		public String getContig() {
			return contig;
		}

		public void setContig(String contig) {
			this.contig = contig;
		}

		public List<String> getAnnotations() {
			return annotations;
		}

		public void setAnnotations(List<String> annotations) {
			this.annotations = annotations;
		}

		private String contig;
		private List<String> annotations;
	}
	*/
}
