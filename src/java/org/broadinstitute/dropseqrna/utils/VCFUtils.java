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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;

import java.io.File;
import java.nio.file.Path;
import java.util.Collection;

public class VCFUtils {
	private static final Log log = Log.getInstance(VCFUtils.class);

	public static boolean hasIndex(final File vcfFile) {
		return hasIndex(vcfFile.toPath());
	}

	public static boolean hasIndex(final Path vcfFile) {
		// validate that there's an index for the VCF file.
		VCFFileReader vcfReader=null;
		try {
			vcfReader = new VCFFileReader(vcfFile, true);
		} catch (TribbleException te) {
			log.error("Must supply an index file for [" + FileUtils.toAbsoluteString(vcfFile)+"]");
			return false;
		} finally {
			CloserUtil.close(vcfReader);
		}
		return true;
	}

	/**
	 * Does the VCF file use genotype qualities?
	 * @param vcfReader
	 * @return
	 */
	public static boolean GQInHeader (final VCFFileReader vcfReader) {
		Collection<VCFFormatHeaderLine> r = vcfReader.getFileHeader().getFormatHeaderLines();
		for (VCFFormatHeaderLine l : r) {
			String id = l.getID();
			if (id.equals("GQ")) return true;
		}
		return false;
	}
}
