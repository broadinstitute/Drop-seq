package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.io.File;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;

import htsjdk.samtools.util.IOUtil;

public class GdacAlleleFrequencyWriter {

	private final DecimalFormat ratioFormat = new DecimalFormat("0.000");
	private final PrintStream out;
	
	
	public GdacAlleleFrequencyWriter (File f) {
		this.out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(f));
	}
	
	public void writeHeader () {		
		List<String> header = Arrays.asList("chromosome", "position", "ref_allele", "alt_allele", "ref_reads", "alt_reads", "ref_umi", "alt_umi", "maf_reads", "maf_umi");
		String h = StringUtils.join(header, "\t");
		out.println(h);
	}
	
	public void writeLine (final DigitalAlleleCounts dac) {
		writeLine(new GdacAlleleFrequency(dac));
	}
	
	public void writeLine(final GdacAlleleFrequency f) {
				
		List<String> line = new ArrayList<>();
		line.addAll(Arrays.asList(f.getSnpInterval().getContig(), Integer.toString(f.getSnpInterval().getStart()), f.getRefAllele()+"", f.getAltAllele()+""));
		
		// if there are 0 reads on the ref/alt allele, get rid of this.
		if (f.getRefReadCount()+f.getAltReadCount()==0)
			return;
		
		line.add(Integer.toString(f.getRefReadCount()));
		line.add(Integer.toString(f.getAltReadCount()));
		
				
		// if there are 0 umis on the ref/alt allele, get rid of this.
		if (f.getRefUmiCount()+f.getAltUmiCount()==0)
			return;
				
		line.add(Integer.toString(f.getRefUmiCount()));
		line.add(Integer.toString(f.getAltUmiCount()));
		
		String readMaf=getRatioStringFormatted(f.getRefReadCount(), f.getAltReadCount(), ratioFormat);		
		String umiMaf=getRatioStringFormatted(f.getRefUmiCount(), f.getAltUmiCount(), ratioFormat);
		line.add(readMaf);
		line.add(umiMaf);
		
		
		String h = StringUtils.join(line, "\t");
		out.println(h);	
		
	}
	
	private String getRatioStringFormatted (int ref, int alt, DecimalFormat ratioFormat) {
		if (ref+alt==0)
			return ("NA");
		double maf = (double) alt/ (ref+alt);
		// if (Double.isInfinite(umiMaf)) 
		return (ratioFormat.format(maf));							
	}
	
	public void close () {
		this.out.close();
	}
	
	
}
