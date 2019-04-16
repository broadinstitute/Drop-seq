/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.BarcodeListRetrieval;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@CommandLineProgramProperties(summary = "Convert from a cell barcode tag to a sample group", oneLineSummary = "Convert from a cell barcode tag to a sample group", omitFromCommandLine = false, programGroup = DropSeq.class)
public class ConvertTagToReadGroup extends CommandLineProgram {

	private static final Log log = Log.getInstance(ConvertTagToReadGroup.class);
	private ProgressLogger progress = new ProgressLogger(log);

	public final String USAGE = "USAGE: ConvertTagToReadGroup.  For a read tag (like a cell barcode) find the top set of tags as defined by reads that have at least some map quality.  "
			+ "Tag reads to have a read group with the same name as the tag, to take advantage of Picard programs that aggregate data by read group.";

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.  Must be coordinate sorted. (???)")
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file of per-cell exonic/intronic/genic/intergenic/rRNA levels.")
	public File OUTPUT;

	@Argument(doc = "The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
	public String CELL_BARCODE_TAG = "XC";

	@Argument(doc = "The map quality of the read to be included for determining which cells will be measured.")
	public Integer READ_MQ = 10;

	@Argument(doc = "Number of cells that you think are in the library. The top NUM_CORE_BARCODES will be reported in the output.", optional=false, mutex={"CELL_BC_FILE"})
	public Integer NUM_CORE_BARCODES = null;

	@Argument(doc="Override NUM_CORE_BARCODES, and process reads that have the cell barcodes in this file instead.  The file has 1 column with no header.", optional=false, mutex={"NUM_CORE_BARCODES"})
	public File CELL_BC_FILE=null;

	@Argument(doc = "The same name for the experiment.  If given, will be concatonated onto the cell barcode in the fashion sample:cell.", optional = false)
	public String SAMPLE_NAME;

	@Argument(shortName="LB", doc="The library name to place into the LB attribute in the read group header", optional=true)
    public String LIBRARY_NAME;

    @Argument(shortName="PU", doc="The platform unit (often run_barcode.lane) to insert into the read group header", optional=true)
    public String PLATFORM_UNIT;

    @Argument(shortName="PL", doc="The platform type (e.g. illumina, solid) to insert into the read group header", optional=true)
    public String PLATFORM;

	@Override
	protected int doWork() {

		IOUtil.assertFileIsWritable(OUTPUT);
		IOUtil.assertFileIsWritable(INPUT);

		SamReader in = SamReaderFactory.makeDefault().open(INPUT);
		final SAMFileHeader inHeader = in.getFileHeader();
		SAMReadGroupRecord readGroupTemplate = getReadGroupTemplate(inHeader);

		Set<String> cellBarcodes = getCellBarcodes ();

		// get read groups
		List<SAMReadGroupRecord> rg = getReadGroups(cellBarcodes, readGroupTemplate);
		final SAMFileHeader outHeader = inHeader.clone();
		outHeader.setReadGroups(rg);
		final SAMFileWriter outWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(outHeader,outHeader.getSortOrder() == inHeader.getSortOrder(), OUTPUT);

		for (SAMRecord r : in) {
			String cellBarcode = r.getStringAttribute(this.CELL_BARCODE_TAG);
			if (cellBarcodes.contains(cellBarcode)) {
				r.setAttribute(SAMTag.RG.name(), cellBarcode);
				outWriter.addAlignment(r);
			}
			progress.record(r);
		}

		CloserUtil.close(in);
		outWriter.close();

		log.info("Tagging finished.");
		return (0);
	}

	private Set<String> getCellBarcodes () {
		BarcodeListRetrieval u = new BarcodeListRetrieval();
		if (this.NUM_CORE_BARCODES!=null) {
			Set<String> cellBarcodes = new HashSet<>(u.getListCellBarcodesByReadCount(this.INPUT, this.CELL_BARCODE_TAG, this.READ_MQ, null, this.NUM_CORE_BARCODES));
			return (cellBarcodes);
		}

		// this must be set.
		List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(this.CELL_BC_FILE);
		log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
		return (new HashSet<>(cellBarcodes));
	}

	/**
	 * If all the read groups have the same template with the same library, platform, platformUnit
	 * you can return this group.  Otherwise returns null.
	 * @param inHeader
	 * @return
	 */
	private SAMReadGroupRecord getReadGroupTemplate(final SAMFileHeader inHeader) {
		List<SAMReadGroupRecord> recs = inHeader.getReadGroups();
		String library = recs.get(0).getLibrary();
		String platform = recs.get(0).getPlatform();
		String platformUmit = recs.get(0).getPlatformUnit();
		boolean multiple=false;
		// checks to see that the params are the same if there are multiple read groups.
		for (SAMReadGroupRecord gr : recs) {
			if (!gr.getLibrary().equals(library)) multiple=true;
			if (!gr.getPlatform().equals(platform)) multiple=true;
			if (!gr.getPlatformUnit().equals(platformUmit)) multiple=true;
		}

		if (inHeader.getReadGroups().size()>1 && multiple)
			log.warn("There is more than one read group in this BAM.  Using the first group, overriding PLATFORM/PLATFORM_UNIT/LIBRARY_NAME with provided values!");
		else
			log.info("Auto-determined Project, Platform, Platform Unit.");

		SAMReadGroupRecord gr = inHeader.getReadGroups().get(0);
		if (this.PLATFORM!=null ) {
			log.info("Overriding Platform with value " + this.PLATFORM);
			gr.setPlatform(this.PLATFORM);
		}
		if (this.PLATFORM_UNIT!=null) {
			log.info("Overriding Platform Unit with value " + this.PLATFORM_UNIT);
			gr.setPlatformUnit(this.PLATFORM_UNIT);
		}
		if (this.LIBRARY_NAME!=null) {
			log.info("Overriding Library Name with value " + this.LIBRARY_NAME);
			gr.setLibrary(this.LIBRARY_NAME);
		}
		return gr;
	}


	public List<SAMReadGroupRecord> getReadGroups(final Set<String> cellBarcodes, final SAMReadGroupRecord template) {
		List<SAMReadGroupRecord> g = new ArrayList<>(
				cellBarcodes.size());
		for (String id : cellBarcodes) {
			SAMReadGroupRecord rg = new SAMReadGroupRecord(id);
			rg.setLibrary(template.getLibrary());
			rg.setPlatform(template.getPlatform());
			String sampleID = this.SAMPLE_NAME + ":" + id;
			rg.setSample(sampleID);
			rg.setPlatformUnit(template.getPlatformUnit());
			g.add(rg);
		}
		return (g);

	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new ConvertTagToReadGroup().instanceMain(args));
	}
}
