package org.broadinstitute.dropseqrna.barnyard;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorStrategy;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import picard.annotation.LocusFunction;
import picard.cmdline.CommandLineProgram;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public abstract class GeneFunctionCommandLineBase extends CommandLineProgram {

    public static final String DEFAULT_GENE_NAME_TAG = "gn";
    public static final String DEFAULT_GENE_STRAND_TAG = "gs";
    public static final String DEFAULT_GENE_FUNCTION_TAG = "gf";
    public static final String DEFAULT_FUNCTION_TAG = "XF";
    public static final StrandStrategy DEFAULT_STRAND_STRATEGY = StrandStrategy.SENSE;

    public static final FunctionalDataProcessorStrategy DEFAULT_FUNCTIONAL_STRATEGY = FunctionalDataProcessorStrategy.DROPSEQ;
    public static final List<LocusFunction> DEFAULT_LOCUS_FUNCTION_LIST = Collections.unmodifiableList(new ArrayList<>(Arrays.asList(LocusFunction.CODING, LocusFunction.UTR)));
    @Argument(doc="Gene Name tag.  Takes on the gene name this read overlaps (if any)")
	public String GENE_NAME_TAG= DEFAULT_GENE_NAME_TAG;

	@Argument(doc="Gene Strand tag.  For a given gene name <GENE_NAME_TAG>, this is the strand of the gene.")
	public String GENE_STRAND_TAG= DEFAULT_GENE_STRAND_TAG;

	@Argument(doc="Gene Function tag.  For a given gene name <GENE_NAME_TAG>, this is the function of the gene at this read's position: UTR/CODING/INTRONIC/...")
	public String GENE_FUNCTION_TAG= DEFAULT_GENE_FUNCTION_TAG;

    @Argument(doc="The strand strategy decides which reads will be used by analysis.  The SENSE strategy requires the read and annotation to have the same strand.  "
    		+ "The ANTISENSE strategy requires the read and annotation to be on opposite strands.  The BOTH strategy is permissive, and allows the read to be on either strand.")
    public StrandStrategy STRAND_STRATEGY= DEFAULT_STRAND_STRATEGY;

    @Argument(doc="A list of functional annotations that reads need to be completely contained by to be considered for analysis.")
    public List<LocusFunction> LOCUS_FUNCTION_LIST= new ArrayList<>(DEFAULT_LOCUS_FUNCTION_LIST);

    @Argument(doc="A strategy for interpreting functional annotations.  DropSeq is the default strategy.  STARSOLO strategy priority is very similar to DropSeq, except" +
            "in cases where a read overlaps both an intron on the sense strand and a coding region on the antisense strand.  In these cases, DropSeq " +
            "favors the intronic interpretation, while STARSolo interprets this as a technical artifact and labels the read as coming from the antisense coding gene, " +
            "and the read does not contribute to the expression counts matrix.")
    public FunctionalDataProcessorStrategy FUNCTIONAL_STRATEGY=DEFAULT_FUNCTIONAL_STRATEGY;
}
