package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import htsjdk.samtools.util.Log;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionProcessor;
import picard.annotation.LocusFunction;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class ValidateAnnotations {
    private static final Log log = Log.getInstance(ValidateAnnotations.class);
    private final ClassifyDropSeqFunctionalData dropseqClassifier;
    private final ClassifyStarSoloFunctionalData starClassifier;

    public ValidateAnnotations(GeneFunctionProcessor gfp) {
        dropseqClassifier = new ClassifyDropSeqFunctionalData(gfp);
        starClassifier = new ClassifyStarSoloFunctionalData(gfp);

    }

    /**
     * For each read, validate that the annotated function is the same for the two sets of tags.
     * StarSolo produces a priority score but for some priority scores like antisense does include the gene name.
     * DropSeq's list of FunctionalData needs to be simplified to a single result - in the case of an ambiguous
     * result (multiple genes with the same priority) the result can be simplified to an INTERGENIC result to match STAR.
     * @param starsoloFD
     * @param fdMap
     * @return A list of invalid read names.  An empty list if all reads pass validation.
     */
    public Map<String, ValidationStatus> validate (final Map<String, FunctionalData> starsoloFD, final Map<String, List<FunctionalData>> fdMap, boolean verbose) {
        Map<String, ValidationStatus> validationResult = new HashMap<>();

        for (String readName: starsoloFD.keySet()) {
//            if (readName.equals("LH00118:69:22CNFNLT3:3:2164:21055:18833"))
//                log.info("STOP");
            ValidationStatus validationStatus = validate(starsoloFD.get(readName), fdMap.get(readName), verbose);
            validationResult.put(readName, validationStatus);
            if (!validationStatus.isValid()) {
                if (verbose) log.info("Read fails validation [" + readName+"]");
            }
        }
        return validationResult;
    }

    public ValidationStatus validate (final FunctionalData starSolo, final List<FunctionalData> dropseq, boolean verbose) {
        GeneWithFunction star = this.starClassifier.convert(starSolo);
        GeneWithFunction ds = this.dropseqClassifier.convert(dropseq);
        ValidationStatus result = new ValidationStatus(star, ds);
        return result;
    }

}
