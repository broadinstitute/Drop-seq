/*
 * MIT License
 *
 * Copyright 2024 Broad Institute
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

package org.broadinstitute.dropseqrna.eqtl;

import java.io.File;
import java.util.*;

public record NonNumericCovariate(File file, String attribute, String donor, String value) {
    public static String generateNonNumericMessage(List<NonNumericCovariate> nonNumericCovariates) {
        Map<File, Map<String, List<NonNumericCovariate>>> groupedCovariates = new HashMap<>();

        // Group the non-numeric covariates by file, then by attribute
        for (NonNumericCovariate covariate : nonNumericCovariates) {
            groupedCovariates
                    .computeIfAbsent(covariate.file(), k -> new LinkedHashMap<>())
                    .computeIfAbsent(covariate.attribute(), k -> new ArrayList<>())
                    .add(covariate);
        }

        return generateNonNumericMessage(groupedCovariates);
    }

    @SuppressWarnings("StringConcatenationInLoop")
    private static String generateNonNumericMessage(
            Map<File, Map<String, List<NonNumericCovariate>>> groupedCovariates
    ) {
        String message = "Non-numeric covariates found:\n";
        for (final Map.Entry<File, Map<String, List<NonNumericCovariate>>> fileEntry : groupedCovariates.entrySet()) {
            message += "File: '" + fileEntry.getKey().getAbsolutePath() + "':\n";
            for (final Map.Entry<String, List<NonNumericCovariate>> attributeEntry : fileEntry.getValue().entrySet()) {
                message += " - Attribute: '" + attributeEntry.getKey() + "':\n";
                for (final NonNumericCovariate covariate : attributeEntry.getValue()) {
                    message += "    - Donor: '" + covariate.donor() + "', Value: '" + covariate.value() + "'\n";
                }
            }
        }
        return message;
    }
}
