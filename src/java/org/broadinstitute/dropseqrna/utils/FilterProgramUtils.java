/*
 * MIT License
 *
 * Copyright 2019 Broad Institute
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

import htsjdk.samtools.util.Log;

public class FilterProgramUtils {
    public static void reportAndCheckFilterResults(
            final String elementType, final long elementsAccepted, final long elementsRejected,
            final Double passingThreshold, final Log log) {
        final long totalElements = elementsAccepted + elementsRejected;
        log.info(String.format("Total %d %s processed.  %d %s accepted; %d %s rejected.",
                totalElements, elementType, elementsAccepted, elementType, elementsRejected, elementType));
        if (passingThreshold != null) {
            if (passingThreshold >= 1) {
                if (elementsAccepted < passingThreshold) {
                    throw new RuntimeException(String.format("Fewer than %d %s passed filters",
                            passingThreshold.intValue(), elementType));
                }
            } else if (elementsAccepted/((double)totalElements) < passingThreshold) {
                throw new RuntimeException(String.format("A smaller fraction than %f %s passed filters",
                        passingThreshold, elementType));
            }
        }
        
    }
}
