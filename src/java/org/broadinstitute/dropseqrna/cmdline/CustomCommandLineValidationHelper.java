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
package org.broadinstitute.dropseqrna.cmdline;

import org.apache.commons.lang3.ArrayUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * Manages annoying null check and conversion between List and primitive array.
 */
public class CustomCommandLineValidationHelper {
    /**
     * @param superErrors may be null
     * @param thisErrors may be null
     * @return combined content of superErrors and thisErrors, or null if there are none.
     */
    public static String[] makeValue(
            final String[] superErrors,
            List<String> thisErrors) {
        if (thisErrors == null || thisErrors.isEmpty()) {
            return superErrors;
        } else {
            if (superErrors != null) {
                // Make a copy because argument may not be mutable.
                thisErrors = new ArrayList<>(thisErrors);
                for (final String msg: superErrors) {
                    thisErrors.add(msg);
                }
            }
        }
        return thisErrors.toArray(new String[thisErrors.size()]);
    }

    public static String[] makeValue(
            final String[] superErrors,
            final String[] thisErrors) {
        if (superErrors == null) {
            return thisErrors;
        } else if (thisErrors == null) {
            return superErrors;
        } else {
            return (String[])ArrayUtils.addAll(superErrors, thisErrors);
        }
    }
}
