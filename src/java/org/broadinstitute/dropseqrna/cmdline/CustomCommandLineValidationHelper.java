/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2016 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.cmdline;

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
}
