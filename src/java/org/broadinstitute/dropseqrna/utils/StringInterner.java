/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2015 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils;

import java.util.HashMap;
import java.util.Map;

/**
 * Replacement for String.intern() that allows for all the interned strings to be GCed when this cache is GCed.
 */
public class StringInterner {

    private Map<String, String> cache = new HashMap<>();

    public String intern(final String str) {
        String cachedValue = cache.get(str);
        if (cachedValue != null) {
            return cachedValue;
        } else {
            cache.put(str, str);
            return str;
        }
    }
}
