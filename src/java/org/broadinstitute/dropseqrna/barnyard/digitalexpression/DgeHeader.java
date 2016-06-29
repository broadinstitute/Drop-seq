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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import java.util.Iterator;
import java.util.LinkedHashMap;

public class DgeHeader {
    public enum ExpressionFormat {raw,log10, log10_normalized, unknown}
    public static String CURRENT_VERSION = "1.1";

    private String version = CURRENT_VERSION;
    private ExpressionFormat expressionFormat = ExpressionFormat.raw;

    private LinkedHashMap<String, DgeHeaderLibrary> libraries = new LinkedHashMap<>();

    public DgeHeader() {
    }

    public String getVersion() {
        return version;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public ExpressionFormat getExpressionFormat() {
        return expressionFormat;
    }

    public void setExpressionFormat(ExpressionFormat expressionFormat) {
        this.expressionFormat = expressionFormat;
    }

    /**
     *
     * @return previous library with same UEI, or null if there wasn't one.
     */
    public DgeHeaderLibrary addLibrary(final DgeHeaderLibrary library) {
        return libraries.put(library.getUei(), library);
    }

    public DgeHeaderLibrary getLibrary(final String uei) {
        return libraries.get(uei);
    }

    public DgeHeaderLibrary getLibrary(int i) {
        if (libraries.size() <= i) {
            throw new IllegalArgumentException("Requested library " + i + " but there are only " + libraries.size() +
            " libraries in DgeHeader");
        }
        final Iterator<DgeHeaderLibrary> it = libraries.values().iterator();
        while (i-- > 0) {it.next();}
        return it.next();
    }

    public DgeHeaderLibrary removeLibrary(final String uei) {
        return libraries.remove(uei);
    }

    public DgeHeaderLibrary removeLibrary(int i) {
        if (libraries.size() <= i) {
            throw new IllegalArgumentException("Requested library " + i + " but there are only " + libraries.size() +
                    " libraries in DgeHeader");
        }
        final Iterator<DgeHeaderLibrary> it = libraries.values().iterator();
        while (i-- > 0) {it.next();}
        final DgeHeaderLibrary ret = it.next();
        it.remove();
        return ret;
    }

    public int getNumLibraries() { return libraries.size(); }

    public Iterator<String> iterateLibraries() {
        return libraries.keySet().iterator();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        DgeHeader dgeHeader = (DgeHeader) o;

        if (version != null ? !version.equals(dgeHeader.version) : dgeHeader.version != null) return false;
        if (expressionFormat != dgeHeader.expressionFormat) return false;
        return libraries.equals(dgeHeader.libraries);

    }

    @Override
    public int hashCode() {
        int result = version != null ? version.hashCode() : 0;
        result = 31 * result + (expressionFormat != null ? expressionFormat.hashCode() : 0);
        result = 31 * result + libraries.hashCode();
        return result;
    }
}