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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;

public class DgeHeader {
    public enum ExpressionFormat {raw,log10, log10_normalized, unknown}
    public static String CURRENT_VERSION = "1.1";

    private String version = CURRENT_VERSION;
    private ExpressionFormat expressionFormat = ExpressionFormat.raw;

    /** key: prefix */
    private LinkedHashMap<String, DgeHeaderLibrary> libraries = new LinkedHashMap<>();
    private ArrayList<DgeHeaderCommand> commands = new ArrayList<>();

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
        return libraries.put(library.getPrefix(), library);
    }

    public DgeHeaderLibrary getLibrary(final String prefix) {
        return libraries.get(prefix);
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

    public DgeHeaderLibrary removeLibrary(final String prefix) {
        return libraries.remove(prefix);
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

    public void addCommand(final DgeHeaderCommand command) {
        commands.add(command);
    }

    /**
     * If a CLP is invoked in a testing context, there may be no command line
     * @param maybeNullCommandLine If null, no DgeHeaderCommand record is added
     */
    public void addCommand(final String maybeNullCommandLine) {
        if (maybeNullCommandLine != null) {
            addCommand(new DgeHeaderCommand(maybeNullCommandLine));
        }
    }

    public DgeHeaderCommand getCommand(final int i) {
        return commands.get(i);
    }

    public DgeHeaderCommand removeCommand(final int i) {
        return commands.remove(i);
    }

    public int getNumCommands() { return commands.size(); }

    public Iterator<DgeHeaderCommand> iterateCommands() { return commands.iterator(); }

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
