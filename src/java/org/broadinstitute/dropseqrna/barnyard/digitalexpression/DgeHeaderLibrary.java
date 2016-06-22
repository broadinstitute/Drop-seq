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

import java.io.File;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

public class DgeHeaderLibrary {
    private File input;
    private File inputDge;
    private File reference;
    final private String uei;
    private String prefix;
    private final LinkedHashMap<String, String> otherTags = new LinkedHashMap<>();

    public DgeHeaderLibrary(final String uei) {
        this.uei = uei;
    }

    public File getInput() {
        return input;
    }

    public void setInput(File input) {
        this.input = input;
    }

    public File getInputDge() {
        return inputDge;
    }

    public void setInputDge(File inputDge) {
        this.inputDge = inputDge;
    }

    public File getReference() {
        return reference;
    }

    public void setReference(File reference) {
        this.reference = reference;
    }

    public String getUei() {
        return uei;
    }

    public String getPrefix() {
        return prefix;
    }

    public void setPrefix(String prefix) {
        this.prefix = prefix;
    }

    public void setTag(final String key, final String value) {
        otherTags.put(key, value);
    }

    public String getTag(final String key) {
        return otherTags.get(key);
    }

    public Iterator<Map.Entry<String, String>> iterateOtherTags() {
        return otherTags.entrySet().iterator();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        DgeHeaderLibrary that = (DgeHeaderLibrary) o;

        if (input != null ? !input.equals(that.input) : that.input != null) return false;
        if (inputDge != null ? !inputDge.equals(that.inputDge) : that.inputDge != null) return false;
        if (reference != null ? !reference.equals(that.reference) : that.reference != null) return false;
        if (!uei.equals(that.uei)) return false;
        if (prefix != null ? !prefix.equals(that.prefix) : that.prefix != null) return false;
        return otherTags.equals(that.otherTags);

    }

    @Override
    public int hashCode() {
        int result = input != null ? input.hashCode() : 0;
        result = 31 * result + (inputDge != null ? inputDge.hashCode() : 0);
        result = 31 * result + (reference != null ? reference.hashCode() : 0);
        result = 31 * result + uei.hashCode();
        result = 31 * result + (prefix != null ? prefix.hashCode() : 0);
        result = 31 * result + otherTags.hashCode();
        return result;
    }
}
