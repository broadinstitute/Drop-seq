package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Iterator;
import java.util.regex.Pattern;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.dropseqrna.utils.TransformingIterator;

public class ReadNameCleanupIterator extends TransformingIterator<SAMRecord, SAMRecord> {

    private final String prefixToRemove;
    private final String prefixToAdd;
    private final String suffixToRemove;
    private final String suffixToAdd;
    private final Pattern patternToRemove;

    public ReadNameCleanupIterator(Iterator<SAMRecord> underlyingIterator, final String prefixToRemove, final String prefixToAdd, String suffixToRemove, String suffixToAdd, Pattern patternToRemove) {
        super(underlyingIterator);
        this.prefixToRemove = prefixToRemove;
        this.prefixToAdd = prefixToAdd;
        this.suffixToRemove = suffixToRemove;
        this.suffixToAdd = suffixToAdd;
        this.patternToRemove = patternToRemove;
    }

    @Override
    public SAMRecord next() {
        SAMRecord r = this.underlyingIterator.next();
        String readName = r.getReadName();

        // Remove prefix if it exists
        if (prefixToRemove != null && readName.startsWith(prefixToRemove)) {
            readName = readName.substring(prefixToRemove.length());
        }

        // Add prefix
        if (prefixToAdd != null) {
            readName = prefixToAdd + readName;
        }

        // Remove suffix if it exists
        if (suffixToRemove != null && readName.endsWith(suffixToRemove)) {
            readName = readName.substring(0, readName.length() - suffixToRemove.length());
        }

        // Add suffix
        if (suffixToAdd != null) {
            readName = readName + suffixToAdd;
        }

        // Remove pattern if it exists
        if (patternToRemove != null) {
            readName = patternToRemove.matcher(readName).replaceAll("");
        }

        r.setReadName(readName);
        return r;
    }

    public static class Builder {
        private final Iterator<SAMRecord> underlyingIterator;
        private String prefixToRemove;
        private String prefixToAdd;
        private String suffixToRemove;
        private String suffixToAdd;
        private Pattern patternToRemove;

        public Builder(Iterator<SAMRecord> underlyingIterator) {
            this.underlyingIterator = underlyingIterator;
        }

        public Builder prefixToRemove(String prefixToRemove) {
            this.prefixToRemove = prefixToRemove;
            return this;
        }

        public Builder prefixToAdd(String prefixToAdd) {
            this.prefixToAdd = prefixToAdd;
            return this;
        }

        public Builder suffixToRemove(String suffixToRemove) {
            this.suffixToRemove = suffixToRemove;
            return this;
        }

        public Builder suffixToAdd(String suffixToAdd) {
            this.suffixToAdd = suffixToAdd;
            return this;
        }

        public Builder patternToRemove(Pattern patternToRemove) {
            this.patternToRemove = patternToRemove;
            return this;
        }

        public ReadNameCleanupIterator build() {
            return new ReadNameCleanupIterator(this.underlyingIterator, this.prefixToRemove, this.prefixToAdd, this.suffixToRemove, this.suffixToAdd, this.patternToRemove);
        }
    }
}