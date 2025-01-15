package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Iterator;
import java.util.regex.Pattern;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.dropseqrna.utils.TransformingIterator;

public class BAMTagCleanupIterator extends TransformingIterator<SAMRecord, SAMRecord> {

    private final String tag;
    private final String prefixToRemove;
    private final String prefixToAdd;
    private final String suffixToRemove;
    private final String suffixToAdd;
    private final Pattern patternToRemove;

    public BAMTagCleanupIterator(Iterator<SAMRecord> underlyingIterator, String tag, final String prefixToRemove, final String prefixToAdd, String suffixToRemove, String suffixToAdd, Pattern patternToRemove) {
        super(underlyingIterator);
        this.tag = tag;
        this.prefixToRemove = prefixToRemove;
        this.prefixToAdd = prefixToAdd;
        this.suffixToRemove = suffixToRemove;
        this.suffixToAdd = suffixToAdd;
        this.patternToRemove = patternToRemove;
    }

    @Override
    public SAMRecord next() {
        SAMRecord r = this.underlyingIterator.next();
        String tagValue = r.getStringAttribute(tag);

        if (tagValue != null) {
            // Remove prefix if it exists
            if (prefixToRemove != null && tagValue.startsWith(prefixToRemove)) {
                tagValue = tagValue.substring(prefixToRemove.length());
            }

            // Add prefix
            if (prefixToAdd != null) {
                tagValue = prefixToAdd + tagValue;
            }

            // Remove suffix if it exists
            if (suffixToRemove != null && tagValue.endsWith(suffixToRemove)) {
                tagValue = tagValue.substring(0, tagValue.length() - suffixToRemove.length());
            }

            // Add suffix
            if (suffixToAdd != null) {
                tagValue = tagValue + suffixToAdd;
            }

            // Remove pattern if it exists
            if (patternToRemove != null) {
                tagValue = patternToRemove.matcher(tagValue).replaceAll("");
            }

            r.setAttribute(tag, tagValue);
        }

        return r;
    }

    public static class Builder {
        private final Iterator<SAMRecord> underlyingIterator;
        private String tag;
        private String prefixToRemove;
        private String prefixToAdd;
        private String suffixToRemove;
        private String suffixToAdd;
        private Pattern patternToRemove;

        public Builder(Iterator<SAMRecord> underlyingIterator) {
            this.underlyingIterator = underlyingIterator;
        }

        public Builder tag(String tag) {
            this.tag = tag;
            return this;
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

        public BAMTagCleanupIterator build() {
            return new BAMTagCleanupIterator(this.underlyingIterator, this.tag, this.prefixToRemove, this.prefixToAdd, this.suffixToRemove, this.suffixToAdd, this.patternToRemove);
        }
    }
}