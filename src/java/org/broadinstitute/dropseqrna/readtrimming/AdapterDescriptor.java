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
package org.broadinstitute.dropseqrna.readtrimming;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTagUtil;
import htsjdk.samtools.util.SequenceUtil;

import java.util.ArrayList;
import java.util.List;

public class AdapterDescriptor {
    /** RCed molecular barcode, RCed cell barcode, plus fixed adapter sequence */
    public static String DEFAULT_ADAPTER= "~XM~XCACGTACTCTGCGTTGCTACCACTG";

    interface AdapterElement {
        String getSequence(SAMRecord rec);
    }

    static class FixedAdapterElement implements AdapterElement {
        private final String sequence;

        FixedAdapterElement(String sequence) { this.sequence = sequence; }

        @Override
        public String getSequence(SAMRecord rec) { return sequence; }

        @Override
        public String toString() {
            return sequence;
        }
    }

    static class TagAdapterElement implements AdapterElement {
        private final String stringTag;
        private final short binaryTag;
        private final boolean reverseComplement;

        TagAdapterElement(String tag, boolean reverseComplement) {
            this.stringTag = tag;
            this.binaryTag = SAMTagUtil.getSingleton().makeBinaryTag(tag);
            this.reverseComplement = reverseComplement;
        }

        @Override
        public String getSequence(SAMRecord rec) {
            return SequenceUtil.reverseComplement((String) rec.getAttribute(binaryTag));
        }

        @Override
        public String toString() {
            return (reverseComplement? RC_TAG_CHAR: TAG_CHAR) + stringTag;
        }
    }

    private static char TAG_CHAR = '^';
    private static char RC_TAG_CHAR = '~';
    private static int TAG_LENGTH = 2; // All tags are 2 chars long

    private static AdapterElement[] parseAdapterDescriptor(final String descriptor) {
        final List<AdapterElement> elements = new ArrayList<>();
        int i = 0;
        while (i < descriptor.length()) {
            // Slurp a contiguous chunk of fixed bases
            final StringBuilder builder = new StringBuilder();
            for(; i < descriptor.length() && descriptor.charAt(i) != TAG_CHAR && descriptor.charAt(i) != RC_TAG_CHAR; ++i) {
                builder.append(descriptor.charAt(i));
            }
            if (builder.length() > 0) {
                elements.add(new FixedAdapterElement(builder.toString()));
            }

            if (i == descriptor.length()) { break; }

            final boolean reverseComplement = (descriptor.charAt(++i) == RC_TAG_CHAR);
            elements.add(new TagAdapterElement(descriptor.substring(i, i+TAG_LENGTH), reverseComplement));
            i += TAG_LENGTH;
        }
        return elements.toArray(new AdapterElement[elements.size()]);
    }

    private final AdapterElement[] elements;

    public AdapterDescriptor(final String descriptor) {
        elements = parseAdapterDescriptor(descriptor);
    }

    public String getAdapterSequence(final SAMRecord rec) {
        final StringBuilder builder = new StringBuilder();
        for (final AdapterElement element: elements) {
            builder.append(element.getSequence(rec));
        }
        return builder.toString();

    }

    @Override
    public String toString() {
        final StringBuilder builder = new StringBuilder();
        for (final AdapterElement element: elements) {
            builder.append(element.toString());
        }
        return builder.toString();
    }
}
