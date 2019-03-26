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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.vcf.VCFFileReader;

import java.util.Collections;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

/**
 * Utility for pulling sequence dictionaries out of a pair of objects, and testing the intersection of sequences names
 * in the 2 dictionaries.  Objects can be any of the following:
 * - {@link SAMSequenceDictionary}
 * - {@link SAMFileHeader}
 * - {@link SamReader}
 * - {@link VCFFileReader}
 * - {@link IntervalList}
 */
public class SequenceDictionaryIntersection {
    private enum ObjectType {DICT, BAM_HEADER, BAM, VCF, INTERVAL_LIST}

    private static ObjectType getObjectType(final Object o, String description) {
        if (o == null) {
            throw new IllegalArgumentException("o==null");
        } else if (o instanceof SAMSequenceDictionary) {
            return ObjectType.DICT;
        } else if (o instanceof SAMFileHeader) {
            return ObjectType.BAM_HEADER;
        } else if (o instanceof SamReader) {
            return ObjectType.BAM;
        } else if (o instanceof VCFFileReader) {
            return ObjectType.VCF;
        } else if (o instanceof IntervalList) {
            return ObjectType.INTERVAL_LIST;
        } else {
            if (description == null) {
                description = o.getClass().getName();
            }
            throw new IllegalArgumentException(String.format("Cannot convert %s into sequence dictionary", description));
        }
    }

    private static SAMSequenceDictionary getSequenceDictionary(final Object o, final ObjectType objectType) {
        switch(objectType) {
            case DICT:
                return (SAMSequenceDictionary)o;
            case BAM_HEADER:
                return ((SAMFileHeader) o).getSequenceDictionary();
            case BAM:
                return ((SamReader)o).getFileHeader().getSequenceDictionary();
            case VCF:
                return ((VCFFileReader)o).getFileHeader().getSequenceDictionary();
            case INTERVAL_LIST:
                return ((IntervalList)o).getHeader().getSequenceDictionary();
        }
        throw new IllegalStateException("unpossible");
    }



    private final String description1;
    private final String description2;
    private final Set<String> sequences1;
    private final Set<String> sequences2;
    private final Set<String> intersection;
    private final Set<String> sequencesOnlyIn1;
    private final Set<String> sequencesOnlyIn2;

    /**
     * Compare sequence dictionaries of 2 objects containing SDs
     * @param o1
     * @param description1 For including in messages
     * @param o2
     * @param description2 For including in messages
     */
    public SequenceDictionaryIntersection(final Object o1, final String description1, final Object o2, final String description2) {
        final ObjectType type1 = getObjectType(o1, description1);
        final ObjectType type2 = getObjectType(o2, description2);
        this.description1 = (description1 != null? description1: type1.name());
        this.description2 = (description2 != null? description2: type2.name());;
        sequences1 = Collections.unmodifiableSet(getSequenceSet(getSequenceDictionary(o1, type1)));
        sequences2 = Collections.unmodifiableSet(getSequenceSet(getSequenceDictionary(o2, type2)));
        final TreeSet<String> intersection = new TreeSet<>(sequences1);
        intersection.retainAll(sequences2);
        this.intersection = Collections.unmodifiableSet(intersection);
        final TreeSet<String> onlyIn1 = new TreeSet<>(sequences1);
        onlyIn1.removeAll(sequences2);
        this.sequencesOnlyIn1 = Collections.unmodifiableSet(onlyIn1);
        final TreeSet<String> onlyIn2 = new TreeSet<>(sequences2);
        onlyIn1.removeAll(sequences1);
        this.sequencesOnlyIn2 = Collections.unmodifiableSet(onlyIn2);
    }

    /**
     * Compare sequence dictionaries of 2 objects containing SDs.
     * In messages, object are referred to by their type, e.g. BAM, VCF, etc.
     * @param o1
     * @param o2
     */
    public SequenceDictionaryIntersection(final Object o1, final Object o2) {
        this(o1, null, o2, null);
    }

    private Set<String> getSequenceSet(final SAMSequenceDictionary sd) {
        return new TreeSet<>(sd.getSequences().stream().map((r) -> r.getSequenceName()).collect(Collectors.toSet()));
    }

    /**
     * @return user-friendly name for first object
     */
    public String getDescription1() {
        return description1;
    }

    /**
     * @return user-friendly name for second object
     */
    public String getDescription2() {
        return description2;
    }

    /**
     * @return sequences from first object, lexically sorted.
     */
    public Set<String> getSequences1() {
        return sequences1;
    }

    /**
     * @return sequences from second object, lexically sorted.
     */
    public Set<String> getSequences2() {
        return sequences2;
    }

    /**
     * @return interesection of sequences, lexically sorted.
     */
    public Set<String> getIntersection() {
        return intersection;
    }

    /**
     * @return sequences from first object that are not in second object, lexically sorted.
     */
    public Set<String> getSequencesOnlyIn1() {
        return sequencesOnlyIn1;
    }

    /**
     * @return sequences from second object that are not in first object, lexically sorted.
     */
    public Set<String> getSequencesOnlyIn2() {
        return sequencesOnlyIn2;
    }

    private String asString(Set<String> set) {
        return String.join(", ", set);
    }

    /**
     * @return comma-separated string of first-object sequences
     */
    public String getSequences1AsString() {
        return asString(sequences1);
    }

    /**
     * @return comma-separated string of second-object sequences
     */
    public String getSequences2AsString() {
        return asString(sequences2);
    }

    /**
     * @return comma-separated string of shared sequences
     */
    public String getIntersectionAsString() {
        return asString(intersection);
    }


    /**
     * @param verbose more info if true
     * @return String suitable for logging about intersection of the 2 sequence dictionaries
     */
    public String message(boolean verbose) {
        final StringBuilder sb = new StringBuilder();
        if (verbose) {
            sb.append(String.format("%s has %d contigs.\n", description1, sequences1.size()));
            sb.append(String.format("%s has %d contigs.\n", description2, sequences2.size()));
        }
        if (intersection.isEmpty()) {
            sb.append("No contigs in common!\n");
        } else {
            sb.append(String.format("Contigs in common: %s.\n", getIntersectionAsString()));
        }
        if (verbose || intersection.isEmpty()) {
            sb.append(String.format("%s contigs: %s\n", description1, getSequences1AsString()));
            sb.append(String.format("%s contigs: %s\n", description2, getSequences2AsString()));
        }
        return sb.toString();
    }
}
