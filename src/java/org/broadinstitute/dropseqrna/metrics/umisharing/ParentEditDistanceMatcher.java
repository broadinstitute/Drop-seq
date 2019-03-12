package org.broadinstitute.dropseqrna.metrics.umisharing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.broadinstitute.dropseqrna.utils.editdistance.MapBarcodesByEditDistance;

import htsjdk.samtools.SAMRecord;

public class ParentEditDistanceMatcher {
	
	private final int[] zeroEditDistanceIndices;
    private final int[] nonZeroEditDistanceIndices;    
    private final List<Integer> editDistance;
    private final List<String> countTag;
    private final boolean findIndels;
    private final MapBarcodesByEditDistance med;
    
    public ParentEditDistanceMatcher(List<String> countTag, List<Integer> editDistance, final boolean findIndels, int numThreads) {
        final ArrayList<Integer> zeroEDs = new ArrayList<>();
        final ArrayList<Integer> nonZeroEDs = new ArrayList<>();
        for (int i = 0; i < editDistance.size(); ++i) {
            if (editDistance.get(i) == 0) {
                zeroEDs.add(i);
            } else {
                nonZeroEDs.add(i);
            }
        }
        zeroEditDistanceIndices = zeroEDs.stream().mapToInt(i->i).toArray();
        nonZeroEditDistanceIndices = nonZeroEDs.stream().mapToInt(i->i).toArray();
        this.editDistance=editDistance;
        this.findIndels=findIndels;
        this.countTag=countTag;
        this.med = new MapBarcodesByEditDistance(false, numThreads, 0);
    }

    public TagValues getValues(final SAMRecord rec) {
        final TagValues ret = new TagValues(zeroEditDistanceIndices.length, nonZeroEditDistanceIndices.length);
        for (int i = 0; i < zeroEditDistanceIndices.length; ++i) {
            ret.zeroEditDistanceValues[i] = rec.getAttribute(countTag.get(zeroEditDistanceIndices[i]));
        }
        for (int i = 0; i < nonZeroEditDistanceIndices.length; ++i) {
            ret.nonZeroEditDistanceValues[i] = rec.getStringAttribute(countTag.get(nonZeroEditDistanceIndices[i]));
        }
        return ret;
    }
    
    public Set<TagValues> getValues (final Collection<SAMRecord> recs) {
    	return recs.stream().map(this::getValues).collect(Collectors.toSet());
    }
    
    public int computeNumShared(final Set<TagValues> parentTuples, final Set<TagValues> childTuples) {
        return (int)childTuples.stream().filter((childTuple) -> this.isShared(childTuple, parentTuples)).count();
    }

    public boolean isShared(final TagValues childTuple, final Set<TagValues> parentTuples) {
        if (nonZeroEditDistanceIndices.length == 0) {
            return parentTuples.contains(childTuple);
        }

        // Accumulate parents for which ED-zero values match the child.
        final List<TagValues> exactMatchParentTuples = parentTuples.stream().
                filter(childTuple::zeroEditDistanceEquals).collect(Collectors.toList());

        if (exactMatchParentTuples.isEmpty()) {
            return false;
        }

        if (nonZeroEditDistanceIndices.length == 1) {
            // Fast way for single tag
            final List<String> parentValues = exactMatchParentTuples.stream().map((tuple) -> tuple.nonZeroEditDistanceValues[0]).collect(Collectors.toList());
            return matchesWithinEditDistance(childTuple.nonZeroEditDistanceValues[0], parentValues, editDistance.get(nonZeroEditDistanceIndices[0]));
        }

        // Unfortunately, if more than one ED-nonzero tag, needs to test each individually.
        for (final TagValues parentTuple : exactMatchParentTuples) {
            boolean matched = true;
            for (int i = 0; i < nonZeroEditDistanceIndices.length; ++i) {
                if (!matchesWithinEditDistance(
                        childTuple.nonZeroEditDistanceValues[i],
                        parentTuple.nonZeroEditDistanceValues[i],
                        editDistance.get(nonZeroEditDistanceIndices[i]))) {
                    matched = false;
                    break;
                }
            }
            // Matched on every ED-nonzero tag
            if (matched) {
                return true;
            }
        }
        return false;
    }

    private boolean matchesWithinEditDistance(final String barcode, final List<String> comparisonBarcodes, final int editDistance) {
        return !med.processSingleBarcodeMultithreaded(barcode, comparisonBarcodes, findIndels, editDistance).isEmpty();
    }

    private boolean matchesWithinEditDistance(final String barcode, final String comparisonBarcode, final int editDistance) {
        return matchesWithinEditDistance(barcode, Collections.singletonList(comparisonBarcode), editDistance);
    }
    
    public static class TagValues {
        final Object[] zeroEditDistanceValues;
        final String[] nonZeroEditDistanceValues;

        TagValues(final int numZeroEditDistanceTags, final int numNonZeroEditDistanceTags) {
            zeroEditDistanceValues = new Object[numZeroEditDistanceTags];
            nonZeroEditDistanceValues = new String[numNonZeroEditDistanceTags];
        }

        boolean zeroEditDistanceEquals(final TagValues other) {
            return Arrays.equals(this.zeroEditDistanceValues, other.zeroEditDistanceValues);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            TagValues tagValues = (TagValues) o;
            return Arrays.equals(zeroEditDistanceValues, tagValues.zeroEditDistanceValues) &&
                    Arrays.equals(nonZeroEditDistanceValues, tagValues.nonZeroEditDistanceValues);
        }

        @Override
        public int hashCode() {
            int result = Arrays.hashCode(zeroEditDistanceValues);
            result = 31 * result + Arrays.hashCode(nonZeroEditDistanceValues);
            return result;
        }
    }

}