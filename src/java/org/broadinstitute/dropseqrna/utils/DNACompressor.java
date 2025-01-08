package org.broadinstitute.dropseqrna.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * This compresses and decompresses DNA sequences into a byte array.
 * Each sequence is compressed into a byte array, where each base is represented by 2 bits.
 * The first byte of the compressed data contains the remainder of the sequence length divided by 4.
 * The remaining bytes contain the bases packed into 2 bits each.
 * The decompressed sequence is reconstructed by extracting the bases from the packed bytes.
 *
 * The key difference between this version and DNACompressorVaryingLengths is that this version
 * does not store the lengths of the original sequences in the compressed data, which can save space if the lengths are short.
 *
 */
public class DNACompressor {

    /**
     * Compresses a list of DNA strings into a single byte array, including a leading bitmask for null values.
     *
     * @param dnaList A list of DNA strings to compress. Each string must contain only 'A', 'C', 'G', 'T' or be null.
     * @return A byte array representing the compressed DNA sequences with a leading bitmask.
     */
    public static byte[] compressList(List<String> dnaList) {
        int numEntries = dnaList.size();
        int bitmaskSize = (numEntries + 7) / 8; // 1 bit per entry, rounded up to the nearest byte
        byte[] bitmask = new byte[bitmaskSize];

        List<byte[]> compressedSequences = new ArrayList<>();
        int totalCompressedSize = bitmaskSize; // Include the bitmask size

        // Create the bitmask and compress non-null sequences
        for (int i = 0; i < numEntries; i++) {
            if (dnaList.get(i) == null) {
                // Mark null entries in the bitmask
                bitmask[i / 8] |= (1 << (i % 8));
            } else {
                // Compress non-null sequences
                byte[] compressed = compress(dnaList.get(i));
                compressedSequences.add(compressed);
                totalCompressedSize += compressed.length;
            }
        }

        // Pack the bitmask and compressed data into a single byte array
        byte[] result = new byte[totalCompressedSize];
        System.arraycopy(bitmask, 0, result, 0, bitmaskSize);

        int offset = bitmaskSize;
        for (byte[] compressed : compressedSequences) {
            System.arraycopy(compressed, 0, result, offset, compressed.length);
            offset += compressed.length;
        }

        return result;
    }

    /**
     * Decompresses a single byte array into a list of DNA strings, including handling null values using the bitmask.
     *
     * @param compressed The byte array containing the bitmask and compressed DNA sequences.
     * @param lengths    An array of lengths of the original DNA sequences.
     * @return A list of decompressed DNA strings, including null values.
     */
    public static List<String> decompressList(byte[] compressed, int[] lengths) {
        int numEntries = lengths.length;
        int bitmaskSize = (numEntries + 7) / 8;

        // Extract the bitmask
        byte[] bitmask = new byte[bitmaskSize];
        System.arraycopy(compressed, 0, bitmask, 0, bitmaskSize);

        List<String> decompressed = new ArrayList<>();
        int offset = bitmaskSize;

        // Decompress each sequence based on the bitmask
        for (int i = 0; i < numEntries; i++) {
            if ((bitmask[i / 8] & (1 << (i % 8))) != 0) {
                // Null value indicated by the bitmask
                decompressed.add(null);
            } else {
                int length = lengths[i];
                int remainder = length % 4;
                int byteLength = (length + 3) / 4;

                // Extract the compressed data for this sequence
                byte[] sequenceData = new byte[1 + byteLength];
                System.arraycopy(compressed, offset, sequenceData, 0, 1 + byteLength);

                // Decompress and add to the list
                decompressed.add(decompress(sequenceData));
                offset += (1 + byteLength);
            }
        }

        return decompressed;
    }


    /**
     * Compresses a single DNA string into a byte array.
     */
    public static byte[] compress(String dna) {
        int length = dna.length();
        int remainder = length % 4;
        int byteLength = (length + 3) / 4;

        byte[] compressed = new byte[1 + byteLength];
        compressed[0] = (byte) remainder;

        for (int i = 0; i < length; i++) {
            int base = switch (dna.charAt(i)) {
                case 'A' -> 0b00;
                case 'C' -> 0b01;
                case 'G' -> 0b10;
                case 'T' -> 0b11;
                default -> throw new IllegalArgumentException("Invalid base: " + dna.charAt(i));
            };

            int byteIndex = 1 + (i / 4);
            int bitPosition = (3 - (i % 4)) * 2;
            compressed[byteIndex] |= (base << bitPosition);
        }

        return compressed;
    }

    /**
     * Decompresses a single compressed DNA sequence back into a string.
     */
    public static String decompress(byte[] compressed) {
        int remainder = compressed[0] & 0xFF;
        int byteLength = compressed.length - 1;
        int totalBases = (byteLength * 4) - ((4 - remainder) % 4);

        StringBuilder dna = new StringBuilder(totalBases);

        for (int i = 0; i < totalBases; i++) {
            int byteIndex = 1 + (i / 4);
            int bitPosition = (3 - (i % 4)) * 2;
            int base = (compressed[byteIndex] >> bitPosition) & 0b11;

            char nucleotide = switch (base) {
                case 0b00 -> 'A';
                case 0b01 -> 'C';
                case 0b10 -> 'G';
                case 0b11 -> 'T';
                default -> throw new IllegalStateException("Invalid base bits: " + base);
            };
            dna.append(nucleotide);
        }

        return dna.toString();
    }

    /**
     * Compares two compressed DNA outputs (byte arrays) lexicographically,
     * skipping the bitmask portion of the array.
     *
     * @param lengths The lengths of the entries in the compressed data.
     * @return Comparator for comparing compressed DNA outputs.
     */
    public static Comparator<byte[]> getDecompressedComparator(int[] lengths) {
        return (a, b) -> {
            // Decompress both byte arrays into lists of strings
            List<String> listA = decompressList(a, lengths);
            List<String> listB = decompressList(b, lengths);

            // Compare the decompressed lists lexicographically
            int len = Math.min(listA.size(), listB.size());
            for (int i = 0; i < len; i++) {
                String valueA = listA.get(i);
                String valueB = listB.get(i);

                // Null values are treated as smaller than non-null values
                if (valueA == null && valueB == null) {
                    continue; // Equal
                } else if (valueA == null) {
                    return -1; // Null < Non-null
                } else if (valueB == null) {
                    return 1; // Non-null > Null
                }

                // Compare non-null strings
                int cmp = valueA.compareTo(valueB);
                if (cmp != 0) {
                    return cmp; // Return as soon as a difference is found
                }
            }

            // If all compared strings are equal, the shorter list is smaller
            return Integer.compare(listA.size(), listB.size());
        };
    }


    public static void main(String[] args) {
        // Example DNA sequences

        List<String> dnaList = List.of("ACGT", "TTACG", "GCGTACGTAC");
        System.out.println("Original DNA List: " + dnaList);

        // Compress the DNA sequences
        byte[] compressed = compressList(dnaList);
        System.out.println("Compressed Data: " + java.util.Arrays.toString(compressed));

        // Lengths of the original sequences
        int[] lengths = dnaList.stream().mapToInt(String::length).toArray();

        // Decompress back into the original DNA sequences
        List<String> decompressed = decompressList(compressed, lengths);
        System.out.println("Decompressed DNA List: " + decompressed);

        // Verify that the decompressed sequences match the originals
        System.out.println("Matches original: " + dnaList.equals(decompressed));

        // additional tests with NULL values
        List<String> dnaList2 = Arrays.asList("ACGT", "TTACG", null, "GCGTACGTAC");
        System.out.println("Original DNA List: " + dnaList2);

        // Compress the DNA sequences
        byte[] compressed2 = compressList(dnaList2);
        System.out.println("Compressed Data: " + java.util.Arrays.toString(compressed2));

        // Lengths of the original sequences are the same as the original.  Null values have length 0.
        lengths = dnaList2.stream().mapToInt(dna -> dna == null ? 0 : dna.length()).toArray();

        // Decompress back into the original DNA sequences
        List<String> decompressed2 = decompressList(compressed2, lengths);
        System.out.println("Decompressed DNA List: " + decompressed2);

        // Verify that the decompressed sequences match the originals
        System.out.println("Matches original: " + dnaList2.equals(decompressed2));


    }
}
