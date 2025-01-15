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

 * The key difference between this version and DNACompressorVaryingLengths is that this version
 * does not store the lengths of the original sequences in the compressed data, which can save space if the lengths are short.
 *
 */
public class DNACompressor {

    private static final int BITS_PER_BASE = 3;
    private static final int BASES_PER_BYTE = 8 / BITS_PER_BASE;

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
                bitmask[i / 8] |= (byte) (1 << (i % 8));
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
                // Calculate byte length and remainder for the sequence
                int length = lengths[i];
                int[] calc = calculateByteLengthAndRemainder(length);
                int byteLength = calc[0];

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
     * Compresses a single DNA string into a byte array, supporting 'A', 'C', 'G', 'T', and 'N'.
     */
    public static byte[] compress(String dna) {
        int[] calc = calculateByteLengthAndRemainder(dna.length());
        int byteLength = calc[0];
        int remainder = calc[1];

        byte[] compressed = new byte[1 + byteLength];
        compressed[0] = (byte) remainder;

        for (int i = 0; i < dna.length(); i++) {
            int base = switch (dna.charAt(i)) {
                case 'A' -> 0b000;
                case 'C' -> 0b001;
                case 'G' -> 0b010;
                case 'T' -> 0b011;
                case 'N' -> 0b100;
                default -> throw new IllegalArgumentException("Invalid base: " + dna.charAt(i));
            };

            int byteIndex = 1 + (i / BASES_PER_BYTE);
            int bitPosition = (BASES_PER_BYTE - 1 - (i % BASES_PER_BYTE)) * BITS_PER_BASE;
            compressed[byteIndex] |= (byte) (base << bitPosition);
        }

        return compressed;
    }

    /**
     * Decompresses a single compressed DNA sequence back into a string, supporting 'A', 'C', 'G', 'T', and 'N'.
     */
    public static String decompress(byte[] compressed) {
        int remainder = compressed[0] & 0xFF;
        int totalBases = (compressed.length - 1) * BASES_PER_BYTE - ((BASES_PER_BYTE - remainder) % BASES_PER_BYTE);

        StringBuilder dna = new StringBuilder(totalBases);

        for (int i = 0; i < totalBases; i++) {
            int byteIndex = 1 + (i / BASES_PER_BYTE);
            int bitPosition = (BASES_PER_BYTE - 1 - (i % BASES_PER_BYTE)) * BITS_PER_BASE;
            int base = (compressed[byteIndex] >> bitPosition) & 0b111;

            char nucleotide = switch (base) {
                case 0b000 -> 'A';
                case 0b001 -> 'C';
                case 0b010 -> 'G';
                case 0b011 -> 'T';
                case 0b100 -> 'N';
                default -> throw new IllegalStateException("Invalid base bits: " + base);
            };
            dna.append(nucleotide);
        }

        return dna.toString();
    }

    /**
     * Calculates the byte length and remainder for a given DNA sequence length.
     *
     * @param sequenceLength The length of the DNA sequence.
     * @return An array where [0] is the byte length and [1] is the remainder.
     */
    private static int[] calculateByteLengthAndRemainder(int sequenceLength) {
        int remainder = sequenceLength % BASES_PER_BYTE;
        int byteLength = (sequenceLength + BASES_PER_BYTE - 1) / BASES_PER_BYTE;
        return new int[]{byteLength, remainder};
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
}
