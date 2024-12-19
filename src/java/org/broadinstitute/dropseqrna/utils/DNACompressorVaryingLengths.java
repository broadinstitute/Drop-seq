package org.broadinstitute.dropseqrna.utils;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.List;

/**
 * The `DNACompressor` class provides utilities for compressing and decompressing DNA sequences.
 *
 * Compression:
 * - DNA sequences are represented using 2 bits per nucleotide ('A', 'C', 'G', 'T').
 * - Compressed data is stored in a single byte array.
 * - The metadata includes:
 *   - Number of sequences (4 bytes).
 *   - Length of each compressed sequence (4 bytes per sequence).
 *   - Compressed DNA data concatenated sequentially.
 *
 * Decompression:
 * - Reads the metadata to determine the number and lengths of sequences.
 * - Extracts the compressed DNA data using the length information.
 * - Reconstructs the original DNA sequences from the compressed representation.
 *
 * Use Cases:
 * - Efficient storage and retrieval of DNA sequences of varying lengths.
 * - Eliminates the need for external length information by embedding metadata.
 *
 * If you know your sequences are of fixed length, you can use `DNACompressor` instead.
 * For example if your sequuences were of length 16,16,10,10 this solution would use 34 bytes, and
 * `DNACompressor` would use 14 bytes.
 */

public class DNACompressorVaryingLengths {

    /**
     * Compresses a list of DNA strings into a single byte array.
     *
     * @param dnaList A list of DNA strings to compress. Each string must contain only 'A', 'C', 'G', 'T'.
     * @return A byte array representing the compressed DNA sequences, including metadata for lengths.
     */
    public static byte[] compressList(List<String> dnaList) {
        List<byte[]> compressedSequences = new ArrayList<>();
        int totalSize = 4; // 4 bytes for storing the number of sequences

        // Compress each sequence and calculate total size
        for (String dna : dnaList) {
            byte[] compressed = compress(dna);
            compressedSequences.add(compressed);
            totalSize += 4 + compressed.length; // 4 bytes for length + compressed data
        }

        // Pack all compressed data into a single byte array
        ByteBuffer buffer = ByteBuffer.allocate(totalSize);
        buffer.putInt(dnaList.size()); // Number of sequences

        for (byte[] compressed : compressedSequences) {
            buffer.putInt(compressed.length); // Length of each compressed sequence
            buffer.put(compressed); // Compressed data
        }

        return buffer.array();
    }

    /**
     * Decompresses a single byte array into a list of DNA strings.
     *
     * @param compressed The byte array containing compressed DNA sequences with metadata.
     * @return A list of decompressed DNA strings.
     */
    public static List<String> decompressList(byte[] compressed) {
        ByteBuffer buffer = ByteBuffer.wrap(compressed);

        // Read the number of sequences
        int numSequences = buffer.getInt();
        List<String> decompressed = new ArrayList<>();

        // Decompress each sequence
        for (int i = 0; i < numSequences; i++) {
            int length = buffer.getInt(); // Length of the compressed sequence
            byte[] sequenceData = new byte[length];
            buffer.get(sequenceData); // Extract the compressed sequence
            decompressed.add(decompress(sequenceData));
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

    public static void main(String[] args) {
        // Example DNA sequences
        List<String> dnaList = List.of("ACGT", "TTACG", "GCGTACGTAC");
        System.out.println("Original DNA List: " + dnaList);

        // Compress the DNA sequences
        byte[] compressed = compressList(dnaList);
        System.out.println("Compressed Data: " + java.util.Arrays.toString(compressed));

        // Decompress back into the original DNA sequences
        List<String> decompressed = decompressList(compressed);
        System.out.println("Decompressed DNA List: " + decompressed);

        // Verify that the decompressed sequences match the originals
        System.out.println("Matches original: " + dnaList.equals(decompressed));
    }
}
