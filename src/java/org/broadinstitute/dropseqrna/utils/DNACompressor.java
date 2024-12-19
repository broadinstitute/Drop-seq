package org.broadinstitute.dropseqrna.utils;

import java.util.ArrayList;
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
 */
public class DNACompressor {

    /**
     * Compresses a list of DNA strings into a single byte array.
     *
     * @param dnaList A list of DNA strings to compress. Each string must contain only 'A', 'C', 'G', 'T'.
     * @return A byte array representing the compressed DNA sequences.
     */
    public static byte[] compressList(List<String> dnaList) {
        List<byte[]> compressedSequences = new ArrayList<>();
        int totalSize = 0;

        // Compress each sequence and calculate total size
        for (String dna : dnaList) {
            byte[] compressed = compress(dna);

            /*
            String test = decompress(compressed);
            if (!test.equals(dna)) {
                throw new IllegalStateException("Decompressed DNA does not match original: " + dna + " -> " + test);
            }
            */
            compressedSequences.add(compressed);
            totalSize += compressed.length;
        }

        // Pack all compressed data into a single byte array
        byte[] result = new byte[totalSize];
        int offset = 0;
        for (byte[] compressed : compressedSequences) {
            System.arraycopy(compressed, 0, result, offset, compressed.length);
            offset += compressed.length;
        }

        return result;
    }

    /**
     * Decompresses a single byte array into a list of DNA strings, given the lengths of the original sequences.
     *
     * @param compressed The byte array containing compressed DNA sequences.
     * @param lengths    An array of lengths of the original DNA sequences.
     * @return A list of decompressed DNA strings.
     */
    public static List<String> decompressList(byte[] compressed, int[] lengths) {
        List<String> decompressed = new ArrayList<>();
        int offset = 0;

        // Decompress each sequence using the provided lengths
        for (int length : lengths) {
            int remainder = length % 4;
            int byteLength = (length + 3) / 4;

            // Extract the compressed data for this sequence
            byte[] sequenceData = new byte[1 + byteLength];
            System.arraycopy(compressed, offset, sequenceData, 0, 1 + byteLength);

            // Decompress and add to the list
            decompressed.add(decompress(sequenceData));
            offset += (1 + byteLength);
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

        // Lengths of the original sequences
        int[] lengths = dnaList.stream().mapToInt(String::length).toArray();

        // Decompress back into the original DNA sequences
        List<String> decompressed = decompressList(compressed, lengths);
        System.out.println("Decompressed DNA List: " + decompressed);

        // Verify that the decompressed sequences match the originals
        System.out.println("Matches original: " + dnaList.equals(decompressed));
    }
}
