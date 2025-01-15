package org.broadinstitute.dropseqrna.utils;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import static org.testng.Assert.*;

public class DNACompressorTest {

    @Test(dataProvider = "dnaStringLists")
    public void testCompressListRoundTrip(String[] dnaList) {
        // Compress the DNA sequences
        byte[] compressed = DNACompressor.compressList(Arrays.asList(dnaList));
        int[] lengths = getListLengths(dnaList);

        // Decompress the DNA sequences
        List<String> result = DNACompressor.decompressList(compressed, lengths);
        // Verify that the decompressed sequences match the original
        assertEquals(result, Arrays.asList(dnaList));
    }

    @Test(dataProvider = "dnaStrings")
    public void testCompressRoundTrip(String dna) {
        // Compress the DNA sequence
        byte[] compressed = DNACompressor.compress(dna);
        // Decompress the DNA sequence
        String result = DNACompressor.decompress(compressed);
        // Verify that the decompressed sequence matches the original
        assertEquals(result, dna);
    }

    @Test
    public void testGetDecompressedComparator() {
        // Sample DNA sequences
        List<String> dnaList1 = Arrays.asList("ACGT", "TTACG", "GCGTACGTAC");
        List<String> dnaList2 = Arrays.asList("ACGT", "TTACG", "GCGTACGTAC");
        List<String> dnaList3 = Arrays.asList("ACGT", "TTACG", "GCGTANGTAC");

        // Compress the DNA sequences
        byte[] compressed1 = DNACompressor.compressList(dnaList1);
        byte[] compressed2 = DNACompressor.compressList(dnaList2);
        byte[] compressed3 = DNACompressor.compressList(dnaList3);

        // Lengths of the original sequences
        int[] lengths = dnaList1.stream().mapToInt(String::length).toArray();

        // Get the comparator
        Comparator<byte[]> comparator = DNACompressor.getDecompressedComparator(lengths);

        // Compare the compressed sequences
        assertEquals(comparator.compare(compressed1, compressed2), 0);
        assertTrue(comparator.compare(compressed1, compressed3) < 0);
        assertTrue(comparator.compare(compressed3, compressed1) > 0);
    }



    @DataProvider(name = "dnaStrings")
    public Object[][] dnaStrings() {
        return new Object[][] {
                {"ACGT"},
                {"TTACG"},
                {"GCGTACGTAC"},
                {"GCGTANGTAC"},
        };
    }

    @DataProvider(name = "dnaStringLists")
    public Object[][] dnaStringLists() {
        return new Object[][] {
                {new String[] {"ACGT", "TTACG", "GCGTACGTAC"}},
                {new String[] {"ACGT", "TTACG", null, "GCGTACGTAC"}},
                {new String[] {"ACGT", "TTACG", null, "GCGTANGTAC"}}
        };
    }

    private int [] getListLengths(String [] dnaList) {
        return Arrays.stream(dnaList).mapToInt(dna -> dna == null ? 0 : dna.length()).toArray();
    }


}